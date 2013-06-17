library(ProjectTemplate)
load.project()

panel <- "Th1/2/17"

path_Lyoplate <- "/loc/no-backup/ramey/Lyoplate"

# There is a bug in the built-in 'list.dirs' function. The argument 'full.names'
# does not work as advertised. After a quick Google search, others recently have
# had similar results.
# To fix this, I manually grab the directory names using 'strsplit'.
centers <- list.dirs(path_Lyoplate, recursive = FALSE, full.names = FALSE)
centers <- sapply(strsplit(centers, split = "/"), tail, n = 1)
centers <- setdiff(centers, "gating-sets")

# Because BSMS apparently did not include CXCR3 in their FCS files, we
# remove them from consideration for now.
centers <- setdiff(centers, "BSMS")

# These are the markers that we will keep after the data have been preprocessed.
markers_of_interest <- c("FSC-A", "SSC-A", "Live", "CXCR3", "CD4", "CCR6", "CD38",
                         "CD8", "CD3", "HLA-DR")

# For each center, we construct a flowSet of FCS files after compensating and
# transforming the flowSet created from the FCS files in the center's
# subdirectory under 'path_Lyoplate'.
message("Compensating and Estimating Transformation Parameters")
lyoplate_list <- lapply(centers, function(center) {
  message("Center: ", center)
  path <- file.path(path_Lyoplate, center)

  # The filename of the manually edited Excel file: assumed to be in getwd()
  xlsx <- dir(pattern = center)

  # Constructs a compensated flowSet object for the current center
  lyoplate_out <- flowset_lyoplate(path = path, xlsx = xlsx,
                                   comp_matrix = compensation_matrices[[center]],
                                   center = center, panel = panel)

  # Maps channels to markers
  lyoplate_out$widths$Marker <- sapply(lyoplate_out$widths$Channel, function(channel) {
    openCyto:::getChannelMarker(lyoplate_out$flow_set[[1]], channel)$desc
  })

  lyoplate_out
})
names(lyoplate_list) <- centers
fs_list <- lapply(lyoplate_list, "[[", "flow_set")

# Saves estimated widths ('w') for FCSTrans
trans_widths <- lapply(lyoplate_list, "[[", "widths")
trans_widths <- do.call(rbind, trans_widths)

# Converts the marker names to a common name
trans_widths$Marker <- marker_conversion(trans_widths$Marker)

# Calculates the median of the widths across centers for each marker
widths_summary <- with(trans_widths, aggregate(width, by = list(Marker), median))
colnames(widths_summary) <- c("Marker", "median_width")
trans_widths <- plyr:::join(trans_widths, widths_summary)
trans_widths$width_applied <- trans_widths$median_width

# Updates data.frame to include only one width...not 3
# This would be confusing otherwise
trans_widths <- subset(trans_widths, select = c(Center, Marker, Channel, width_applied))
colnames(trans_widths)[4] <- "Width"

# Saves transformations to a CSV file
write.csv(trans_widths, file = "data/thelper-transformations.csv", row.names = FALSE)

# Applies the FCStrans transformation to each center
message ("Transforming flowSets for each center")
dev_null <- with(trans_widths, mapply(function(center, channel, width) {
  message("Center: ", center, " -- Channel: ", channel)
  trans_channel <- transformList(from = channel,
                                 tfun = FCSTransTransform(w = width))
  fs_list[[center]] <<- transform(fs_list[[center]], trans_channel)

  NULL
}, Center, Channel, Width))


# Swaps the channels and markers for the current 'flowSet' object. This ensures
# that we can 'rbind2' the 'GatingSetList' below because the stain names do not
# match otherwise.
message ("Swapping flowSet channels and markers")
fs_list <- lapply(centers, function(center) {
  message("Center: ", center)

  fsApply(fs_list[[center]], preprocess_flowframe,
          markers_keep = markers_of_interest)
})
names(fs_list) <- centers

# Merges the list of flowSet objects into a single flowSet object. This code is
# verbose but it circumvents an issue introduced recently in flowIncubator.
flow_set <- fs_list[[1]]
for (i in seq.int(2, length(lyoplate_list))) {
  flow_set <- rbind2(flow_set, fs_list[[i]])
}

gs_thelper <- GatingSet(flow_set)

# Archives the results
save_gs(gs_thelper, path = file.path(path_Lyoplate, "gating-sets/gs-thelper"))
