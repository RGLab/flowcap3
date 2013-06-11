library(ProjectTemplate)
load.project()

panel <- "B cell"

path_Lyoplate <- "/loc/no-backup/ramey/Lyoplate"

# There is a bug in the built-in 'list.dirs' function. The argument 'full.names'
# does not work as advertised. After a quick Google search, others recently have
# had similar results.
# To fix this, I manually grab the directory names using 'strsplit'.
centers <- list.dirs(path_Lyoplate, recursive = FALSE, full.names = FALSE)
centers <- sapply(strsplit(centers, split = "/"), tail, n = 1)
centers <- setdiff(centers, "gating-sets")

# Because BSMS apparently did not include IgD in their FCS files, we remove them
# from consideration for now.
centers <- setdiff(centers, "BSMS")

# These are the markers that we will keep after the data have been preprocessed.
markers_of_interest <- c("FSC-A", "SSC-A", "Live", "CD3", "CD19", "CD20", "IgD",
                         "CD27", "CD38", "CD24")

# For each center, we construct a flowSet of FCS files after compensating and
# transforming the flowSet created from the FCS files in the center's
# subdirectory under 'path_Lyoplate'.
lyoplate_list <- lapply(centers, function(center) {
  message("Center: ", center)
  center <- "NHLBI"
  path <- file.path(path_Lyoplate, center)

  # The filename of the manually edited Excel file: assumed to be in getwd()
  xlsx <- dir(pattern = center)

  # Constructs a flowSet object for the current center. The flowSet will be
  # compensated and transformed.
  lyoplate_out <- flowset_lyoplate(path = path, xlsx = xlsx,
                                   comp_matrix = compensation_matrices[[center]],
                                   transform = "FCSTrans", center = center, panel = panel)

  # Swaps the channels and markers for the current 'flowSet' object. This ensures
  # that we can 'rbind2' the 'GatingSetList' below because the stain names do not
  # match otherwise.
  lyoplate_out$flow_set <- fsApply(lyoplate_out$flow_set, preprocess_flowframe,
                                   markers_keep = markers_of_interest)

  w <- melt(lyoplate_out$w)
  colnames(w) <- c("w", "Channel")
  w$Marker <- openCyto:::channels2markers(lyoplate_out$flow_set[[1]], channels = w$Channel)

  lyoplate_out$w <- w

  lyoplate_out
})
names(lyoplate_list) <- centers

# Maps channels to markers
lyoplate_list <- lapply(centers, function(center) {
  channels <- lyoplate_list[[center]]$w$Channel

  markers <- sapply(channels, function(channel) {
    openCyto:::getChannelMarker(lyoplate_list[[center]]$flow_set[[1]], channel)$name
  })

  lyoplate_list[[center]]$w$Marker <- markers
  
  lyoplate_list[[center]]
})
names(lyoplate_list) <- centers

# Merges the list of flowSet objects into a single flowSet object. This code is
# verbose but it circumvents an issue introduced recently in flowIncubator.
flow_set <- lyoplate_list[[1]]$flow_set
for (i in seq.int(2, length(lyoplate_list))) {
  flow_set <- rbind2(flow_set, lyoplate_list[[i]]$flow_set)
}

gs_bcell <- GatingSet(flow_set)

# Saves estimated values of 'w' for FCSTrans
w_FCSTrans <- lapply(centers, function(center) {
  cbind(Center = center, lyoplate_list[[center]]$w)
})
w_FCSTrans <- do.call(rbind, w_FCSTrans)
save(w_FCSTrans, file = "FCSTrans_transformation_estimates.RData")

# Archives the results
save_gs(gs_bcell, path = file.path(path_Lyoplate, "gating-sets/gs-bcell"))

# Plots FCSTrans estimates
p <- ggplot(w_FCSTrans, aes(x = Center, weight = w)) + geom_bar() + facet_wrap(~ Marker)
p + ylab("Estimate of FCSTrans Width, w")
