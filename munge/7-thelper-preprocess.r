library(ProjectTemplate)
load.project()

panel <- "Th1/2/17"
path_Lyoplate <- "/loc/no-backup/ramey/Lyoplate"
plot <- FALSE

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
                         "CD8", "CD3", "HLADR")

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

# To better estimate transformation parameters, we first apply three gates:
# 1. Boundary on FSC-A and SSC-A
# 2. Debris
# 3. Lymphocytes
gs_pregate <- GatingSet(flow_set)

# Creates the gating-template object from a CSV file
gt_csv <- "gt-preprocess.csv"
gating_template <- gatingTemplate(gt_csv, panel)

# Boundary gate
boundary_gate <- rectangleGate(filterId = "boundary", "FSC-A" = c(0, 2.5e5),
                               "SSC-A" = c(0, 2.5e5))
add(gs_pregate, boundary_gate)
recompute(gs_pregate)

# Applies OpenCyto to GatingSet
gating(gating_template, gs_pregate, mc.cores = 8, parallel_type = "multicore",
       prior_group = "Center")

# Next, we extract the flowSet the lymphocyte subpopulation
fs_lymph <- getData(gs_pregate, "lymph")

# For each center, we compute the 5th percentile of each marker's negative values
# These percentiles determine the center-specific transformation for each marker
fs_split <- split(fs_lymph, pData(fs_lymph)$Center)
widths_quantiles <- lapply(fs_split, function(fs) {
  quantiles <- fsApply(fs, each_col, quantile_negatives, probs = 0.05)
  # Removes FSC-A and SSC-A
  quantiles <- quantiles[, -(1:2)]
  quantiles <- apply(quantiles, 2, median, na.rm = TRUE)
  widths <- quantile2width(quantiles)

  list(quantiles = quantiles, widths = widths)
})
quantiles <- do.call(rbind, lapply(widths_quantiles, "[[", "quantiles"))
quantiles <- melt(quantiles)
colnames(quantiles) <- c("Center", "Marker", "Quantile")

widths <- do.call(rbind, lapply(widths_quantiles, "[[", "widths"))
widths <- melt(widths)
colnames(widths) <- c("Center", "Marker", "Width")

if (plot) {
  p <- ggplot(quantiles, aes(Marker, Quantile)) + geom_boxplot()
  p + ylab("Center-specific Quantiles") + ggtitle("Median Negative Quantiles -- Lymphocytes")
  p <- ggplot(widths, aes(Marker, Width)) + geom_boxplot()
  p + ylab("Center-specific Widths") + ggtitle("Median Widths -- Lymphocytes")
}

# Next, we apply the center-specific transformations to each marker.
message("Applying transformations for each center...")
for (center in names(fs_split)) {
  message("Center: ", center)
  widths_center <- subset(widths, Center == center)
  widths_center$Marker <- as.character(widths_center$Marker)

  for (i in seq_len(nrow(widths_center))) {
    trans_marker <- with(widths_center, transformList(Marker[i], FCSTransTransform(w = Width[i])))
    fs_split[[center]] <- transform(fs_split[[center]], trans_marker)
  }
}

# Removes median widths from 'widths' to avoid confusion when transformations are saved.
# widths <- subset(widths, select = -median_width)

# Saves transformations to a CSV file
write.csv(widths, file = "data/thelper-transformations.csv", row.names = FALSE)

# Merges the list of flowSet objects into a single flowSet object. This code is
# verbose but it circumvents an issue introduced recently in flowIncubator.
flow_set <- fs_split[[1]]
for (i in seq_along(fs_split)[-1]) {
  flow_set <- rbind2(flow_set, fs_split[[i]])
}

# Creates a new GatingSet without the upstream gates. We copy those from the
# previous gating set to the new GatingSet, which has the transformed data.
gs_thelper <- GatingSet(flow_set)

boundary_gates <- getGate(gs_pregate, "boundary")
add(gs_thelper, boundary_gates)
recompute(gs_thelper)

nonDebris_gates <- getGate(gs_pregate, "nonDebris")
add(gs_thelper, nonDebris_gates, parent = "boundary")
recompute(gs_thelper)

lymph_gates <- getGate(gs_pregate, "lymph")
add(gs_thelper, lymph_gates, parent = "nonDebris")
recompute(gs_thelper)

# Archives the results
save_gs(gs_thelper, path = file.path(path_Lyoplate, "gating-sets/gs-thelper"))

