library(ProjectTemplate)
load.project()

panel <- "B cell"
path_Lyoplate <- "/loc/no-backup/ramey/Lyoplate"
plot <- FALSE

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
gs_bcell <- GatingSet(flow_set)

# Creates the gating-template object from a CSV file
gt_csv <- "gt-preprocess.csv"
gating_template <- gatingTemplate(gt_csv, panel)

# Boundary gate
boundary_gate <- rectangleGate(filterId = "boundary", "FSC-A" = c(0, 2.5e5),
                               "SSC-A" = c(0, 2.5e5))
add(gs_bcell, boundary_gate)
recompute(gs_bcell)

# Applies OpenCyto to GatingSet
gating(gating_template, gs_bcell, mc.cores = 8, parallel_type = "multicore",
       prior_group = "Center")

# Next, we extract the flowSet the lymphocyte subpopulation
fs_lymph <- getData(gs_bcell, "lymph")

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

# For the following centers/markers, we have discovered that aggregation of the
# widths across all centers improves the transformations of the markers listed.
# Baylor - CD27
# NHLBI - CD27
# UCLA - CD19, CD27
# Yale - CD19, CD27, Live
widths <- ddply(widths, .(Marker), transform, median_width = median(Width))
widths$Marker <- as.character(widths$Marker)
widths$Center <- as.character(widths$Center)

improve_indices <- with(widths,
                        (Marker == "CD27" & Center %in% c("Baylor", "NHLBI", "UCLA", "Yale")) |
                        (Marker == "CD19" & Center %in% c("UCLA", "Yale")) |
                        (Marker == "Live" & Center == "Yale"))
widths <- within(widths, Width[improve_indices] <- median_width[improve_indices])

# TODO: Problem markers that persist. Need to update these manually
# NHLBI - CD27 (Subjective value previously was 0.8)
# UCLA - CD27 (Much improved, but sample 4 has lots of zeros. Subjective value previously was 1.6)
# UCLA - CD3 (Double-check this)
# Yale - CD3 (Double-check this) - May be fixed by improved compensation
# Yale - CD19 (Double-check this) - May be fixed by improved compensation
# Yale - CD27 - May be fixed by improved compensation

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
widths <- subset(widths, select = -median_width)

# Saves transformations to a CSV file
write.csv(widths, file = "data/bcell-transformations.csv", row.names = FALSE)

# Merges the list of flowSet objects into a single flowSet object. This code is
# verbose but it circumvents an issue introduced recently in flowIncubator.
flow_set <- fs_split[[1]]
for (i in seq_along(fs_split)[-1]) {
  flow_set <- rbind2(flow_set, fs_split[[i]])
}

# Creates a new GatingSet
gs_bcell <- GatingSet(flow_set)

# TODO:
# By creating a second GatingSet, we have removed the gates previously applied.
# This conflicts with the GatingTemplate CSV file. These gates will be reapplied
# and likely cause an issue. With this in mind, I have updated the CSV files
# and removed the 3 gates. The CD3 gate's parent is 'root'. This needs to be
# fixed.


# Archives the results
save_gs(gs_bcell, path = file.path(path_Lyoplate, "gating-sets/gs-bcell"))

