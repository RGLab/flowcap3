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
markers_of_interest <- c("FSC-A", "SSC-A", "CXCR3", "CD4", "CCR6", "CD38", "CD8",
                         "CD3", "HLA-DR")

# For each center, we construct a flowSet of FCS files after compensating and
# transforming the flowSet created from the FCS files in the center's
# subdirectory under 'path_Lyoplate'.
fs_list <- lapply(centers, function(center) {
  message("Center: ", center)
  path <- file.path(path_Lyoplate, center)

  # The filename of the manually edited Excel file: assumed to be in getwd()
  xlsx <- dir(pattern = center)

  # Constructs a flowSet object for the current center. The flowSet will be
  # compensated and transformed.
  flow_set <- flowset_lyoplate(path = path, xlsx = xlsx,
                               comp_matrix = compensation_matrices[[center]],
                               transform = TRUE, center = center, panel = panel)

  # Swaps the channels and markers for the current 'flowSet' object. This ensures
  # that we can 'rbind2' the 'GatingSetList' below because the stain names do not
  # match otherwise.
  flow_set <- fsApply(flow_set, preprocess_flowframe,
                      markers_keep = markers_of_interest)

  flow_set
})

# Merges the list of flowSet objects into a single flowSet object. This code is
# verbose but it circumvents an issue introduced recently in flowIncubator.
flow_set <- fs_list[[1]]
for (i in seq.int(2, length(fs_list))) {
  flow_set <- rbind2(flow_set, fs_list[[i]])
}

gs_th1_2_17 <- GatingSet(flow_set)

# Archives the results
save_gs(gs_th1_2_17, path = file.path(path_Lyoplate, "gating-sets/gs-th1-2-17"))

