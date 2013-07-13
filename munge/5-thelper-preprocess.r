library(ProjectTemplate)
load.project()

panel <- "Thelper"
path_workspaces <- "workspaces/XML"
path_lyoplate <- "/loc/no-backup/ramey/Lyoplate"

# Should we extract the spillover matrices for each center?
# By default, no because the data are already compensated and transformed.
extract_spillover <- FALSE

# These are the markers that we will keep after the data have been preprocessed.
markers_of_interest <- c("FSC-A", "SSC-A", "Live", "CXCR3", "CD4", "CCR6", "CD38",
                         "CD8", "CD3", "HLADR")

# There is a bug in the built-in 'list.dirs' function. The argument 'full.names'
# does not work as advertised. After a quick Google search, others recently have
# had similar results.
# To fix this, I manually grab the directory names using 'strsplit'.
centers <- list.dirs(path_workspaces, recursive = FALSE, full.names = FALSE)
centers <- sapply(strsplit(centers, split = "/"), tail, n = 1)

# If specified, extracts spillover matrices for each center as a list
if (extract_spillover) {
  message("Extracting spillover matrices...")
  spillover_matrices <- sapply(centers, function(center) {
    spillover_file <- file.path(path_workspaces, center, "spillover-matrices", panel)
    read.table(spillover_file, skip = 2, nrows = 8, header = TRUE, check.names = FALSE, sep = "\t")
  }, simplify = FALSE)
}

# We parse each center's workspace and then extract its flowSet.
message("Parsing Lyoplate workspaces...")
fs_list <- sapply(centers, function(center) {
  message("Center: ", center)
  fcs_path <- file.path(path_lyoplate, center)
  ws <- openWorkspace(file.path(path_workspaces, center, "workspace.xml"))
  gating_set <- parseWorkspace(ws, name = panel, path = fcs_path, isNcdf = FALSE)
  closeWorkspace(ws)
  getData(gating_set)
}, simplify = FALSE)

# Swaps the channels and markers for the current 'flowSet' object. This ensures
# that we can 'rbind2' the 'GatingSetList' below because the stain names do not
# match otherwise.
message ("Swapping flowSet channels and markers")
fs_list <- sapply(centers, function(center) {
  message("Center: ", center)

  fsApply(fs_list[[center]], preprocess_flowframe,
          markers_keep = markers_of_interest)
}, simplify = FALSE)

# Merges the list of flowSet objects into a single flowSet object. This code is
# verbose but it circumvents an issue introduced recently in flowIncubator.
flow_set <- fs_list[[1]]
for (i in seq.int(2, length(fs_list))) {
  flow_set <- rbind2(flow_set, fs_list[[i]])
}

# Creates GatingSet from the flowSet
gs_thelper <- GatingSet(flow_set)

# Archives the results
save_gs(gs_thelper, path = file.path(path_lyoplate, "gating-sets/gs-thelper"))

