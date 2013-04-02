library(ProjectTemplate)
load.project()

panel <- "Bcell"
path <- "/loc/no-backup/ramey/cytotrol"

# For now, we manually remove Miami and Yale because they do not contain the FSC-H marker.
# In FlowCAP 3 we used a special singlet gate for both centers,
# but currently OpenCyto cannot handle this use case.
centers <- c("Baylor", "BSMS", "CIMR", "Kings", "Miami", "NHLBI", "Stanford", "UCLA", "Yale")

fs_list <- lapply(centers, function(center) {
  data_path <- file.path(path, panel, center)
  fcs_files <- list.files(data_path, pattern = ".fcs", full.names = TRUE)
  read.flowSet(fcs_files, min.limit = -100)
})

# These are the markers that we will keep after the data have been preprocessed.
markers_of_interest <- c("FSC-A", "SSC-A", "CD3", "CD19", "CD20", "IgD", "CD27",
                         "CD38", "CD24")

fs_list <- lapply(seq_along(centers), function(i) {
  fs_center <- fs_list[[i]]

  # Compensates the flowSet
  comp_mat <- keyword(fs_center[[1]])$SPILL
  fs_center <- compensate(fs_center, comp_mat)

  # Transforms the flowSet
  transformation <- estimateMedianLogicle(fs_center, channels = colnames(comp_mat))
  fs_center <- transform(fs_center, transformation)
  
  flow_set <- fsApply(fs_center, preprocess_flowframe,
                      markers_keep = markers_of_interest)
  flow_set <- update_phenodata(flow_set)

  flow_set
})

gs_list <- lapply(fs_list, GatingSet)
gs_list <- GatingSetList(gs_list)

gs_bcell <- rbind2(gs_list)

archive(gs_bcell, file = "/loc/no-backup/ramey/cytotrol/Bcell/Cytotrol-Bcell.tar")

