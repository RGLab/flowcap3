library(ProjectTemplate)
load.project()

path_Lyoplate <- "/loc/no-backup/ramey/Lyoplate/"
panel <- "Bcell"

centers <- dir(path_Lyoplate, include.dirs = TRUE)

# For each center, we construct a GatingSet of FCS files after compensating and
# transforming the flowSet created from the FCS files in the center's
# subdirectory under 'path_Lyoplate'.
gs_list <- lapply(centers, function(center) {
  message("Center: ", center)
  path <- file.path(path_Lyoplate, center)

  # The filename of the manually edited Excel file: assumed to be in getwd()
  xlsx <- dir(pattern = center)
  comp_matrix <- compensation_lyoplate(path = path, xlsx = xlsx, panel = panel)

  gatingset_lyoplate(path = path, xlsx = xlsx, comp_matrix = comp_matrix,
                     center = center, panel = panel)
})

# Creates a GatingSetList and then rbind2 to form a single GatingSet
gs_bcell <- GatingSetList(gs_list)
gs_bcell <- rbind2(gs_bcell)

# TODO: Fix error when 'rbind2' is called
# Error: "The individual flowFrames do not contain identical stains."

# Creates the gating-template object from a CSV file
# Applies OpenCyto to GatingSet
gt_csv <- "/loc/no-backup/ramey/cytotrol/Bcell/Cytotrol_Bcell_GatingTemplate.csv"
gating_template <- gatingTemplate(gt_csv, "Bcell")
gating(gating_template, gs_bcell, num_nodes = 10, parallel_type = "multicore")

# TODO: Add cytokine prior elicitation to OpenCyto with artificial peaks option
