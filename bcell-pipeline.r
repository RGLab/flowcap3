library(ProjectTemplate)
load.project()

path_Lyoplate <- "/loc/no-backup/ramey/Lyoplate/"
panel <- "Bcell"

gs_file <- file.path(path_Lyoplate, "gs-Bcell.tar")

# Loads the archived gatingSet object
gs_bcell <- unarchive(file = gs_file, path = path_Lyoplate)

# Creates the gating-template object from a CSV file
# Applies OpenCyto to GatingSet
gt_csv <- "gt-bcell.csv"
gating_template <- gatingTemplate(gt_csv, panel)
gating(gating_template, gs_bcell, num_nodes = 12, parallel_type = "multicore")

archive(gs_bcell, file = gs_file)

plotGate(gs_bcell, 2, lattice = TRUE, xbin = 128,
         cond = "factor(Center):factor(Sample):factor(Replicate)")

# Manually applies more accurate IgD gate that OpenCyto does not currently support.
Rm("IgD-cd27+", gs_bcell)
Rm("IgD+cd27+", gs_bcell)
Rm("IgD+cd27-", gs_bcell)
Rm("IgD-cd27-", gs_bcell)
gating.bsub(gs_bcell, parent = "cd19 & cd20")

# TODO: Add cytokine prior elicitation to OpenCyto with artificial peaks option
#   This would ensure that OpenCyto can deliver the same results as the less automatic IgD gate above.
