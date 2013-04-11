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


plotGate(gs_bcell, 2, lattice = TRUE, xbin = 128,
                  cond = "factor(Center):factor(Sample):factor(Replicate)")


# TODO: This should not have to be done manually. OpenCyto should do this.
# Related: https://github.com/RGLab/openCyto/issues/30
#  and https://github.com/RGLab/openCyto/issues/31
# Manually applies more accurate IgD gate that OpenCyto does not currently
# support.
Rm("IgD-cd27+", gs_bcell)
Rm("IgD+cd27+", gs_bcell)
Rm("IgD+cd27-", gs_bcell)
Rm("IgD-cd27-", gs_bcell)
gating.bsub(gs_bcell, parent = "cd19 & cd20")


archive(gs_bcell, file = gs_file)
