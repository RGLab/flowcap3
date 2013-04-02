library(ProjectTemplate)
load.project()

gs_file <- "/loc/no-backup/ramey/cytotrol/Bcell/Cytotrol-Bcell.tar"

# Loads the archived gatingSet object
gs_bcell <- unarchive(file = gs_file, path = "/loc/no-backup/ramey/cytotrol/Bcell/")

# Creates the gating-template object from a CSV file
gt_csv <- "/loc/no-backup/ramey/cytotrol/Bcell/Cytotrol_Bcell_GatingTemplate.csv"
gating_template <- gatingTemplate(gt_csv, "Bcell")
gating(gating_template, gs_bcell, num_nodes = 10, parallel_type = "multicore")

archive(gs_bcell, file = gs_file)

plotGate(gs_bcell, 13, lattice = TRUE, xbin = 128, cond = "factor(Center):factor(Replicate)")


Rm("IgD-cd27+", gs_bcell)
Rm("IgD+cd27+", gs_bcell)
Rm("IgD+cd27-", gs_bcell)
Rm("IgD-cd27-", gs_bcell)


gating.bsub(gs_bcell, parent = "cd19 & cd20")
