library(ProjectTemplate)
load.project()

panel <- "Bcell"
gs_path <- "/shared/silo_researcher/Gottardo_R/ramey_working/Lyoplate/gs-bcell"

# Loads the archived gatingSet object
library(flowIncubator)
gs_bcell <- load_gs(gs_path)

# Creates the gating-template object from a CSV file
gt_csv <- "gt-bcell.csv"
gating_template <- gatingTemplate(gt_csv, panel)

# Applies OpenCyto to GatingSet
gating(gating_template, gs_bcell, num_nodes = 12, parallel_type = "multicore")

# Archives the GatingSet
save_gs(gs_bcell, path = gs_path, overwrite = TRUE)

plotGate(gs_bcell, 4, lattice = TRUE, xbin = 128,
                  cond = "factor(Center):factor(Sample):factor(Replicate)")



