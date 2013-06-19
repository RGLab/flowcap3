library(ProjectTemplate)
load.project()

panel <- "Bcell"
gs_path <- "/loc/no-backup/ramey/Lyoplate/gating-sets/gs-bcell"

# Loads the archived gatingSet object
gs_bcell <- load_gs(gs_path)

# Creates the gating-template object from a CSV file
gt_csv <- "gt-bcell.csv"
gating_template <- gatingTemplate(gt_csv, panel)

# TEMP
Rm("boundary", gs_bcell)
Rm("nonDebris", gs_bcell)
Rm("lymph", gs_bcell)

# Applies OpenCyto to GatingSet
gating(gating_template, gs_bcell, mc.cores = 8, parallel_type = "multicore", prior_group = "Center")
# gating(gating_template, gs_bcell, prior_group = "Center")

# Archives the GatingSet
save_gs(gs_bcell, path = gs_path, overwrite = TRUE)

plotGate(gs_bcell, "transitional", lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)")

