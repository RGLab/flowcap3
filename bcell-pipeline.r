library(ProjectTemplate)
load.project()
set.seed(42)

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
gating(gating_template, gs_bcell, mc.cores = 8, parallel_type = "multicore")

# Archives the GatingSet
save_gs(gs_bcell, path = gs_path, overwrite = TRUE)

plotGate(gs_bcell, "lymph", lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)")

png("bcell-gates/transitional.png", width = 960, height = 960)
print(plotGate(gs_bcell, "transitional", lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)"))
dev.off()

png("bcell-gates/IgD.png", width = 960, height = 960)
plotGate(gs_bcell, 13:16, lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)")
dev.off()

png("bcell-gates/plasmablast.png", width = 960, height = 960)
plotGate(gs_bcell, "plasmablast", lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)")
dev.off()

