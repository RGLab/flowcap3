library(ProjectTemplate)
load.project()
set.seed(42)

panel <- "Tcell"
gs_path <- "/loc/no-backup/ramey/Lyoplate/gating-sets/gs-tcell"

# Loads the archived gatingSet object
gs_tcell <- load_gs(gs_path)

# Creates the gating-template object from a CSV file
gt_csv <- "gt-tcell.csv"
gating_template <- gatingTemplate(gt_csv, panel)

# Applies OpenCyto to GatingSet
gating(gating_template, gs_tcell, mc.cores = 8, parallel_type = "multicore") #, prior_group = "Center")

# gating(gating_template, gs_tcell)

png("tcell-gates/cd4cd8.png", width = 960, height = 960)
print(plotGate(gs_tcell, 5:8, lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)"))
dev.off()

png("tcell-gates/activated cd4.png", width = 960, height = 960)
plotGate(gs_tcell, "activated cd4", lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)")
dev.off()

png("tcell-gates/activated cd8.png", width = 960, height = 960)
plotGate(gs_tcell, "activated cd8", lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)")
dev.off()

png("tcell-gates/cd8_ccr7_cd45ra.png", width = 960, height = 960)
plotGate(gs_tcell, 18:21, lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)")
dev.off()

png("tcell-gates/cd4_ccr7_cd45ra.png", width = 960, height = 960)
plotGate(gs_tcell, 23:26, lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)")
dev.off()

system("scp -r tcell-gates/ 10.6.156.144:/Users/jramey/Dropbox/rglab/flowcap3/")


# Archives the GatingSet
save_gs(gs_tcell, path = gs_path, overwrite = TRUE)
