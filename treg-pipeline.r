library(ProjectTemplate)
load.project()
set.seed(42)

panel <- "Treg"
gs_path <- "/loc/no-backup/ramey/Lyoplate/gating-sets/gs-treg"

# Loads the archived gatingSet object
gs_treg <- load_gs(gs_path)

# Creates the gating-template object from a CSV file
gt_csv <- "gt-treg.csv"
gating_template <- gatingTemplate(gt_csv, panel)

# Applies OpenCyto to GatingSet
gating(gating_template, gs_treg, mc.cores = 10, parallel_type = "multicore") #, prior_group = "Center")

png("treg-gates/CCR4_CD45RO.png", width = 1200, height = 1200)
print(plotGate(gs_treg, 14:17, lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)"))
dev.off()

png("treg-gates/CCR4_HLADR.png", width = 1200, height = 1200)
print(plotGate(gs_treg, 18:21, lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)"))
dev.off()

system("scp -r treg-gates/ 10.6.156.144:/Users/jramey/Dropbox/rglab/flowcap3/")

# Archives the GatingSet
save_gs(gs_treg, path = gs_path, overwrite = TRUE)

