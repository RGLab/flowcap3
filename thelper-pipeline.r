library(ProjectTemplate)
load.project()

panel <- "Thelper"
gs_path <- "/loc/no-backup/ramey/Lyoplate/gating-sets/gs-thelper"

# Loads the archived gatingSet object
gs_thelper <- load_gs(gs_path)

# Creates the gating-template object from a CSV file
gt_csv <- "gt-thelper.csv"
gating_template <- gatingTemplate(gt_csv, panel)

# Applies OpenCyto to GatingSet
set.seed(42)
gating(gating_template, gs_thelper, mc.cores = 12, parallel_type = "multicore") #prior_group = "Center")

debug(openCyto:::.gating_gtMethod)
gating(gating_template, gs_thelper)

png("thelper-gates/cd4cd8.png", width = 1200, height = 1200)
print(plotGate(gs_thelper, c(6,8), lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)"))
dev.off()

png("thelper-gates/cd4_cd38_HLADR.png", width = 1200, height = 1200)
print(plotGate(gs_thelper, 15:18, lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)"))
dev.off()

png("thelper-gates/cd4_CXCR3_CCR6.png", width = 1200, height = 1200)
print(plotGate(gs_thelper, 17:20, lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)"))
dev.off()


png("thelper-gates/cd8_cd38_HLADR.png", width = 1200, height = 1200)
print(plotGate(gs_thelper, 23:26, lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)"))
dev.off()

# [21] "CXCR3+CCR6+"    "CXCR3-CCR6+"    "CXCR3+CCR6-"    "CXCR3-CCR6-"   
png("thelper-gates/cd8_CXCR3_CCR6.png", width = 1200, height = 1200)
print(plotGate(gs_thelper, 29:32, lattice = TRUE, xbin = 128, margin = TRUE,
         cond = "factor(Center):factor(Sample):factor(Replicate)"))
dev.off()


system("scp -r thelper-gates/ 10.6.156.144:/Users/jramey/Dropbox/rglab/flowcap3/")


# Archives the GatingSet
save_gs(gs_thelper, path = gs_path, overwrite = TRUE)
