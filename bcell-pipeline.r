library(ProjectTemplate)
load.project()

panel <- "Bcell"
gs_path <- "/loc/no-backup/ramey/Lyoplate/gating-sets/gs-bcell"

# Loads the archived gatingSet object
gs_bcell <- load_gs(gs_path)

# Creates the gating-template object from a CSV file
gt_csv <- "gt-bcell.csv"
gating_template <- gatingTemplate(gt_csv, panel)

# Temp
Rm("nonDebris", gs_bcell)

# Boundary gate
boundary_gate <- rectangleGate(filterId = "boundary", "FSC-A" = c(0, 2.5e5),
                               "SSC-A" = c(0, 2.5e5))
add(gs_bcell, boundary_gate)
recompute(gs_bcell)

# Applies OpenCyto to GatingSet
gating(gating_template, gs_bcell, mc.cores = 8, parallel_type = "multicore", stop.at = "cd3",
       prior_group = "Center")

# Archives the GatingSet
save_gs(gs_bcell, path = gs_path, overwrite = TRUE)

plotGate(gs_bcell, "lymph", lattice = TRUE, xbin = 128, margin = TRUE,
                  cond = "factor(Center):factor(Sample):factor(Replicate)")



# Temp
fs_lymph <- getData(gs_bcell, "lymph")
center <- "Yale"
fcs_file <- subset(pData(gs_bcell), Center == center & Sample == 12828 & Replicate == 1)$name
fcs_file <- subset(pData(gs_bcell), Center == center & Sample == 1369 & Replicate == 1)$name
fcs_file <- subset(pData(gs_bcell), Center == center & Sample == 1349 & Replicate == 1)$name

gate <- openCyto:::mindensity(fs_root[[fcs_file]], channel = "FSC-A", gate_range = c(5e4, 1e5), min = 0)


x <- exprs(fs_lymph[[fcs_file]])[, "CD3"]
plot(density(x[x > 0], adjust = 1.5), main = paste("Center:", center))
abline(v = gate@min)

plot(hexbin(exprs(fs_lymph[[fcs_file]])[, c("IgD", "CD27")]))
