library(flowWorkspace)

# Creates directories where plots will be placed
dir.create("gates")
dir.create("gates/bcell")
dir.create("gates/tcell")
dir.create("gates/treg")
dir.create("gates/thelper")
dir.create("gates/DC")

# Constants
plot_width <- plot_height <- 800
the_answer <- 42

# For each panel, we create a subfolder for each pair of markers that are plotted.
# We plot each sample individually.
# The markers are specified in a named list:
# 1. The folder takes the list element's name
# 2. The markers plotted are determined by the list element

# B-cells
plots_path <- "gates/bcell"
gs_path <- "/loc/no-backup/ramey/Lyoplate/gating-sets/gs-bcell"
gs_bcell <- load_gs(gs_path)

markers <- list(Lymphocytes = "lymph", Live = "Live", CD3_CD19 = c("CD3", "CD19"),
                CD20 = "CD20", Plasmablasts = "Plasmablasts", Transitional = "Transitional",
                IgD_CD27 = c("IgD+CD27+", "IgD+CD27-", "IgD-CD27+", "IgD-CD27-"))

for (marker_i in seq_along(markers)) {
  # Creates subfolder for plots
  marker_path <- file.path(plots_path, names(markers)[marker_i])
  dir.create(marker_path)

  for (sample_j in seq_along(gs_bcell)) {
    # The plot has the following format: Center-Sample-Replicate.png
    # Example: Baylor-12828-1.png
    pData_j <- pData(gs_bcell)[sample_j, ]
    plot_filename <- with(pData_j, paste(Center, Sample, Replicate, sep = "-"))
    plot_filename <- paste0(file.path(marker_path, plot_filename), ".png")

    png(file = plot_filename, width = plot_width, height = plot_height)
    plotGate(gs_bcell[[pData_j$name]], markers[[marker_i]], lattice = TRUE, xbin = 128)
    dev.off()
  }
}

rm(gs_bcell)

# T-cells
plots_path <- "gates/tcell"
gs_path <- "/loc/no-backup/ramey/Lyoplate/gating-sets/gs-tcell"
gs_tcell <- load_gs(gs_path)

markers <- list(Lymphocytes = "lymph", Live = "Live", CD3 = "CD3", CD4_CD8 = c("CD4", "CD8"),
                CD4_CD38_HLADR = "CD4/CD38+HLADR+",
                CD8_CD38_HLADR = "CD8/CD38+HLADR+",
                CD4_CCR7_CD45RA = c("CD4/CCR7+CD45RA+", "CD4/CCR7-CD45RA+",
                  "CD4/CCR7+CD45RA-", "CD4/CCR7-CD45RA-"),
                CD8_CCR7_CD45RA = c("CD8/CCR7+CD45RA+", "CD8/CCR7-CD45RA+",
                  "CD8/CCR7+CD45RA-", "CD8/CCR7-CD45RA-"))

for (marker_i in seq_along(markers)) {
  # Creates subfolder for plots
  marker_path <- file.path(plots_path, names(markers)[marker_i])
  dir.create(marker_path)

  for (sample_j in seq_along(gs_tcell)) {
    # The plot has the following format: Center-Sample-Replicate.png
    # Example: Baylor-12828-1.png
    pData_j <- pData(gs_tcell)[sample_j, ]
    plot_filename <- with(pData_j, paste(Center, Sample, Replicate, sep = "-"))
    plot_filename <- paste0(file.path(marker_path, plot_filename), ".png")

    png(file = plot_filename, width = plot_width, height = plot_height)
    plotGate(gs_tcell[[pData_j$name]], markers[[marker_i]], lattice = TRUE, xbin = 128)
    dev.off()
  }
}

rm(gs_tcell)

# T-reg
plots_path <- "gates/treg"
gs_path <- "/loc/no-backup/ramey/Lyoplate/gating-sets/gs-treg"
gs_treg <- load_gs(gs_path)

markers <- list(Lymphocytes = "lymph", Live = "Live", CD3 = "CD3", CD4 = "CD4",
                CD25_CD127 = "CD25+CD127-",
                CCR4_CD45RO = c("CCR4-CD45RO+", "CCR4+CD45RO-", "CCR4-CD45RO-"),
                CCR4_HLADR = c("CCR4+HLADR+", "CCR4-HLADR+", "CCR4+HLADR-", "CCR4-HLADR-"))

for (marker_i in seq_along(markers)) {
  # Creates subfolder for plots
  marker_path <- file.path(plots_path, names(markers)[marker_i])
  dir.create(marker_path)

  for (sample_j in seq_along(gs_treg)) {
    # The plot has the following format: Center-Sample-Replicate.png
    # Example: Baylor-12828-1.png
    pData_j <- pData(gs_treg)[sample_j, ]
    plot_filename <- with(pData_j, paste(Center, Sample, Replicate, sep = "-"))
    plot_filename <- paste0(file.path(marker_path, plot_filename), ".png")

    png(file = plot_filename, width = plot_width, height = plot_height)
    plotGate(gs_treg[[pData_j$name]], markers[[marker_i]], lattice = TRUE, xbin = 128)
    dev.off()
  }
}

rm(gs_treg)

# T-helper
plots_path <- "gates/thelper"
gs_path <- "/loc/no-backup/ramey/Lyoplate/gating-sets/gs-thelper"
gs_thelper <- load_gs(gs_path)

markers <- list(Lymphocytes = "lymph", Live = "Live", CD3 = "CD3", CD4_CD8 = c("CD4", "CD8"),
                CD4_CD38_HLADR = c("CD4/CD38+HLADR+", "CD4/CD38-HLADR+",
                  "CD4/CD38+HLADR-", "CD4/CD38-HLADR-"),
                CD8_CD38_HLADR = c("CD8/CD38+HLADR+", "CD8/CD38-HLADR+",
                  "CD8/CD38+HLADR-", "CD8/CD38-HLADR-"),
                CD4_CXCR3_CCR6 = c("CD4/CXCR3+CCR6+", "CD4/CXCR3-CCR6+",
                  "CD4/CXCR3+CCR6-", "CD4/CXCR3-CCR6-"),
                CD8_CXCR3_CCR6 = c("CD8/CXCR3+CCR6+", "CD8/CXCR3-CCR6+",
                  "CD8/CXCR3+CCR6-", "CD8/CXCR3-CCR6-"))

for (marker_i in seq_along(markers)) {
  # Creates subfolder for plots
  marker_path <- file.path(plots_path, names(markers)[marker_i])
  dir.create(marker_path)

  for (sample_j in seq_along(gs_thelper)) {
    # The plot has the following format: Center-Sample-Replicate.png
    # Example: Baylor-12828-1.png
    pData_j <- pData(gs_thelper)[sample_j, ]
    plot_filename <- with(pData_j, paste(Center, Sample, Replicate, sep = "-"))
    plot_filename <- paste0(file.path(marker_path, plot_filename), ".png")

    png(file = plot_filename, width = plot_width, height = plot_height)
    plotGate(gs_thelper[[pData_j$name]], markers[[marker_i]], lattice = TRUE, xbin = 128)
    dev.off()
  }
}

rm(gs_thelper)


# DC/Mono/NK
plots_path <- "gates/DC"
gs_path <- "/loc/no-backup/ramey/Lyoplate/gating-sets/gs-DC"
gs_DC <- load_gs(gs_path)

markers <- list(Monocytes = "Monocytes", Live = "Live", HLADR = "HLADR+",
                CD14_Lineage = c("CD14-Lineage-", "CD14+Lineage-"), CD14_CD16 = "CD14+CD16+",
                CD16_CD56 = c("CD16+CD56+", "CD16+CD56-", "CD16-CD56+", "CD16-CD56-"),
                CD11c_CD123 = c("CD11c+CD123+", "CD11c+CD123-", "CD11c-CD123+", "CD11c-CD123-"))
               
for (marker_i in seq_along(markers)) {
  # Creates subfolder for plots
  marker_path <- file.path(plots_path, names(markers)[marker_i])
  dir.create(marker_path)

  for (sample_j in seq_along(gs_DC)) {
    # The plot has the following format: Center-Sample-Replicate.png
    # Example: Baylor-12828-1.png
    pData_j <- pData(gs_DC)[sample_j, ]
    plot_filename <- with(pData_j, paste(Center, Sample, Replicate, sep = "-"))
    plot_filename <- paste0(file.path(marker_path, plot_filename), ".png")

    png(file = plot_filename, width = plot_width, height = plot_height)
    plotGate(gs_DC[[pData_j$name]], markers[[marker_i]], lattice = TRUE, xbin = 128)
    dev.off()
  }
}

rm(gs_DC)
