#' # Pairwise Plots of Markers of Interest after Applying FCStrans to Lyoplate 3.0
#' ## B-Cell Panel
#'
#' The data have been preprocessed with debris and lymphocyte gates.

#+ setup, include=FALSE, cache=FALSE, echo=FALSE, warning=FALSE
opts_chunk$set(fig.align = 'default', dev = 'png', message = FALSE, warning = FALSE,
               cache = TRUE, echo = FALSE, fig.path = 'figure/FCStrans-',
               cache.path = 'cache/FCStrans-', fig.width = 12, fig.height = 12,
               results = 'hide')

#+ load_data  
setwd("..")
library(ProjectTemplate)
load.project()

path_Lyoplate <- "/loc/no-backup/ramey/Lyoplate/"
panel <- "Bcell"

gs_file <- file.path(path_Lyoplate, "gs-Bcell.tar")

# Loads the archived gatingSet object
gs_bcell <- unarchive(file = gs_file, path = path_Lyoplate)

#+ setup_data
fs <- getData(gs_bcell, 3)

fs_miami <- fs[1:9]
fs_nhlbi <- fs[10:18]
fs_stanford <- fs[19:27]
fs_ucla <- fs[28:36]

# QUESTION: Did the 4 centers record FSC-H, so that we can apply a singlet gate?

#' ## Miami

#+ miami_scatter
par(mfrow = c(3, 3))
lapply(seq_along(fs_miami), function(i) {
  plot(exprs(fs_miami[[i]][, c("FSC-A", "SSC-A")]))
})

#+ miami_cd3
par(mfrow = c(3, 3))
lapply(seq_along(fs_miami), function(i) {
  plot(exprs(fs_miami[[i]][, c("CD3", "SSC-A")]))
})

#+ miami_cd19_cd20
par(mfrow = c(3, 3))
lapply(seq_along(fs_miami), function(i) {
  plot(exprs(fs_miami[[i]][, c("CD19", "CD20")]))
})

#+ miami_cd38_cd27
par(mfrow = c(3, 3))
lapply(seq_along(fs_miami), function(i) {
  plot(exprs(fs_miami[[i]][, c("CD38", "CD27")]))
})

#+ miami_IgD_cd27
par(mfrow = c(3, 3))
lapply(seq_along(fs_miami), function(i) {
  plot(exprs(fs_miami[[i]][, c("IgD", "CD27")]))
})

#+ miami_cd38_cd24
par(mfrow = c(3, 3))
lapply(seq_along(fs_miami), function(i) {
  plot(exprs(fs_miami[[i]][, c("CD38", "CD24")]))
})

#' ## NHLBI

#+ nhlbi_scatter
par(mfrow = c(3, 3))
lapply(seq_along(fs_nhlbi), function(i) {
  plot(exprs(fs_nhlbi[[i]][, c("FSC-A", "SSC-A")]))
})

#+ nhlbi_cd3
par(mfrow = c(3, 3))
lapply(seq_along(fs_nhlbi), function(i) {
  plot(exprs(fs_nhlbi[[i]][, c("CD3", "SSC-A")]))
})

#+ nhlbi_cd19_cd20
par(mfrow = c(3, 3))
lapply(seq_along(fs_nhlbi), function(i) {
  plot(exprs(fs_nhlbi[[i]][, c("CD19", "CD20")]))
})

#+ nhlbi_cd38_cd27
par(mfrow = c(3, 3))
lapply(seq_along(fs_nhlbi), function(i) {
  plot(exprs(fs_nhlbi[[i]][, c("CD38", "CD27")]))
})

#+ nhlbi_IgD_cd27
par(mfrow = c(3, 3))
lapply(seq_along(fs_nhlbi), function(i) {
  plot(exprs(fs_nhlbi[[i]][, c("IgD", "CD27")]))
})

#+ nhlbi_cd38_cd24
par(mfrow = c(3, 3))
lapply(seq_along(fs_nhlbi), function(i) {
  plot(exprs(fs_nhlbi[[i]][, c("CD38", "CD24")]))
})

#' ## Stanford

#+ stanford_scatter
par(mfrow = c(3, 3))
lapply(seq_along(fs_stanford), function(i) {
  plot(exprs(fs_stanford[[i]][, c("FSC-A", "SSC-A")]))
})

#+ stanford_cd3
par(mfrow = c(3, 3))
lapply(seq_along(fs_stanford), function(i) {
  plot(exprs(fs_stanford[[i]][, c("CD3", "SSC-A")]))
})

#+ stanford_cd19_cd20
par(mfrow = c(3, 3))
lapply(seq_along(fs_stanford), function(i) {
  plot(exprs(fs_stanford[[i]][, c("CD19", "CD20")]))
})

#+ stanford_cd38_cd27
par(mfrow = c(3, 3))
lapply(seq_along(fs_stanford), function(i) {
  plot(exprs(fs_stanford[[i]][, c("CD38", "CD27")]))
})

#+ stanford_IgD_cd27
par(mfrow = c(3, 3))
lapply(seq_along(fs_stanford), function(i) {
  plot(exprs(fs_stanford[[i]][, c("IgD", "CD27")]))
})

#+ stanford_cd38_cd24
par(mfrow = c(3, 3))
lapply(seq_along(fs_stanford), function(i) {
  plot(exprs(fs_stanford[[i]][, c("CD38", "CD24")]))
})

#' ## UCLA

#+ ucla_scatter
par(mfrow = c(3, 3))
lapply(seq_along(fs_ucla), function(i) {
  plot(exprs(fs_ucla[[i]][, c("FSC-A", "SSC-A")]))
})

#+ ucla_cd3
par(mfrow = c(3, 3))
lapply(seq_along(fs_ucla), function(i) {
  plot(exprs(fs_ucla[[i]][, c("CD3", "SSC-A")]))
})

#+ ucla_cd19_cd20
par(mfrow = c(3, 3))
lapply(seq_along(fs_ucla), function(i) {
  plot(exprs(fs_ucla[[i]][, c("CD19", "CD20")]))
})

#+ ucla_cd38_cd27
par(mfrow = c(3, 3))
lapply(seq_along(fs_ucla), function(i) {
  plot(exprs(fs_ucla[[i]][, c("CD38", "CD27")]))
})

#+ ucla_IgD_cd27
par(mfrow = c(3, 3))
lapply(seq_along(fs_ucla), function(i) {
  plot(exprs(fs_ucla[[i]][, c("IgD", "CD27")]))
})

#+ ucla_cd38_cd24
par(mfrow = c(3, 3))
lapply(seq_along(fs_ucla), function(i) {
  plot(exprs(fs_ucla[[i]][, c("CD38", "CD24")]))
})


