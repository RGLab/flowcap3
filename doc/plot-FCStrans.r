#' # Pairwise Plots of Markers of Interest after Applying FCStrans to Lyoplate 3.0
#' ## B-Cell Panel
#'
#' The data have been preprocessed with debris and lymphocyte gates.

#+ setup, include=FALSE, cache=FALSE, echo=FALSE, warning=FALSE
opts_chunk$set(fig.align = 'default', dev = 'png', message = FALSE, warning = FALSE,
               cache = FALSE, echo = FALSE, fig.path = 'figure/FCStrans-',
               cache.path = 'cache/FCStrans-', fig.width = 12, fig.height = 12,
               results = 'hide')

#+ load_data  
setwd("..")
library(ProjectTemplate)
load.project()

# Loads the archived gatingSet object
library(flowIncubator)
gs_bcell <- load_gs("/shared/silo_researcher/Gottardo_R/ramey_working/Lyoplate/gs-bcell")

#+ setup_data
fs <- getData(gs_bcell, 3)

fs_baylor <- fs[1:9]
fs_miami <- fs[10:18]
fs_nhlbi <- fs[19:27]
fs_stanford <- fs[28:36]
fs_ucla <- fs[37:45]


#' ## Baylor

#' ### FSC-A vs SSC-A
#+ baylor_scatter
lapply(seq_along(fs_baylor), function(i) {
  marginal_gating_plot(data = exprs(fs_baylor[[i]]), feature_pairs = c("FSC-A", "SSC-A"))
})

#' ### CD3
#+ baylor_cd3
lapply(seq_along(fs_baylor), function(i) {
  marginal_gating_plot(data = exprs(fs_baylor[[i]]), feature_pairs = c("CD3", "SSC-A"))
})

#' ### CD19 vs CD20
#+ baylor_cd19_cd20
lapply(seq_along(fs_baylor), function(i) {
  marginal_gating_plot(data = exprs(fs_baylor[[i]]), feature_pairs = c("CD19", "CD20"))
})

#' ### CD38 vs CD27
#+ baylor_cd38_cd27
lapply(seq_along(fs_baylor), function(i) {
  marginal_gating_plot(data = exprs(fs_baylor[[i]]), feature_pairs = c("CD38", "CD27"))
})

#' ### IgD vs CD27
#+ baylor_IgD_cd27
lapply(seq_along(fs_baylor), function(i) {
  marginal_gating_plot(data = exprs(fs_baylor[[i]]), feature_pairs = c("IgD", "CD27"))
})

#' ### CD38 vs CD24
#+ baylor_cd38_cd24
lapply(seq_along(fs_baylor), function(i) {
  marginal_gating_plot(data = exprs(fs_baylor[[i]]), feature_pairs = c("CD38", "CD24"))
})


#' ## Miami

#' ### FSC-A vs SSC-A
#+ miami_scatter
lapply(seq_along(fs_miami), function(i) {
  marginal_gating_plot(data = exprs(fs_miami[[i]]), feature_pairs = c("FSC-A", "SSC-A"))
})

#' ### CD3
#+ miami_cd3
lapply(seq_along(fs_miami), function(i) {
  marginal_gating_plot(data = exprs(fs_miami[[i]]), feature_pairs = c("CD3", "SSC-A"))
})

#' ### CD19 vs CD20
#+ miami_cd19_cd20
lapply(seq_along(fs_miami), function(i) {
  marginal_gating_plot(data = exprs(fs_miami[[i]]), feature_pairs = c("CD19", "CD20"))
})

#' ### CD38 vs CD27
#+ miami_cd38_cd27
lapply(seq_along(fs_miami), function(i) {
  marginal_gating_plot(data = exprs(fs_miami[[i]]), feature_pairs = c("CD38", "CD27"))
})

#' ### IgD vs CD27
#+ miami_IgD_cd27
lapply(seq_along(fs_miami), function(i) {
  marginal_gating_plot(data = exprs(fs_miami[[i]]), feature_pairs = c("IgD", "CD27"))
})

#' ### CD38 vs CD24
#+ miami_cd38_cd24
lapply(seq_along(fs_miami), function(i) {
  marginal_gating_plot(data = exprs(fs_miami[[i]]), feature_pairs = c("CD38", "CD24"))
})


#' ## NHLBI

#' ### FSC-A vs SSC-A
#+ nhlbi_scatter
lapply(seq_along(fs_nhlbi), function(i) {
  marginal_gating_plot(data = exprs(fs_nhlbi[[i]]), feature_pairs = c("FSC-A", "SSC-A"))
})

#' ### CD3
#+ nhlbi_cd3
lapply(seq_along(fs_nhlbi), function(i) {
  marginal_gating_plot(data = exprs(fs_nhlbi[[i]]), feature_pairs = c("CD3", "SSC-A"))
})

#' ### CD19 vs CD20
#+ nhlbi_cd19_cd20
lapply(seq_along(fs_nhlbi), function(i) {
  marginal_gating_plot(data = exprs(fs_nhlbi[[i]]), feature_pairs = c("CD19", "CD20"))
})

#' ### CD38 vs CD27
#+ nhlbi_cd38_cd27
lapply(seq_along(fs_nhlbi), function(i) {
  marginal_gating_plot(data = exprs(fs_nhlbi[[i]]), feature_pairs = c("CD38", "CD27"))
})

#' ### IgD vs CD27
#+ nhlbi_IgD_cd27
lapply(seq_along(fs_nhlbi), function(i) {
  marginal_gating_plot(data = exprs(fs_nhlbi[[i]]), feature_pairs = c("IgD", "CD27"))
})

#' ### CD38 vs CD24
#+ nhlbi_cd38_cd24
lapply(seq_along(fs_nhlbi), function(i) {
  marginal_gating_plot(data = exprs(fs_nhlbi[[i]]), feature_pairs = c("CD38", "CD24"))
})


#' ## Stanford

#' ### FSC-A vs SSC-A
#+ stanford_scatter
lapply(seq_along(fs_stanford), function(i) {
  marginal_gating_plot(data = exprs(fs_stanford[[i]]), feature_pairs = c("FSC-A", "SSC-A"))
})

#' ### CD3
#+ stanford_cd3
lapply(seq_along(fs_stanford), function(i) {
  marginal_gating_plot(data = exprs(fs_stanford[[i]]), feature_pairs = c("CD3", "SSC-A"))
})

#' ### CD19 vs CD20
#+ stanford_cd19_cd20
lapply(seq_along(fs_stanford), function(i) {
  marginal_gating_plot(data = exprs(fs_stanford[[i]]), feature_pairs = c("CD19", "CD20"))
})

#' ### CD38 vs CD27
#+ stanford_cd38_cd27
lapply(seq_along(fs_stanford), function(i) {
  marginal_gating_plot(data = exprs(fs_stanford[[i]]), feature_pairs = c("CD38", "CD27"))
})

#' ### IgD vs CD27
#+ stanford_IgD_cd27
lapply(seq_along(fs_stanford), function(i) {
  marginal_gating_plot(data = exprs(fs_stanford[[i]]), feature_pairs = c("IgD", "CD27"))
})

#' ### CD38 vs CD24
#+ stanford_cd38_cd24
lapply(seq_along(fs_stanford), function(i) {
  marginal_gating_plot(data = exprs(fs_stanford[[i]]), feature_pairs = c("CD38", "CD24"))
})


#' ## UCLA

#' ### FSC-A vs SSC-A
#+ ucla_scatter
lapply(seq_along(fs_ucla), function(i) {
  marginal_gating_plot(data = exprs(fs_ucla[[i]]), feature_pairs = c("FSC-A", "SSC-A"))
})

#' ### CD3
#+ ucla_cd3
lapply(seq_along(fs_ucla), function(i) {
  marginal_gating_plot(data = exprs(fs_ucla[[i]]), feature_pairs = c("CD3", "SSC-A"))
})

#' ### CD19 vs CD20
#+ ucla_cd19_cd20
lapply(seq_along(fs_ucla), function(i) {
  marginal_gating_plot(data = exprs(fs_ucla[[i]]), feature_pairs = c("CD19", "CD20"))
})

#' ### CD38 vs CD27
#+ ucla_cd38_cd27
lapply(seq_along(fs_ucla), function(i) {
  marginal_gating_plot(data = exprs(fs_ucla[[i]]), feature_pairs = c("CD38", "CD27"))
})

#' ### IgD vs CD27
#+ ucla_IgD_cd27
lapply(seq_along(fs_ucla), function(i) {
  marginal_gating_plot(data = exprs(fs_ucla[[i]]), feature_pairs = c("IgD", "CD27"))
})

#' ### CD38 vs CD24
#+ ucla_cd38_cd24
lapply(seq_along(fs_ucla), function(i) {
  marginal_gating_plot(data = exprs(fs_ucla[[i]]), feature_pairs = c("CD38", "CD24"))
})

