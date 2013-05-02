#' # Lyoplate 3.0 Gates: B-Cell Panel
#'
#' These results are for the first 4 centers that submitted data. We are in the process of incorporating the other centers.

#+ setup, include=FALSE, cache=FALSE, echo=FALSE, warning=FALSE
opts_chunk$set(fig.align = 'default', dev = 'png', message = FALSE, warning = FALSE,
               cache = FALSE, echo = FALSE, fig.path = 'figure/gates-bcell-',
               cache.path = 'cache/gates-bcell-', fig.width = 14, fig.height = 14,
               results = 'hide')

#+ load_data  
setwd("..")
library(ProjectTemplate)
load.project()

panel <- "Bcell"
gs_path <- "/shared/silo_researcher/Gottardo_R/ramey_working/Lyoplate/gs-bcell"

# Loads the archived gatingSet object
library(flowIncubator)
gs_bcell <- load_gs(gs_path)


#' ## Lymph
#+ lymph
plotGate(gs_bcell, 3, lattice = TRUE, xbin = 128,
                  cond = "factor(Center):factor(Sample):factor(Replicate)")

#' ## CD3
#+ cd3
plotGate(gs_bcell, 4, lattice = TRUE, xbin = 128,
                  cond = "factor(Center):factor(Sample):factor(Replicate)")

#' ## CD19, CD20
#+ cd19_cd20
plotGate(gs_bcell, 18:21, lattice = TRUE, xbin = 128,
                  cond = "factor(Center):factor(Sample):factor(Replicate)")

#' ## CD27+IgD+
#+ IgD
plotGate(gs_bcell, 14:17, lattice = TRUE, xbin = 128,
                  cond = "factor(Center):factor(Sample):factor(Replicate)")

#' ## Plasmablast
#+ plasmablast
plotGate(gs_bcell, 11, lattice = TRUE, xbin = 128,
                  cond = "factor(Center):factor(Sample):factor(Replicate)")

#' ## Transitional
#+ transitional
plotGate(gs_bcell, 22, lattice = TRUE, xbin = 128,
                  cond = "factor(Center):factor(Sample):factor(Replicate)")


#+ plot_gate, fig.keep='all',eval=FALSE
plotGate(gs_bcell, lattice = TRUE, xbin = 128,
                  cond = "factor(Center):factor(Sample):factor(Replicate)")
