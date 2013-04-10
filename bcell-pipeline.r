library(ProjectTemplate)
load.project()

path_Lyoplate <- "/loc/no-backup/ramey/Lyoplate/"
panel <- "Bcell"

gs_file <- file.path(path_Lyoplate, "gs-Bcell.tar")

# Loads the archived gatingSet object
gs_bcell <- unarchive(file = gs_file, path = path_Lyoplate)

# Creates the gating-template object from a CSV file
# Applies OpenCyto to GatingSet
gt_csv <- "gt-bcell.csv"
gating_template <- gatingTemplate(gt_csv, panel)

# TODO: Ensure that prior arguments are passed to prior_flowClust
# TODO: Relax the requirement that K must be specified



fs <- getData(gs_bcell)




trace("gating", signature = c("gtMethod", "GatingSet"), tracer = browser)
untrace("gating", signature = c("gtMethod", "GatingSet"))


set.seed(42)
prior <- prior_flowClust1d(fs, channel = "FSC-A", min = 0)


gates <- fsApply(fs, flowClust.1d, params = "FSC-A", filterId = "nonDebris",
                                  prior = prior, K = length(prior$Mu0), min = 0)
sapply(gates, function(x) x@min)

densityplot(~ `FSC-A`, fs, filter = gates)


plotGate(gs_bcell, 2, lattice = TRUE, xbin = 128,
                  cond = "factor(Center):factor(Sample):factor(Replicate)")














# Manually applies more accurate IgD gate that OpenCyto does not currently
# support.
Rm("IgD-cd27+", gs_bcell)
Rm("IgD+cd27+", gs_bcell)
Rm("IgD+cd27-", gs_bcell)
Rm("IgD-cd27-", gs_bcell)
gating.bsub(gs_bcell, parent = "cd19 & cd20")


gating(gating_template, gs_bcell, num_nodes = 12, parallel_type = "multicore")
archive(gs_bcell, file = gs_file)
