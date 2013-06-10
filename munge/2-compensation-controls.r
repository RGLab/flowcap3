library(ProjectTemplate)
load.project()

set.seed(42)
path_Lyoplate <- "/loc/no-backup/ramey/Lyoplate"

save_CSV <- FALSE

# There is a bug in the built-in 'list.dirs' function. The argument 'full.names'
# does not work as advertised. After a quick Google search, others recently have
# had similar results.
# To fix this, I manually grab the directory names using 'strsplit'.
centers <- list.dirs(path_Lyoplate, recursive = FALSE, full.names = FALSE)
centers <- sapply(strsplit(centers, split = "/"), tail, n = 1)
centers <- setdiff(centers, "gating-sets")

compensation_matrices <- lapply(centers, function(center) {
  message("Center: ", center)
  center <- "Yale"
  path <- file.path(path_Lyoplate, center)

  # The filename of the manually edited Excel file: assumed to be in getwd()
  xlsx <- dir(pattern = center)
  compensation_lyoplate(path = path, xlsx = xlsx, plot = FALSE, pregate = TRUE, K = 4:5, cells_kept = 0.8,
                        min = c(2e4, 5e4), level = 0.9, plot_markers = FALSE, method = "median", plot_spillover = FALSE)
})
names(compensation_matrices) <- centers

# Manual updating of Yale spillover matrices

# FITC-A
# The compensation control has very few positively expressed cells. Among those
# that are positively expressed, several them have very large SSC-A, and the
# spillover obtained from the FITC-A control does not appear reliable. Instead,
# we use the spillover provided in the FCS keywords.
compensation_matrices$Yale["FITC-A",] <-
  c(0.00604867206091464, 0.0702000061250049, 0, 0.00231576416476302,
  0.0584999991763228, 0.00710027118842876, 1, 0.0050999691448011)

# PerCP-Cy5-5-A
# The compensation control has a small amount of cells expressed at 7.5K units.
# We manually used the small peak at 7.5K to elicit the spillover below.
# Before the manual correction, the spillover matrices greatly differed from
# those in the FCS keywords.
compensation_matrices$Yale["PerCP-Cy5-5-A",] <-
  c(0.00853999938964844, 0.0170800008138021, 0.0506706685384115, 
  0.117954665374756, 0.0141146662394206, 0.0579773328463236, 0.000607999928792318, 1)


# If selected, we save the spillover matrices for each center to CSV files.
if (save_CSV) {
  for (i in seq_along(compensation_matrices)) {
    write.csv(compensation_matrices[[i]],
              file = paste0("compensation/", names(compensation_matrices)[i], ".csv"))
  }
}

# Archives the results
save(compensation_matrices, file = "cache/compensation-matrices.RData")


trace("spillover", signature = "flowSet", tracer = browser)
