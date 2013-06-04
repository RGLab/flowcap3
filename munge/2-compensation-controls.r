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
  path <- file.path(path_Lyoplate, center)

  # The filename of the manually edited Excel file: assumed to be in getwd()
  xlsx <- dir(pattern = center)
  compensation_lyoplate(path = path, xlsx = xlsx, plot = FALSE, pregate = TRUE,
                        min = c(2e4, 5e4), level = 0.99)
})
names(compensation_matrices) <- centers

# If selected, we save the spillover matrices for each center to CSV files.
if (save_CSV) {
  for (i in seq_along(compensation_matrices)) {
    write.csv(compensation_matrices[[i]],
              file = paste0("compensation/", names(compensation_matrices)[i], ".csv"))
  }
}

# Archives the results
save(compensation_matrices, file = "cache/compensation-matrices.RData")

