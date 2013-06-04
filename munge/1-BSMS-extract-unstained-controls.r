library(ProjectTemplate)
load.project()

# Path of the FCS files for the center selected
path_Lyoplate <- "/loc/no-backup/ramey/Lyoplate"
center <- "BSMS"
path_center <- file.path(path_Lyoplate, center)

# The filenames of the compensation control FCS files
compensation_files <- dir(path_center, "Compensation.*fcs", full.names = TRUE)
comp_flowset <- read.flowSet(compensation_files)

# Applies flowClust to SSC-A and FSC-A with K = 3 to remove outliers and debris
# Keeps the densest cluster and gates out the rest of the cells.
set.seed(42)
flow_frame <- comp_flowset[[1]]

tmix_filter <- tmixFilter(filterId = "Lymph", parameters = c("FSC-A", "SSC-A"), K = 3)
tmix_results <- filter(flow_frame, tmix_filter)

flow_frame <- split(x = flow_frame, f = tmix_results, population = which.max(tmix_results@w))[[1]]

# For these compensation controls, the single-stained controls can double as
# unstained controls because the positive (stained) and negative peaks are well
# defined. We extract the negative peaks from one of the controls and use it as
# an unstained control.
channel <- "Blue 5-A"
gate_mindensity <- openCyto:::mindensity(flow_frame, channel, positive = FALSE)
flow_frame <- Subset(flow_frame, gate_mindensity)

# Writes the extracted unstained compensation control flowFrame to a FCS file
FCS_filename <- file.path(path_center, "Compensation Controls_Unstained Control.fcs")
write.FCS(x = flow_frame, filename = FCS_filename)

# Updates file permissions and group
Sys.chmod(FCS_filename, mode = "0775")
system(paste0("chgrp hyrax_gottardo '", FCS_filename, "'"))
