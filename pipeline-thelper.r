library(ProjectTemplate)
load.project()

# Cellular subpopulations of interest
subpopulations <- c("lymph", "CD3", "CD4", "CD8",
                    "CD4/CD38+HLADR+", "CD4/CD38-HLADR+", "CD4/CD38+HLADR-", "CD4/CD38-HLADR-",
                    "CD8/CD38+HLADR+", "CD8/CD38-HLADR+", "CD8/CD38+HLADR-", "CD8/CD38-HLADR-",
                    "CD4/CXCR3+CCR6+", "CD4/CXCR3-CCR6+", "CD4/CXCR3+CCR6-", "CD4/CXCR3-CCR6-",
                    "CD8/CXCR3+CCR6+", "CD8/CXCR3-CCR6+", "CD8/CXCR3+CCR6-", "CD8/CXCR3-CCR6-")

panel <- "Thelper"
gs_path <- "/loc/no-backup/ramey/Lyoplate/gating-sets/gs-thelper"

# Loads the archived gatingSet object
gs_thelper <- load_gs(gs_path)

# Creates the gating-template object from a CSV file
gt_csv <- "gt-thelper.csv"
gating_template <- gatingTemplate(gt_csv, panel)

# Applies OpenCyto to GatingSet
set.seed(42)
gating(gating_template, gs_thelper, mc.cores = 10, parallel_type = "multicore") #prior_group = "Center")

# Archives the GatingSet
save_gs(gs_thelper, path = gs_path, overwrite = TRUE)

# Extracts the proportions of the subpopulations of interest and saves them to a CSV file
# For each file, we also add the given pData.
pop_stats <- reshape2::melt(getPopStats(gs_thelper))
colnames(pop_stats) <- c("Population", "Filename", "Proportion")
pop_stats$Population <- as.character(pop_stats$Population)
pop_stats$Filename <- as.character(pop_stats$Filename)

# Splits the strings to extract the populations and their parents
population_split <- strsplit(pop_stats$Population, "/")
parent_pop <- sapply(population_split, function(x) head(tail(x, n = 2), n = 1))
parent_pop <- replace(parent_pop, parent_pop == "root", NA)
parent_pop <- replace(parent_pop, parent_pop == "", "root")
pop_stats$Population <- sapply(population_split, tail, n = 1)
pop_stats$Parent <- parent_pop

which_CD4 <- which(pop_stats$Parent == "CD4")
which_CD8 <- which(pop_stats$Parent == "CD8")
pop_stats <- within(pop_stats,
                    Population[which_CD4] <- paste0("CD4/", Population[which_CD4]))
pop_stats <- within(pop_stats,
                    Population[which_CD8] <- paste0("CD8/", Population[which_CD8]))

# Grabs the counts for each population and its parent
# Then appends the counts to the population statistics
pop_counts <- reshape2::melt(getPopStats(gs_thelper, stat = "count"))
colnames(pop_counts) <- c("Population", "Filename", "Count")
pop_counts$Population <- as.character(pop_counts$Population)
pop_counts$Filename <- as.character(pop_counts$Filename)

# Splits the strings to extract the populations and their parents
population_split <- strsplit(pop_counts$Population, "/")
parent_pop <- sapply(population_split, function(x) head(tail(x, n = 2), n = 1))
parent_pop <- replace(parent_pop, parent_pop == "root", NA)
parent_pop <- replace(parent_pop, parent_pop == "", "root")
pop_counts$Population <- sapply(population_split, tail, n = 1)
pop_counts$Parent <- parent_pop

which_CD4 <- which(pop_counts$Parent == "CD4")
which_CD8 <- which(pop_counts$Parent == "CD8")
pop_counts <- within(pop_counts,
                    Population[which_CD4] <- paste0("CD4/", Population[which_CD4]))
pop_counts <- within(pop_counts,
                    Population[which_CD8] <- paste0("CD8/", Population[which_CD8]))

pop_stats <- merge(pop_stats, pop_counts)

parent_match <- match(pop_stats$Parent, pop_stats$Population)

pop_stats$Count_Parent <- pop_stats$Count[parent_match]

# Subsets and polishes the population statistics
pop_stats <- subset(pop_stats, Population %in% subpopulations)
pop_stats <- merge(pop_stats, pData(gs_thelper), by.x = "Filename", by.y = "name", sort = FALSE)
pop_stats <- subset(pop_stats, select = c(Population, Filename, Center, Sample, Replicate, Proportion, Count, Parent, Count_Parent))

pop_stats <- within(pop_stats, Population[Population == "lymph"] <- "Lymphocytes")

# Replaces CD4 and CD8 with the population names
pop_stats$Population <- sapply(strsplit(pop_stats$Population, "/"), tail, n = 1)
pop_stats <- within(pop_stats, Population[Population == "CD4"] <- "CD4+CD8-")
pop_stats <- within(pop_stats, Population[Population == "CD8"] <- "CD4-CD8+")
pop_stats <- within(pop_stats, Parent[Parent == "CD4"] <- "CD4+CD8-")
pop_stats <- within(pop_stats, Parent[Parent == "CD8"] <- "CD4-CD8+")

write.csv(pop_stats, file = "popstats/popstats-thelper.csv", row.names = FALSE)
