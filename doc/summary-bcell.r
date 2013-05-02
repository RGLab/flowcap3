#' # Lyoplate 3.0 Summary: B-Cell Panel
#'
#' These results are for the first 4 centers that submitted data. We are in the process of incorporating the other centers.

#+ setup, include=FALSE, cache=FALSE, echo=FALSE, warning=FALSE
opts_chunk$set(fig.align = 'default', dev = 'png', message = FALSE, warning = FALSE,
               cache = TRUE, echo = FALSE, fig.path = 'figure/gates-bcell-',
               cache.path = 'cache/gates-bcell-', fig.width = 14, fig.height = 14,
               results = 'hide')

#+ load_data  
setwd("..")
library(ProjectTemplate)
load.project()

# Loads the archived gatingSet object
library(flowIncubator)
gs_bcell <- load_gs("/shared/silo_researcher/Gottardo_R/ramey_working/Lyoplate/gs-bcell")

#+ popstats_summary
pop_stats <- getPopStats(gs_bcell)
pop_stats <- melt(pop_stats)
colnames(pop_stats) <- c("Marker", "Filename", "Proportion")
pop_stats <- merge(x = pop_stats, y = pData(gs_bcell), by.x = "Filename",
                   by.y = "name")
pop_stats <- subset(pop_stats, Marker != "root")
pop_stats$Marker <- factor(as.character(pop_stats$Marker))

# Summarizes the proportions for each pairing of center and marker.
center_summary <- ddply(pop_stats, .(Center, Marker), summarize,
                        Mean = mean(Proportion), Median = median(Proportion),
                        CV = sd(Proportion) / mean(Proportion))

#' ## Coefficient of Variation
#+ marker_boxplot

# Calculates the CV by marker across all centers and plots them.
CV_by_marker <- ddply(center_summary, .(Marker), summarize, CV = sd(Median) / mean(Median))

p <- ggplot(pop_stats, aes(x = Center, y = Proportion)) + geom_boxplot()
p <- p + facet_wrap( ~ Marker, scales = "free_y") + theme(axis.text.x = element_text(angle = 90))
p


#+ CV_by_center
p <- ggplot(center_summary, aes(x = Center, y = CV)) + geom_bar()
p <- p + facet_wrap( ~ Marker, scales = "free_y")
p <- p + theme(axis.text.x = element_text(angle = 90)) + ylab("Coefficient of Variation")
p


#+ CV_by_marker
p <- ggplot(CV_by_marker, aes(x = Marker, y = CV)) + geom_bar()
p <- p + theme(axis.text.x = element_text(angle = 90)) + ylab("Coefficient of Variation")
p

