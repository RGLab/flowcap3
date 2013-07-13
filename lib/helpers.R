#' Updates the pData for the Cytotrol data.
#'
#' We expect the \code{name} list item in the \code{pData} information in the
#' given \code{flow_set} to have the following format:
#' Center-Bcell-#.fcs
#' For example, the second sample for NHLBI is named NHLBI-Bcell-2.fcs.
#' @param flow_set a \code{flowSet} object
#' @return a \code{flowSet} object with updated \code{pData} information.
update_phenodata <- function(flow_set) {
  pheno_data <- pData(flow_set)

  split_name <- strsplit(pheno_data$name, "-")
  parsed_data <- do.call(rbind.data.frame, split_name)
  colnames(parsed_data) <- c("Center", "Panel", "Replicate")
  pheno_data <- cbind(pheno_data, parsed_data)

  # The Center and Panel information are stored as factors. We update them to
  # character vectors.
  pheno_data$Center <- as.character(pheno_data$Center)
  pheno_data$Panel <- as.character(pheno_data$Panel)

  # The Replicate information are currently stored as the replicate number with
  # ".fcs" appended. For instance, the second replicate is "2.fcs". We remove the
  # ".fcs" suffix and convert the resulting Replicates to integers.
  replicate <- strsplit(as.character(pheno_data$Replicate), ".", fixed = TRUE)
  pheno_data$Replicate <- as.integer(do.call(rbind, replicate)[, 1])

  pData(flow_set) <- pheno_data
  
  var_meta <- varMetadata(phenoData(flow_set))
  var_meta[-1, ] <- rownames(var_meta)[-1]
  varMetadata(phenoData(flow_set)) <- var_meta

  flow_set
}

#' Preprocesses a Cytotrol flowFrame object
#'
#' Our goal here is to use swap the marker names and the channel names within a
#' \code{flowFrame} object to ensure that the \code{flowFrame} objects across
#' centers can be merged into a single \code{flowSet}.
#'
#' We also preprocess the marker names to strip out any additional information
#' added to the marker name. For instance, NHLBI uses "IgD V500", which we reduce
#' to "IgD".
#'
#' @param flow_frame the \code{flowFrame} object to preprocess
#' @param markers_keep a character vector containing the markers to keep
#' @return the updated \code{flowFrame} object containing only the markers of
#' interest
preprocess_flowframe <- function(flow_frame, markers_keep) {
  if (missing(markers_keep)) {
    stop("The marker to keep must be specified.")
  }

  fr_rownames <- rownames(parameters(flow_frame)@data)
  
  # Preprocesses each of the columns in the flow_frame
  for (j in seq_len(ncol(flow_frame))) {
    marker_idx <- paste0(fr_rownames[j], "S")
    channel_idx <- paste0(fr_rownames[j], "N")

    marker <- flow_frame@description[[marker_idx]]
    channel <- flow_frame@description[[channel_idx]]

    # In the case the marker name is given, we swap the marker and channel
    # names.
    if (!is.null(marker)) {
      # Converts the marker names to a common name
      marker <- marker_conversion(marker)

      # Updates the channel information in the flow_frame with the marker
      flow_frame@description[[channel_idx]] <- marker
      flow_frame@parameters@data$name[j] <- marker

      # Updates the marker information in the flow_frame with the channel
      flow_frame@description[[marker_idx]] <- channel
      flow_frame@parameters@data$desc[j] <- channel
    }
  }
  colnames(exprs(flow_frame)) <- colnames(flow_frame)
  
  # Subset to markers of interest
  flow_frame[, markers_keep]
}

#' Converts the Lyoplate marker names to a common name
#'
#' For the following list of marker names, we manually update the names so
#' that they are standard across centers.
marker_conversion <- Vectorize(function(marker) {
  # If marker name contains additional info, remove everything after the
  # space. (e.g., "IgD V500" to "IgD")
  marker <- strsplit(marker, " ")[[1]][1]

  if (marker == "19") {
    marker <- "CD19"
  } else if (marker %in% c("LIVE", "LIVE_GREEN", "Live/Dead")) {
    marker <- "Live"
  } else if (marker == "IGD") {
    marker <- "IgD"
  } else if (marker %in% c("HLA", "HLADR", "HLA-DR")) {
    marker <- "HLADR"
  } else if (marker == "CD197") {
    marker <- "CCR7"
  } else if (marker == "CD194") {
    marker <- "CCR4"
  } else if (marker == "CD11C") {
    marker <- "CD11c"
  } else if (marker %in% c("CD3CD19CD20", "CD3+19+20", "CD3_CD19_CD20",
                           "CD3+CD19+CD20+", "CD3+CD19+CD20", "CD3+19+20")) {
    marker <- "Lineage"
  } else if (marker == "CD196") {
    marker <- "CCR6"
  } else if (marker == "CD183") {
    marker <- "CXCR3"
  }

  marker
})


#' Prettifies the table returned by population stats for the pipeline
#'
#' This function updates the rownames to retain only the number of nodes
#' specified. For example, suppose the full path to the node is
#' /singlet/viable/Lymph/CD3/CD8/IFNg. By default, this is updated to IFNg.
#' If instead we specify \code{nodes = 2}, then the corresponding row is
#' CD8/IFNg.
#' 
#' @param population_stats a data frame containing the population summaries,
#' which is returned from the \code{getPopStats} function in the
#' \code{flowWorkspace} package.
#' @param nodes numeric. How many nodes should be preserved in the prettified
#' population statistics data frame returned? Default: 1
#' @return a data frame containing the population statistics but with prettified
#' row names
pretty_popstats <- function(population_stats, nodes = 1) {
  node_names <- sapply(strsplit(rownames(population_stats), split = "/"), tail, n = nodes)
  rownames(population_stats) <- sapply(node_names, paste, collapse = "/")

  population_stats
}

#' Constructs compensation matrix from controls for Lyoplate 3.0
#'
#' From the compensation control files listed in the given Excel file, we
#' construct a spillover matrix.
#'
#' Given that the compensation controls likely have debris, we gate these out
#' using \code{\link{flowClust}}. 
#'
#' @param path the path that contains the compensation control FCS files
#' @param xlsx the name of the XLSX (Excel) file (full path)
#' @param pregate Should the compensation controls be pregating using
#' \code{flowClust} before constructing the spillover matrices? (Default: Yes)
#' @param plot Should the \code{flowClust} summary figures be plotted?
#' @param method the summary statistic to use in constructing the spillover
#' matrix. This value is passed along to \code{\link[flowCore]{spillover}}.
#' @param additional arguments passed to \code{pregate_flowset}
#' @return the compensation matrix
compensation_lyoplate <- function(path, xlsx, pregate = TRUE, plot = FALSE,
                                  method = c("mode", "median", "mean"),
                                  plot_debris = FALSE, plot_markers = FALSE,
                                  plot_spillover = FALSE, ...) {
                                  
  method <- match.arg(method)
  
  comp_controls <- read.xlsx2(file = xlsx, sheetName = "Comp controls",
                              startRow = 3, stringsAsFactors = FALSE)

  # In the manually edited Excel files, we use the format "Marker:Channel". For
  # instance, we might have CD8:APC-H7. We extract the entire string after the
  # colon.
  comp_controls$Marker <- sapply(strsplit(comp_controls$Marker, split = ":"),
                                 tail, n = 1)

  # For some centers, such as Stanford, generic color controls were used for
  # compensation. In some of these cases, a single FCS file was used for all
  # markers sharing the same channel. We handle this case by replacing the empty
  # FCS filenames with a filename that corresponds to the same channel.
  which_empty <- which(comp_controls$FCS.file.name == "")
  generic_FCS_filenames <- sapply(which_empty, function(i) {
    which_generic <- which(with(comp_controls,
                           Marker[i] == Marker & FCS.file.name != ""))

    # In the case that there are multiple matches, we take the first.
    comp_controls$FCS.file.name[which_generic[1]]
  })
  if (length(generic_FCS_filenames) != 0) {
    comp_controls$FCS.file.name <- replace(comp_controls$FCS.file.name,
                                           which_empty, generic_FCS_filenames)
  }

  # For some markers (e.g., APC-H7 from Stanford), there are no FCS files for a
  # generic marker. In this case, the above generic names result in 'NA', in
  # which case, we omit them from the compensation controls data.frame.
  comp_controls <- subset(comp_controls, !is.na(FCS.file.name))

  # Some centers use one compensation control for several markers that share the
  # same channel (e.g., NHLBI), while other centers use one compensation control
  # per marker (e.g., CIMR). To handle this, we use only the first FCS file if
  # the channel names are non-unique.
  comp_controls <- subset(comp_controls, !duplicated(comp_controls$Marker))

  # Constructs a flowSet of the compensation control FCS files along with the
  # unstained control FCS file
  FCS_files <- comp_controls$FCS.file.name
  unstained_control <- dir(path, pattern = "Compensation.*Unstained")

  if (length(unstained_control) == 0) {
    stop("No unstained compensation control was provided. Not supported by flowCore")
  }

  FCS_files <- file.path(path, c(FCS_files, unstained_control))

  comp_flowset <- read.flowSet(FCS_files)

  oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))

  if (plot_debris) {
    marginal_plots <- lapply(seq_along(comp_flowset), function(i) {
      message("Compensation Control: ", i)
      marginal_gating_plot(data = exprs(comp_flowset[[i]]), feature_pairs = c("FSC-A", "SSC-A"))
    })
  }

  # Applies flowClust to SSC-A and FSC-A to remove outliers and debris before
  # constructing a compensation matrix. Keeps the densest clusters and gates out
  # the rest of the cells.
  if (pregate) {
    comp_flowset <- pregate_flowset(comp_flowset, plot = plot, ...)
  }

  if (plot_markers) {
    for (i in seq_along(comp_flowset)) {
      filename <- keyword(comp_flowset[[i]])$GUID
      message("Compensation Control Filename: ", filename)
      x <- exprs(comp_flowset[[i]])
      m_x <- melt(x)[, -1]
      colnames(m_x) <- c("Marker", "value")
      p <- ggplot(m_x, aes(x = value)) + geom_density() + facet_wrap(~ Marker, scales = "free")
      plot(p + ggtitle(filename))
    }
  }
  
  # Constructs the compensation (spillover) matrix
  # The index of 'unstained' is the last of the FCS files read in.
  # A regular expression in 'patt' indicates the markers that we want to compensate.
  #    These are given in the XLSX file.

  # Use all columns except for FSCA, SSC, and Time for compensation
  comp_colnames <- colnames(comp_flowset)
  comp_colnames <- comp_colnames[!grepl("FSC|SSC|Time", comp_colnames)]

  spillover(comp_flowset, unstained = length(comp_flowset),
            patt = paste(comp_controls$Marker, collapse = "|"),
            useNormFilt = FALSE, method = method, stain_match = "regexpr",
            pregate = TRUE, plot = plot_spillover)
}

#' Constructs a flowSet for the specified center for Lyoplate 3.0
#'
#' Creates a \code{flowSet} by reading the FCS files for the current center given
#' in the Excel file. Then, compensates and transforms them. Finally, constructs
#' and returns the resulting \code{flowSet} object.
#'
#' @param path the path that contains the sample FCS files
#' @param xlsx the name of the XLSX (Excel) file (full path)
#' @param comp_matrix a compensation matrix
#' @param center a character string denoting the current center
#' @param panel the panel type
#' @return a \code{flowSet} object
flowset_lyoplate <- function(path, xlsx, comp_matrix,
                             transform = c("FCSTrans", "estimateLogicle"), center,
                             panel = c("B cell", "T cell", "Treg", "DC/mono/NK", "Th1/2/17"),
                             channels_remove = NULL) {
  panel <- match.arg(panel)
  transform <- match.arg(transform)

  exp_samples <- read.xlsx2(file = xlsx, sheetName = "Exp samples",
                            startRow = 3, stringsAsFactors = FALSE)
  colnames(exp_samples) <- c("name", "Institution", "Panel", "Replicate",
                             "Sample")
  exp_samples <- subset(exp_samples, select = -Institution)
  exp_samples$Replicate <- as.character(as.numeric(exp_samples$Replicate))
  exp_samples$Sample <- as.character(as.numeric(exp_samples$Sample))
  exp_samples$Center <- center

  # Subsets the data for the specified panel
  exp_samples <- exp_samples[exp_samples$Panel == panel, ]

  # Reads in the FCS files for the specified panel
  fcs_files <- file.path(path, exp_samples$name)

  if (!is.null(channels_remove)) {
    channels_regexp <- paste0(channels_remove, collapse = "|")
  } else {
    channels_regexp <- NULL
  }
  exp_flowset <- read.flowSet(fcs_files, column.pattern = channels_regexp, invert = TRUE)

  # Updates the flowSet's pData
  exp_pdata <- pData(exp_flowset)
  exp_pdata <- merge(exp_pdata, exp_samples)
  rownames(exp_pdata) <- rownames(pData(exp_flowset))
  pData(exp_flowset) <- exp_pdata 
  varM <- varMetadata(phenoData(exp_flowset))
  varM[-1,] <- rownames(varM)[-1]
  varMetadata(phenoData(exp_flowset)) <- varM

  # Compensates and transforms flowSet for current center
  if (is.null(comp_matrix)) {
    comp_matrix <- keyword(exp_flowset[[1]])$SPILL
  }
  exp_flowset <- compensate(exp_flowset, comp_matrix)

  # Saves list of estimated FCSTrans widths
  width_list <- list()

  for (channel in colnames(comp_matrix)) {
    # Based on Parks et al. (2006) in Cytometry A, we select "the fifth
    # percentile of all events that are below zero as this reference value."
    quantiles_channel <- fsApply(exp_flowset, function(flow_frame) {
      quantile_negatives(exprs(flow_frame)[, channel], probs = 0.05, names = FALSE)
    })
    quantile_median <- median(quantiles_channel, na.rm = TRUE)
    w_channel <- abs(0.5 * (4.5 * log(10) - log(2^18 / abs(quantile_median))))

    width_list[[channel]] <- w_channel
  }

  widths <- melt(width_list)
  colnames(widths) <- c("width", "Channel")
  widths$Center <- center

  list(flow_set = exp_flowset, widths = widths)
}

# Code modified from:
# http://stackoverflow.com/questions/11546256/two-way-density-plot-combined-with-one-way-density-plot-with-selected-regions-in
marginal_gating_plot <- function(data, feature_pairs, bins = 30) {
  require('gtable')
  require('ggplot2')
  data <- as.data.frame(data)
  
  # To automatically plot specified features, we use the `ggplot2:::aes_string` function below.
  # This function does not handle hyphenated column names well. To handle this, we apply the
  # `make.names` function to convert to `ggplot2`-friendly names. However,
  fp_orig <- feature_pairs
  colnames(data) <- make.names(colnames(data))
  feature_pairs <- make.names(feature_pairs)

  # Hexbin plot
  p1 <- ggplot(data, aes_string(x = feature_pairs[1], y = feature_pairs[2]))
  p1 <- p1 + stat_binhex(alpha = 1, bins = bins)
  p1 <- p1 + scale_fill_gradientn(colours = c("darkblue", "lightblue", "green", "yellow", "orange", "red"))
  p1 <- p1 + coord_cartesian(range(data[, feature_pairs[1]]), range(data[, feature_pairs[2]]))
  p1 <- p1 + opts(legend.position = "none") + xlab(fp_orig[1]) + ylab(fp_orig[2])

  # Marginal plot for first feature (x-axis)
  p2 <- ggplot(data, aes_string(x = feature_pairs[1]))
  p2 <- p2 + stat_density(fill = "darkblue")
  p2 <- p2 + coord_cartesian(range(data[, feature_pairs[1]]))

  # Marginal plot for second feature (y-axis)
  p3 <- ggplot(data, aes_string(x = feature_pairs[2]))
  p3 <- p3 + stat_density(fill = "darkblue")
  p3 <- p3 + coord_flip(range(data[, feature_pairs[2]]))
  
  gt <- ggplot_gtable(ggplot_build(p1))
  gt2 <- ggplot_gtable(ggplot_build(p2))
  gt3 <- ggplot_gtable(ggplot_build(p3))
  
  
  gt1 <- gtable_add_cols(gt, unit(0.3, "null"), pos = -1)
  gt1 <- gtable_add_rows(gt1, unit(0.3, "null"), pos = 0)
  
  gt1 <- gtable_add_grob(gt1, gt2$grobs[[which(gt2$layout$name == "panel")]],
                         1, 4, 1, 4)
  gt1 <- gtable_add_grob(gt1, gt2$grobs[[which(gt2$layout$name == "axis-l")]],
                         1, 3, 1, 3, clip = "off")
  
  gt1 <- gtable_add_grob(gt1, gt3$grobs[[which(gt3$layout$name == "panel")]],
                         4, 6, 4, 6)

  gt1 <- gtable_add_grob(gt1, gt3$grobs[[which(gt3$layout$name == "axis-b")]],
                         5, 6, 5, 6, clip = "off")
  
  grid.newpage()
  grid.draw(gt1)
  
  invisible(gt1)
}

#' Pregates a flowSet before compensation or transformation
#'
#' We consider several candidate values for the number of clusters to use in the
#' fitted model and select the value that optimizes the \code{criterion}
#' selected. After the best model is selected, we remove clusters having a
#' relatively small number of cells. To do this, we keep the largest clusters so
#' that the total proportion of cells exceeds \code{cells_kept}.
#'
#' @param flow_set a \code{flowSet} object
#' @param cells_kept the proportion of cells that are guaranteed to be kept.
#' See details.
#' @param K the candidate number of clusters to consider via \code{flowClust}
#' @param trans the transformation parameter that is passed to \code{flowClust}
#' @param nu.est the degrees of freedom estimation code that is passed to
#' \code{flowClust}
#' @param criterion the criterion to use for the automatic selection of \code{K}.
#' By default, the \code{ICL} is used
#' @param ... additional arguments passed to \code{flowClust}
pregate_flowset <- function(flow_set, cells_kept = 0.8, K = 2:4, trans = 0,
                            nu.est = 2, criterion = c("ICL", "BIC"),
                            verbose = TRUE, plot = TRUE, ...) {

  criterion <- match.arg(criterion)

  fsApply(flow_set, function(flow_frame) {
    if (verbose) {
      message("Filename: ", keyword(flow_frame)$FIL)
    }

    fc_out <- flowClust(x = flow_frame, varNames = c("FSC-A", "SSC-A"), K = K,
                        trans = trans, nu.est = nu.est, criterion = criterion, ...)
    optimal_K <- K[fc_out@index]
    cluster_proportions <- fc_out[[fc_out@index]]@w

    if (verbose) {
      message("Optimal Value of K: ", optimal_K)
      message("Cluster Proportions:")
      print(cluster_proportions)
    }

    # Here, we compute the cumulative proportions of the clusters fitted in
    # order from the largest cluster to the smallest cluster. We then keep only
    # the largest clusters such that the sum of the proportions exceeds the
    # cutoff specified by the user.
    prop_cumsum <- cumsum(sort(cluster_proportions, decreasing = TRUE))
    num_clusters <- which(prop_cumsum > cells_kept)[1]
    clusters_kept <- order(cluster_proportions, decreasing = TRUE)[seq_len(num_clusters)]

    if (plot) {
      # These two lines are from 'stats:::plot.lm'. They prompt the user for
      # input before proceeding with the next plot.
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
      
      plot(x = K, y = sapply(fc_out, function(x) x@ICL), type = "l",
           xlab = "Candiate Values of K", ylab = "ICL")
      plot(flow_frame, fc_out, main = "flowClust Fit for Optimal Value of K")
    }

    # Extract a flowFrame object containing the largest clusters
    split(x = flow_frame, fc_out, population = list(clusters_kept))[[1]]
  })
}

##' Clean up the T-cell centralized manual gating table from Marc
##'
##' This function will read and clean up the manual gating table from Marc for the T-cells.
##' Unfortunately all the panels are formatted differently and need to be handled separately.
##' 
##' @title prepareManualTcellGates
##' @param manual A data.table of the manual gates for T-cells are passed in, as returned by read.xlsx with check.names=FALSE
##' @return A data frame of just the T-cell statistics
##' @author Greg Finak
prepareManualTcellGates <- function(manual){
                                        #Clean up column names for the T-cells
  Tcell <- manual
  colnames(Tcell)<-(gsub("\\n"," ",colnames(Tcell))) #remove newlines from column headers
  Tcell<-subset(Tcell,!File%in%c("mean","SD","CV")) # remove the estimates of the mean, CV and SD.
  Tcell$File <- factor(Tcell$File)

                                        #What's in the remaining NA columns? The first column is "Center", the rest we can toss out.
  Tcell<-(Tcell[,-grep("NA",colnames(Tcell))[-1L]]) #drop irrelevant columns
  colnames(Tcell)[grep("NA",colnames(Tcell))] <- "Center" #Second column is center
  colnames(Tcell) <- gsub("\\.([0-9])$","\\.CD8",colnames(Tcell)) # rename the populations that seem to have been calculated an alternate way (some other parent)
  Tcell <- Tcell[-1L,] # drop the first row, it is empty

                                        #Fill in blank rows for sample and center
  Tcell[is.na(Tcell[,1]),1]<-""
  Tcell[,1]<-factor(fill(Tcell[,1]))
  Tcell[,2] <- as.character(Tcell[,2])
  Tcell[is.na(Tcell[,2]),2] <- ""
  Tcell[,2] <- factor(fill(Tcell[,2]))
  Tcell$Sample <- factor(as.numeric(as.character(factor(Tcell$Sample))))
                                        #Remove missing data from centers that were excluded
                                        #Tcell<-na.omit(Tcell)   #actually don't. We'll keep the excluded centers.
  Tcell <- rename(melt(Tcell,id=c("Sample","Center","File")),c(variable="population",value="proportion"))
  Tcell$proportion <- as.numeric(as.character(Tcell$proportion))/100
  invisible(Tcell)
}

##' Read the T-helper manual gates statistics
##'
##' Cleans up the data for comparisons
##' @title prepareManualThelperGates
##' @param manual 
##' @return data.frame
##' @author Greg Finak
prepareManualThelperGates <- function(manual){
  Thelper <- manual
  Thelper <- Thelper[,-ncol(Thelper)]
  Thelper <- Thelper[,apply(Thelper,2,function(x)!(all(is.na(x))))]
  Thelper$Sample <- fill(Thelper$Sample,NA)
  cn <- colnames(Thelper)
  cn[2] <- "Center"
  cn <- gsub("\\.1","\\.CD8",gsub("\\n","",cn))
  colnames(Thelper) <- cn
  Thelper$Center <- fill(Thelper$Center,NA)
  Thelper <- subset(Thelper,!File%in%c("mean","SD","CV"))
  Thelper$File <- factor(Thelper$File)
  Thelper$Sample <- factor(Thelper$Sample)
  Thelper <- melt(Thelper,id=c("Sample","Center","File"))
  setnames(Thelper,c("variable","value"),c("population","proportion"))
  Thelper$proportion <- as.numeric(as.character(Thelper$proportion))/100
  invisible(Thelper)
}

##' Parse the manual gating statistics for Treg panel
##'
##' 
##' @title prepareManualTregGates
##' @param manual 
##' @return data.frame
##' @author Greg Finak
prepareManualTregGates <- function(manual){
  Treg <- manual
                                        #Clean up the columns
  cn <- colnames(Treg)
  head <- Treg[1:2,]
  cn <- fill(cn,"NA")
  cn <- gsub(" $","",gsub("^ ","",gsub("\\n","",gsub("NA","",paste(cn,apply(head,2,function(x)paste(x,collapse=" ")))))))
  Treg <- Treg[-c(1:2),]
  colnames(Treg) <- cn
  Treg <- Treg[,-ncol(Treg)]
  colnames(Treg)[2] <- "Center"
  Treg <- Treg[,!apply(Treg,2,function(x)all(is.na(x)))]
  Treg$Sample <- fill(Treg$Sample,NA)
  Treg$Center <- fill(Treg$Center,NA)
  Treg <- subset(Treg,!File%in%c("mean","SD","CV"))
  Treg$File <- factor(Treg$File)
  Treg$Sample <- factor(as.numeric(as.character(factor(Treg$Sample))))
  Treg <- melt(Treg,id=c("Sample","Center","File"))
  setnames(Treg,c("variable","value"),c("population","proportion"))
  Treg$proportion <- as.numeric(as.character(Treg$proportion))/100
  invisible(Treg)
}



##' Clean up the B-cell gates
##'
##' Reads and cleans up the manual gates from Marc.
##' @title prepareManualBcellGates
##' @param manual a data.frame read in by read.xlsx with check.names=FALSE
##' @return data.frame of the B-cell statistics.
##' @author Greg Finak
prepareManualBcellGates <- function(manual){
  Bcell <- manual
                                        #Parse rows one and two as a new header
  h1 <- colnames(Bcell)
  h2 <- Bcell[1,]

  Bcell <- Bcell[-1L,]#Drop the first row, which  we have stored.
  colnames(Bcell) <- gsub("\\n"," ",as.matrix(h2)) #Remove newlines from column names
  colnames(Bcell)[2] <- c("Center") #name the column which contains the center
  Bcell <- subset(Bcell,!File%in%c("mean","CV","SD")) #Drop the mean, SD and CV
  Bcell$File <- factor(Bcell$File)
  tmp <- as.character(Bcell[,1]) #Fill in missing entries
  tmp[is.na(tmp)] <- ""
  Bcell[,1] <- factor(fill(tmp))
  tmp <- as.character(Bcell[,2])
  tmp[is.na(tmp)] <- ""
  Bcell[,2] <- factor(fill(tmp))

                                        #drop NA.3
  Bcell <- Bcell[,-grep("NA\\.3",colnames(Bcell))]
                                        #NA.1 and NA.2 separate the populatoins that are parented by h1
                                        #prefix the populations with the appropriate parent name
  cur<-eval(as.call(c(`seq`,as.list(grep("NA",colnames(Bcell))+c(1,-1)))))
  colnames(Bcell)[cur] <- paste(h1[11],colnames(Bcell)[cur],sep="/")
  colnames(Bcell)[16:18] <- paste(h1[16],colnames(Bcell)[16:18],sep="/")
                                        #remove NA columns we don't need
  Bcell <- Bcell[,-grep("NA",colnames(Bcell))]
  Bcell <- Bcell[,-4]
                                        # Bcell <- na.omit(Bcell) #Don't drop NA columns, we'll keep factors for the missing centers.
  Bcell$factor <- factor(as.numeric(as.character(factor(Bcell$Sample))))
  Bcell <- rename(melt(Bcell,id=c("Sample","Center","File")),c(variable="population",value="proportion"))
  Bcell$proportion <- as.numeric(as.character(Bcell$proportion))/100
  invisible(Bcell)
}

##' Fill blank elements in intervals of a vector up to the next populated cells
##'
##' @name fill
##' @title fill
##' @export fill
fill <- function(x,y=''){
  wh <- which(!x%in%y)
  wh <- c(wh,length(x)+1)
  for(i in 1:((length(wh)-1))){
    x[wh[i]:(wh[i+1]-1)] <- x[wh[i]]
  }
  x
}

quantile_negatives <- function(x, ...) {
  x <- x[x < 0]
  quantile(x, ...)
}


#' Creates hexbin objects of FCSTrans transformations for two markers of interest
#' based on the transformations widths given
hexbin_transformations <- function(center, widths, marker_pairs) {

  temp_flowset <- fs_list[[center]]

  first_markers <- unique(sapply(marker_pairs, head, n = 1))
  second_markers <- unique(sapply(marker_pairs, tail, n = 1))

  # Grabs channels for first and second markers
  pData_flowset <- na.omit(pData(parameters(temp_flowset[[1]])))
  pData_flowset$desc <- marker_conversion(pData_flowset$desc)
  first_channels <- with(pData_flowset, name[match(first_markers, desc)])
  second_channels <- with(pData_flowset, name[match(second_markers, desc)])
  
  # Transforms first_channels by widths[1]
  trans_channel1 <- lapply(first_channels, function(channel) {
    FCSTransTransform(w = widths[1])
  })
  trans_channel1 <- transformList(first_channels, trans_channel1)
  temp_flowset <- transform(temp_flowset, trans_channel1)

  # Transforms second_channels by widths[2]
  trans_channel2 <- lapply(second_channels, function(channel) {
    FCSTransTransform(w = widths[2])
  })
  trans_channel2 <- transformList(second_channels, trans_channel2)
  temp_flowset <- transform(temp_flowset, trans_channel2)

  # For each marker pair, creates hexbin for each sample within the flowSet
  hexbin_results <- lapply(marker_pairs, function(marker_pair) {
    channel_pair <- with(pData_flowset, name[match(marker_pair, desc)])

    lapply(seq_along(temp_flowset), function(j) {
      x <- exprs(temp_flowset[[j]])
      x <- x[x[, "SSC-A"] > 0, ]
      hexbin(x[, channel_pair])
    })
  })
  # Naming convention, e.g., IgD_CD27
  names(hexbin_results) <- sapply(marker_pairs, paste, collapse = "_")
  
  hexbin_results
}

#' Converts a marker quantile to the FCSTrans width
quantile2width <- function(quantile) {
  abs(0.5 * (4.5 * log(10) - log(2^18 / abs(quantile))))
}

##'
##' Reads the Monocyte Centralized gating, cleaning up the headers.
##'
##' 
##' @title prepareManualMonocyteGates
##' @param manual 
##' @return A melted data frame with the monocyte stats 
##' @author Greg Finak
prepareManualMonocyteGates <- function(manual){
Mono <- manual
                                        #First three rows are header
  cn<-colnames(Mono)
                                        #clean up the header
  head<-gsub(" $","",gsub("^ ","",gsub("\\n","",gsub("NA","",apply(Mono[1:3,],2,function(x)paste0(as.character(x),collapse=" "))))))
  cn<-gsub(" $","",gsub("^ ","",paste(fill(gsub("NA","",names(head))),head,sep=" ")))
  cn[2] <- "Center"
  cn[cn==""] <- NA
  Mono<-Mono[-c(1:3),]
  colnames(Mono) <- cn
                                        #Fill empy rows
  Mono$File <- fill(Mono$File,NA)
  Mono$Sample <- fill(Mono$Sample,NA)
  Mono$Center <- fill(Mono$Center,NA)
                                        #Drop NA columns
  Mono<-Mono[,!(is.na(colnames(Mono))|c(rep(FALSE,length(cn)-1),TRUE))]
                                        #Remove Mean, SD,CV
  Mono <- subset(Mono,!File%in%c("mean","SD","CV"))
  Mono$File <- factor(Mono$File)
  Mono$Sample <- factor(as.numeric(as.character(factor(Mono$Sample))))
  Mono<-melt(Mono,id=c("Sample","Center","File"))
  setnames(Mono,c("variable","value"),c("population","proportion"))
  Mono$proportion <- as.numeric(as.character(Mono$proportion))/100
  Mono
}

##' Returns the levels of a column for a list of data frames
##'
##' 
##' @title get the levels of a specific column
##' @param A list of data frames
##' @param col name of the column 
##' @return 
##' @author Greg Finak
getLevels <- function(A,col){
  lapply(A,function(x)if(!is.null(x)){levels(get(col,x))})
}
