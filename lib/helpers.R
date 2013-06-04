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
  
  # Preprocesses each of the columns in the flow_frame
  for(j in seq_len(ncol(flow_frame))) {
    marker_idx <- paste0("$P", j, "S")
    channel_idx <- paste0("$P", j, "N")

    marker <- flow_frame@description[[marker_idx]]
    channel <- flow_frame@description[[channel_idx]]

    # In the case the marker name is given, we swap the marker and channel
    # names.
    if (!is.null(marker)) {
      # For the following list of marker names, we manually update the names so
      # that they are standard across centers.
      if (marker == "19") {
        marker <- "CD19"
      } else if (marker == "LIVE_GREEN") {
        marker <- "LIVE GREEN"
      } else if (marker == "IGD") {
        marker <- "IgD"
      }

      # If marker name contains additional info, remove everything after the
      # space. (e.g., "IgD V500" to "IgD")
      marker <- strsplit(marker, " ")[[1]][1]

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
                                  method = c("median", "mean", "mode"), ...) {
                                  
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

#  marginal_plots <- lapply(seq_along(comp_flowset), function(i) {
#    message("Compensation Control: ", i)
#    marginal_gating_plot(data = exprs(comp_flowset[[i]]), feature_pairs = c("FSC-A", "SSC-A"))
#  })

  # Applies flowClust to SSC-A and FSC-A to remove outliers and debris before
  # constructing a compensation matrix. Keeps the densest clusters and gates out
  # the rest of the cells.
  if (pregate) {
    comp_flowset <- pregate_flowset(comp_flowset, plot = plot, ...)
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
            useNormFilt = FALSE, method = method, regexpr = TRUE)
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
#' @param min_limit the minimum limit to be passed to \code{read.flowSet}
#' @return a \code{flowSet} object
flowset_lyoplate <- function(path, xlsx, comp_matrix, transform = TRUE, center,
                             panel = c("B cell", "T cell", "Treg", "DC/mono/NK", "Th1/2/17"),
                             min_limit = -100) {
  panel <- match.arg(panel)

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
  exp_flowset <- read.flowSet(fcs_files, min.limit = min_limit)

  # Updates the flowSet's pData
  exp_pdata <- pData(exp_flowset)
  exp_pdata <- merge(exp_pdata, exp_samples)
  rownames(exp_pdata) <- rownames(pData(exp_flowset))
  pData(exp_flowset) <- exp_pdata 
  varM <- varMetadata(phenoData(exp_flowset))
  varM[-1,] <- rownames(varM)[-1]
  varMetadata(phenoData(exp_flowset)) <- varM

  # Compensates and transforms flowSet for current center
  exp_flowset <- compensate(exp_flowset, comp_matrix)

  # Applies the FCStrans transformation to the experimental samples.
  # The channels transformed correspond to the column names of the compensation
  # matrix constructed from the compensation controls. However, note that
  # FCStrans compensates the samples using the spillover matrix stored in the
  # keywords.
  if (transform){
    trans <- transformList(from = colnames(comp_matrix), tfun = FCSTransTransform())
    exp_flowset <- transform(exp_flowset, trans)
  }
  
  exp_flowset
}

##################################################
## B subpopulations
## Gates IgD vs CD27
##################################################
gating.bsub <- function(gs_bcell, parent = "cd19 & cd20") {

  parent_data <- getData(gs_bcell, parent)

  xChannel <- "IgD"
  yChannel <- "CD27"

  prior_CD27 <- list(CD27 = openCyto:::prior_flowClust(parent_data, "CD27", nu0 = 30, min = 0.01))

  fs_list <- split(parent_data, sampleNames(parent_data))

  ## Apply the CD27 gate first.
  cd27.filterList <- mclapply(fs_list, function(fs_sample) {
    openCyto:::flowClust.1d(fs_sample, yChannel = yChannel, prior = prior_CD27, usePrior = "yes", nu = 30, min = 0.01)
  }, mc.cores = 12)


  # Then gate IgD on the cd27+
  cd27_positive <- flowSet(lapply(sampleNames(parent_data), function(curSample) {
    split(parent_data[[curSample]], cd27.filterList[[curSample]])[[1]]
  }))
  sampleNames(cd27_positive) <- sampleNames(parent_data)

  # Here, we collapse all samples. Then, we take the prior of the negative
  # population to be the largest peak after applying a lot of smoothing.
  # Generally, there is only one visible peak anyways.
  IgD_collapsed <- exprs(as(cd27_positive, "flowFrame"))[, xChannel]
  IgD_negative_peak <- openCyto:::find_peaks(IgD_collapsed, adjust = 3)

  # We use three components to model the IgD marker. We use two components for
  # the negative peak and one component for the positive peak.
  prior_IgD <- list()
  huber_s <- huber(IgD_collapsed)$s
  prior_IgD$Mu0 <- matrix(c(IgD_negative_peak - 3.5 * huber_s, IgD_negative_peak, IgD_negative_peak + 2.75 * huber_s), nrow = 3)
  prior_IgD$Omega0 <- array(rep(0.05^2, 3), dim = c(3, 1, 1))
  prior_IgD$Lambda0 <- array(rep(0.5^2, 3), dim = c(3, 1, 1))
  prior_IgD$w0 <- rep(30, 3)
  prior_IgD$nu0 <- rep(30, 3)

  fs_list <- split(cd27_positive, sampleNames(cd27_positive))

  IgD.filterList <- mclapply(fs_list, function(fs_sample) {
    openCyto:::.flowClust.1d(fs_sample, yChannel = xChannel, prior = list(IgD = prior_IgD),
                             usePrior = "yes", nu = 30, K = 3, cutpoint_method = "quantile", quantile = 0.95)
  }, mc.cores = 10)

  # construct quadGate afterwards
  bsub.filterList <- sapply(names(IgD.filterList), function(curSample) {
    xfilter <- IgD.filterList[[curSample]]
    
    # change from neg to pos for quadGate only recognize (xx,Inf)
    yfilter <- cd27.filterList[[curSample]]
    coord <- list(unname(xfilter@min), unname(yfilter@min))
    names(coord) <- c(parameters(xfilter), parameters(yfilter))
    quadGate(coord)
  })

  nodeID <- add(gs_bcell, filterList(bsub.filterList), parent = parent)
  recompute(gs_bcell, nodeID)
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
pregate_flowset <- function(flow_set, cells_kept = 0.7, K = 2:6, trans = 0,
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
