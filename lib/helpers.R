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

    if (is.null(marker)) {
      # In the case that the marker name is not given, we set the marker and channel
      # names both to that given in the channel.
      marker <- channel

      # Updates the marker information in the flow_frame with the channel
      flow_frame@description[[marker_idx]] <- channel
      flow_frame@parameters@data$desc[j] <- channel
    } else {
      # If marker name contains additional info, remove everything after the space. (e.g., "IgD V500" to "IgD")
      marker <- strsplit(marker, " ")[[1]][1]

      # In the special case of Yale, they used "19" instead of "CD19". In this case,
      # we manually update the marker name.
      if (marker == "19") {
        marker <- "CD19"
      }

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
#' @param path the path that contains the compensation control FCS files
#' @param xlsx the name of the XLSX (Excel) file (full path)
#' @param panel the panel type
#' @return the compensation matrix
compensation_lyoplate <- function(path, xlsx, panel = c("Bcell", "Tcell")) {
  require('xlsx')

  panel <- match.arg(panel)
  comp_controls <- read.xlsx2(file = xlsx, sheetName = "Comp controls",
                              startRow = 3, stringsAsFactors = FALSE)

  # First, we remove the spaces in the tube names
  comp_controls$Applies.to.tube <- gsub(" ", "", comp_controls$Applies.to.tube)
  tubes_split <- strsplit(comp_controls$Applies.to.tube, ",")

  if (panel == "Bcell") {
    which_files <- sapply(tubes_split, function(x) any(x %in% c("all", "B", "Bcell")))
  } else if (panel == "Tcell") {
    which_files <- sapply(tubes_split, function(x) any(x %in% c("all", "T", "Tcell")))
  }

  # We exclude the entries that have empty filenames.
  which_files <- which_files & comp_controls$FCS.file.name != ""

  # Constructs a flowSet of the compensation control FCS files along with the
  # unstained control FCS file
  FCS_files <- comp_controls$FCS.file.name[which_files]
  unstained_control <- dir(path, pattern = "Compensation.*Unstained")
  FCS_files <- file.path(path, c(FCS_files, unstained_control))

  comp_flowset <- read.flowSet(FCS_files)

  # Constructs the compensation (spillover) matrix
  # The index of 'unstained' is the last of the FCS files read in.
  # A regular expression in 'patt' indicates the markers that we want to compensate.
  #    These are given in the XLSX file.
  spillover(comp_flowset, unstained = length(comp_flowset),
            patt = paste(comp_controls$Marker[which_files], collapse = "|"),
            useNormFilt = TRUE, method = "mean")
}

#' Constructs a GatingSet for the specified center for Lyoplate 3.0
#'
#' Creates a \code{flowSet} by reading the FCS files for the current center given
#' in the Excel file. Then, compensates and transforms them. Finally, constructs
#' and returns a \code{GatingSet} from the \code{flowSet} object.
#'
#' @param path the path that contains the sample FCS files
#' @param xlsx the name of the XLSX (Excel) file (full path)
#' @param comp_matrix a compensation matrix
#' @param center a character string denoting the current center
#' @param panel the panel type
#' @return a \code{GatingSet} object
gatingset_lyoplate <- function(path, xlsx, comp_matrix, center,
                               panel = c("Bcell", "Tcell", "Treg")) {
  require('xlsx')
  panel <- match.arg(panel)

  exp_samples <- read.xlsx2(file = xlsx, sheetName = "Exp samples",
                            startRow = 3, stringsAsFactors = FALSE)
  colnames(exp_samples) <- c("name", "Institution", "Panel", "Replicate",
                             "Sample")
  exp_samples <- subset(exp_samples, select = -Institution)
  exp_samples$Replicate <- as.character(as.numeric(exp_samples$Replicate))
  exp_samples$Sample <- as.character(as.numeric(exp_samples$Sample))

  # Removes whitespace from Panel names
  exp_samples$Panel <- gsub(" ", "", exp_samples$Panel)

  # Subsets the data for the specified panel
  exp_samples <- exp_samples[exp_samples$Panel == panel, ]

  # Reads in the FCS files for the specified panel
  fcs_files <- file.path(path, exp_samples$name)
  exp_flowset <- read.flowSet(fcs_files, min.limit = -100)

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
  trans <- openCyto:::estimateMedianLogicle(exp_flowset,
                                            channels = colnames(comp_matrix))
  exp_flowset <- transform(exp_flowset, trans)

  GatingSet(exp_flowset)
}




##################################################
## B subpopulations
## Gates IgD vs CD27
##################################################
gating.bsub <- function(gs_bcell, parent = "cd19 & cd20") {

  parent_data <- getData(gs_bcell, parent)

  xChannel <- "IgD"
  yChannel <- "CD27"

  prior_CD27 <- list(CD27 = openCyto:::prior_flowClust(parent_data, "CD27", nu0 = 30))

  fs_list <- split(parent_data, sampleNames(parent_data))

  ## Apply the CD27 gate first.
  cd27.filterList <- mclapply(fs_list, function(fs_sample) {
    openCyto:::.flowClust.1d(fs_sample, yChannel = yChannel, prior = prior_CD27, usePrior = "yes", nu = 30)
  }, mc.cores = 10)


  # Then gate IgD on the cd27+
  cd27_positive <- flowSet(lapply(sampleNames(parent_data), function(curSample) {
    split(parent_data[[curSample]], cd27.filterList[[curSample]])[[1]]
  }))
  sampleNames(cd27_positive) <- sampleNames(parent_data)

  # Here, we collapse all samples. Then, we take the prior of the negative
  # population to be the largest peak after applying a lot of smoothing.
  # Generally, there is only one visible peak anyways.
  IgD_collapsed <- exprs(as(cd27_positive, "flowFrame"))[, xChannel]
  IgD_negative_peak <- openCyto:::find_peaks(IgD_collapsed, peaks = 1, adjust = 3)

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
