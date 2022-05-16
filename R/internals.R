#' @title Rarefy ampvis2 object
#' @description This is just a wrapper of \code{\link[vegan]{rrarefy}} with convenient error messages and adjusted to work with ampvis2 objects.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param rarefy (\emph{required}) Passed directly to \code{\link[vegan]{rrarefy}}.
#'
#' @return An ampvis2 object with rarefied OTU abundances.
#' @importFrom vegan rrarefy
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @keywords internal
amp_rarefy <- function(data, rarefy) {
  ### Data must be in ampvis2 format
  if (class(data) != "ampvis2") {
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis2 functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)", call. = FALSE)
  }

  set.seed(0) # for reproducibility

  reads <- colSums(data$abund)
  maxreads <- max(reads)
  minreads <- min(reads)
  if (is.numeric(rarefy)) {
    if (rarefy > maxreads) {
      stop("The chosen rarefy size is larger than the largest amount of reads in any sample (", as.character(maxreads), ").", call. = FALSE)
    } else if (rarefy < minreads) {
      data$abund <- suppressWarnings(vegan::rrarefy(t(data$abund), sample = rarefy)) %>%
        t() %>%
        as.data.frame()
      warning("The chosen rarefy size (", as.character(rarefy), ") is smaller than the smallest amount of reads in any sample (", as.character(minreads), ").", call. = FALSE)
    } else {
      data$abund <- suppressWarnings(vegan::rrarefy(t(data$abund), sample = rarefy)) %>%
        t() %>%
        as.data.frame()
      if (minreads < rarefy) {
        message("The following sample(s) have not been rarefied (less than ", as.character(rarefy), " reads):\n", paste(rownames(data$metadata[which(reads < rarefy), ]), collapse = ", "))
      }
    }
  } else if (!is.numeric(rarefy)) {
    stop("Argument rarefy must be numerical.", call. = FALSE)
  }
  return(data)
}

#' @title Tidy up taxonomy
#' @description Used internally in other ampvis functions.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible.
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param tax_level The taxonomic level to remove OTUs with no assigned taxonomy, only used when \code{tax_empty = "remove"}. (\emph{default:} \code{"Genus"})
#'
#' @return A list with 3 dataframes (4 if reference sequences are provided).
#' @importFrom dplyr mutate
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' @keywords internal
amp_rename <- function(data, tax_class = NULL, tax_empty = "best", tax_level = "Genus") {
  tax <- data[["tax"]]

  ## First make sure that all entries are strings
  tax[] <- lapply(tax, as.character)

  ## Change a specific phylum to class level
  if (!is.null(tax_class)) {
    for (i in 1:nrow(tax)) {
      if (!is.na(tax$Phylum[i]) & tax$Phylum[i] %in% tax_class) {
        tax$Phylum[i] <- tax$Class[i]
      }
    }
  }

  ## Remove QIIME taxonomy format prefixes across all levels
  tax[] <- lapply(tax, gsub, pattern = "^[dkpcofgs]_+", replacement = "")

  ## Check if there is a species level otherwise add it for consistency
  if (is.null(tax$Species)) {
    tax$Species <- ""
  }

  ## NA's as empty strings
  tax[is.na(tax)] <- ""

  ## How to handle empty taxonomic assignments
  if (tax_empty == "OTU") {
    for (i in 1:nrow(tax)) {
      if (tax[i, "Species"] == "") {
        tax[i, "Species"] <- rownames(tax)[i]
      }
      if (tax[i, "Genus"] == "") {
        tax[i, "Genus"] <- rownames(tax)[i]
      }
      if (tax[i, "Family"] == "") {
        tax[i, "Family"] <- rownames(tax)[i]
      }
      if (tax[i, "Order"] == "") {
        tax[i, "Order"] <- rownames(tax)[i]
      }
      if (tax[i, "Class"] == "") {
        tax[i, "Class"] <- rownames(tax)[i]
      }
      if (tax[i, "Phylum"] == "") {
        tax[i, "Phylum"] <- rownames(tax)[i]
      }
    }
  }

  ## Handle empty taxonomic strings
  rn <- rownames(tax) # damn rownames are silently dropped by mutate()
  if (tax_empty == "best") {
    tax <- mutate(tax, Kingdom, Kingdom = ifelse(Kingdom == "", "Unclassified", Kingdom)) %>%
      mutate(Phylum, Phylum = ifelse(Phylum == "", paste("k__", Kingdom, "_", rownames(tax), sep = ""), Phylum)) %>%
      mutate(Class, Class = ifelse(Class == "", ifelse(grepl("__", Phylum), Phylum, paste("p__", Phylum, "_", rownames(tax), sep = "")), Class)) %>%
      mutate(Order, Order = ifelse(Order == "", ifelse(grepl("__", Class), Class, paste("c__", Class, "_", rownames(tax), sep = "")), Order)) %>%
      mutate(Family, Family = ifelse(Family == "", ifelse(grepl("__", Order), Order, paste("o__", Order, "_", rownames(tax), sep = "")), Family)) %>%
      mutate(Genus, Genus = ifelse(Genus == "", ifelse(grepl("__", Family), Family, paste("f__", Family, "_", rownames(tax), sep = "")), Genus)) %>%
      mutate(Species, Species = ifelse(Species == "", ifelse(grepl("__", Genus), Genus, paste("g__", Genus, "_", rownames(tax), sep = "")), Species))
  }
  rownames(tax) <- rn

  if (tax_empty == "remove") {
    abund <- data[["abund"]]
    tax <- subset(tax, tax[, tax_level] != "")
    abund <- subset(abund, rownames(abund) %in% rownames(tax))
    data[["abund"]] <- abund
  }
  data[["tax"]] <- tax
  rownames(data[["tax"]]) <- rownames(tax)

  return(data)
}

#' @title Get data from the MiDAS field guide API
#' @description Gets all fields from the MiDAS field guide (\url{https://midasfieldguide.org}) returned in a list.
#'
#' @return A list with all fields.
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @keywords internal
getMiDASFGData <- function() {
  checkReqPkg("jsonlite")
  jsonlite::read_json("http://midasfieldguide.org/api/microbes/fieldguide")
}

#' @title Extract functional information about Genera from the MiDAS field guide
#' @description Extract field related to properties and metabolism of all Genera from a list obtained by \code{\link{getMiDASFGData}} and return in a data frame.
#'
#' @param FGList Data list obtained by \code{\link{getMiDASFGData}}.
#'
#' @importFrom data.table rbindlist
#'
#' @return A data frame where each row is a Genus and each column is a "function".
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @keywords internal
extractFunctions <- function(FGList) {
  functions <- lapply(FGList, function(x) {
    outList <- lapply(c(x[["properties"]], x[["metabolism"]]), function(y) {
      if (y[["In situ"]] == "Positive" |
        (any(y[["In situ"]] %in% c("Variable", "Not Assessed")) &
          y[["Other"]] == "Positive")) {
        return("POS")
      }
      if (y[["In situ"]] == "Negative" |
        (y[["In situ"]] == "Not Assessed" & y[["Other"]] == "Negative")) {
        return("NEG")
      }
      if (y[["In situ"]] == "Variable" |
        (y[["In situ"]] == "Not Assessed" & y[["Other"]] == "Variable")) {
        return("VAR")
      }
      if (all(c(y[["In situ"]], y[["Other"]]) %in% "Not Assessed")) {
        return("NT")
      }
    })
    c(
      "Genus" = gsub("^Ca ", "Ca_", x[["name"]]),
      "MiDAS" = "POS",
      outList
    )
  })
  outDF <- as.data.frame(data.table::rbindlist(functions, fill = TRUE))
  outDF[is.na(outDF)] <- "NT"
  return(outDF)
}

#' @title Calculate weighted or unweighted UniFrac distances. Adopted from fastUniFrac() from phyloseq
#'
#' @param abund Abundance table with OTU counts, in \code{ampvis2} objects it is available with simply data$abund
#' @param tree Phylogenetic tree (rooted and with branch lengths) as loaded with \code{\link[ape]{read.tree}}.
#' @param weighted Calculate weighted or unweighted UniFrac distances.
#' @param normalise Should the output be normalised such that values range from 0 to 1 independent of branch length values? Note that unweighted UniFrac is always normalised. (\emph{default:} \code{TRUE})
#' @param num_threads The number of threads to be used for calculating UniFrac distances. If set to more than \code{1} then this is set by using \code{\link[doParallel]{registerDoParallel}}. (\emph{default:} \code{1})
#'
#' @importFrom ape is.rooted node.depth node.depth.edgelength reorder.phylo
#' @return A distance matrix of class \code{dist}.
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @keywords internal
unifrac <- function(abund,
                    tree,
                    weighted = FALSE,
                    normalise = TRUE,
                    num_threads = 1L) {
  if (num_threads == 1L) {
    checkReqPkg("foreach")
  } else if (num_threads > 1L) {
    checkReqPkg("doParallel")
  }

  # check tree
  if (!class(tree) == "phylo") {
    stop("The provided phylogenetic tree must be of class \"phylo\" as loaded with the ape::read.tree() function.", call. = FALSE)
  }
  if (!ape::is.rooted(tree)) {
    stop("Tree is not rooted!", call. = FALSE)
    # message("Tree is not rooted, performing a midpoint root")
    # tree <- phytools::midpoint.root(tree)
  }
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute UniFrac", call. = FALSE)
  }
  OTU <- as.matrix(abund)
  ntip <- length(tree$tip.label)
  if (ntip != nrow(OTU)) {
    stop("OTU table and phylogenetic tree do not match", call. = FALSE)
  }
  # if(!all(rownames(OTU) == tree$tip.label))
  #  OTU <- OTU[tree$tip.label, , drop = FALSE]
  node.desc <- matrix(tree$edge[order(tree$edge[, 1]), ][, 2], byrow = TRUE, ncol = 2)
  edge_array <- matrix(0,
    nrow = ntip + tree$Nnode,
    ncol = ncol(OTU),
    dimnames = list(NULL, sample_names = colnames(OTU))
  )
  edge_array[1:ntip, ] <- OTU
  ord.node <- order(ape::node.depth(tree))[(ntip + 1):(ntip + tree$Nnode)]
  for (i in ord.node) {
    edge_array[i, ] <- colSums(edge_array[node.desc[i - ntip, ], , drop = FALSE], na.rm = TRUE)
  }
  edge_array <- edge_array[tree$edge[, 2], ]
  rm(node.desc)
  if (isFALSE(weighted)) {
    edge_occ <- (edge_array > 0) - 0
  }
  if (isTRUE(weighted) & isTRUE(normalise)) {
    z <- ape::reorder.phylo(tree, order = "postorder")
    tipAges <- ape::node.depth.edgelength(tree)
    tipAges <- tipAges[1:length(tree$tip.label)]
    names(tipAges) <- z$tip.label
    tipAges <- tipAges[rownames(OTU)]
  }
  samplesums <- colSums(OTU)
  if (num_threads == 1L) {
    foreach::registerDoSEQ()
  } else if (num_threads > 1) {
    doParallel::registerDoParallel(num_threads)
  }
  spn <- combn(colnames(OTU), 2, simplify = FALSE)
  distlist <- foreach::foreach(i = spn) %dopar% {
    A <- i[1]
    B <- i[2]
    AT <- samplesums[A]
    BT <- samplesums[B]
    if (isTRUE(weighted)) {
      wUF_branchweight <- abs(edge_array[, A] / AT - edge_array[, B] / BT)
      numerator <- sum(
        {
          tree$edge.length * wUF_branchweight
        },
        na.rm = TRUE
      )
      if (isFALSE(normalise)) {
        return(numerator)
      } else {
        denominator <- sum(
          {
            tipAges * (OTU[, A] / AT + OTU[, B] / BT)
          },
          na.rm = TRUE
        )
        return(numerator / denominator)
      }
    } else {
      edge_occ_AB <- edge_occ[, c(A, B)]
      edge_uni_AB_sum <- sum((tree$edge.length * edge_occ_AB)[rowSums(edge_occ_AB, na.rm = TRUE) < 2, ], na.rm = TRUE)
      uwUFpairdist <- edge_uni_AB_sum / sum(tree$edge.length[rowSums(edge_occ_AB, na.rm = TRUE) > 0])
      return(uwUFpairdist)
    }
  }
  UniFracMat <- matrix(NA_real_, ncol(OTU), ncol(OTU))
  rownames(UniFracMat) <- colnames(UniFracMat) <- colnames(OTU)
  matIndices <- do.call(rbind, spn)[, 2:1]
  if (!is.matrix(matIndices)) matIndices <- matrix(matIndices, ncol = 2)
  UniFracMat[matIndices] <- unlist(distlist)
  return(as.dist(UniFracMat))
}

#' Calculate Jensen-Shannon Divergence distances
#'
#' @param abund Abundance table with OTU counts, in \code{ampvis2} objects it is available with simply data$abund.
#' @param pseudocount Pseudocount to use instead of 0-divisions. (\emph{default:} \code{0.000001})
#'
#' @return A distance matrix of class \code{dist}.
# This is based on http://enterotype.embl.de/enterotypes.html
# Abundances of 0 will be set to the pseudocount value to avoid 0-value denominators
# Unfortunately this code is SLOOOOOOOOW
#' @keywords internal
dist.JSD <- function(abund, pseudocount = 0.000001) {
  inMatrix <- t(abund)
  KLD <- function(x, y) sum(x * log(x / y))
  JSD <- function(x, y) sqrt(0.5 * KLD(x, (x + y) / 2) + 0.5 * KLD(y, (x + y) / 2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)

  inMatrix <- apply(inMatrix, 1:2, function(x) ifelse(x == 0, pseudocount, x))

  for (i in 1:matrixColSize) {
    for (j in 1:matrixColSize) {
      resultsMatrix[i, j] <- JSD(
        as.vector(inMatrix[, i]),
        as.vector(inMatrix[, j])
      )
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix) -> resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix)
}

#' @title Find lowest taxonomic level
#' @description Finds the lowest taxonomic level of two given levels in \code{tax}. The hierarchic order of taxonomic levels is simply taken from the order of column names in \code{tax}
#'
#' @param tax_aggregate A character vector of one or more taxonomic levels (exactly as in the column names of \code{ampvis2obj$tax}), fx Genus or Species. (\emph{default:} \code{NULL})
#' @param tax_add A second character vector similar to \code{tax_aggregate}. (\emph{default:} \code{NULL})
#' @param tax The taxonomy table from an ampvis2 object (\code{ampvis2obj$tax})
#'
#' @return A length one character vector with the lowest taxonomic level.
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @keywords internal
getLowestTaxLvl <- function(tax, tax_aggregate = NULL, tax_add = NULL) {
  if (is.null(tax_aggregate) & is.null(tax_add)) {
    tax_aggregate <- colnames(tax)[ncol(tax)]
  }
  # find the lowest taxonomic level of tax_aggregate and tax_add
  taxlevels <- factor(
    x = colnames(tax),
    levels = colnames(tax)
  )
  lowestlevel <- as.character(taxlevels[max(as.numeric(c(
    taxlevels[which(taxlevels %in% tax_aggregate)],
    taxlevels[which(taxlevels %in% tax_add)]
  )))])
  return(lowestlevel)
}

#' @title Aggregate OTUs to a specific taxonomic level
#' @description Calculates the sum of OTUs per taxonomic level
#'
#' @param abund The OTU abundance table from an ampvis2 object (\code{ampvis2obj$abund})
#' @param tax The OTU abundance table from an ampvis2 object (\code{ampvis2obj$tax})
#' @param tax_aggregate Aggregate (sum) OTU's to a specific taxonomic level. (\emph{default:} \code{"OTU"})
#' @param tax_add Add additional (higher) taxonomic levels to the taxonomy string. The OTU's will be aggregated to whichever level of the \code{tax_aggregate} and \code{tax_add} vectors is the lowest. (\emph{default:} \code{NULL})
#' @param calcSums Whether to include the sums of read counts for each sample and taxonomic group. (\emph{default:} \code{TRUE})
#' @param format Output format, \code{"long"} or \code{"abund"}. \code{"abund"} corresponds to that of a read counts table with samples as columns and the aggregated taxa as rows.
#'
#' @importFrom data.table data.table melt
#' @return A data.table.
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @keywords internal
aggregate_abund <- function(abund,
                            tax,
                            tax_aggregate = "OTU",
                            tax_add = NULL,
                            calcSums = TRUE,
                            format = "long") {
  if (any(colnames(abund) %in% "Display")) {
    stop("A column is named \"Display\" in the OTU abundance table, please change it to continue", call. = FALSE)
  }

  # make sure all tax columns are of type character (nchar() does not allow factors)
  tax[] <- lapply(tax, as.character)

  # find the lowest taxonomic level of tax_aggregate and tax_add
  lowestTaxLevel <- getLowestTaxLvl(
    tax = tax,
    tax_aggregate = tax_aggregate,
    tax_add = tax_add
  )

  # Remove all OTUs that are not assigned at the chosen taxonomic level
  # and print a status message with the number of removed OTUs
  newtax <- tax[which(nchar(tax[[lowestTaxLevel]]) > 1 &
    !is.na(tax[[lowestTaxLevel]]) &
    !grepl("^\\b[dkpcofgs]*[_:;]*\\b$", tax[[lowestTaxLevel]])), ]
  newabund <- abund[rownames(newtax), , drop = FALSE]
  if (nrow(newtax) != nrow(tax)) {
    message(paste0(
      nrow(tax) - nrow(newtax),
      " OTUs (out of ",
      nrow(tax),
      ") with no assigned taxonomy at ",
      lowestTaxLevel,
      " level were removed before aggregating OTUs"
    ))
  }

  abundTax <- data.table::data.table(
    newabund,
    Display = apply(
      newtax[, c(tax_add, tax_aggregate), drop = FALSE],
      1,
      paste,
      collapse = "; "
    )
  )
  abundAggr <- data.table::melt(
    abundTax,
    id.vars = "Display",
    variable.name = "Sample",
    value.name = "Abundance",
    variable.factor = FALSE
  )
  if (isTRUE(calcSums)) {
    abundAggr[,
      Sum := sum(Abundance),
      keyby = .(Display, Sample)
    ]
  }
  if (format == "long") {
    out <- abundAggr
  } else if (format == "abund") {
    out <- as.data.frame(
      data.table::dcast(abundAggr,
        Display ~ Sample,
        value.var = "Abundance",
        fun.aggregate = sum
      )
    )
    rownames(out) <- as.character(out[[1]])
    out <- out[, -1, drop = FALSE]
  } else {
    stop("format must be either \"long\" or \"abund\"")
  }
  return(out)
}

#' @title Check whether abundances appear to be read counts (i.e. not normalised)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#'
#' @return Logical
#' @keywords internal
abundAreCounts <- function(data) {
  ### Data must be in ampvis2 format
  is_ampvis2(data)

  # check if $abund contains read counts and not normalised counts or any decimals
  all(
    !isTRUE(attributes(data)$normalised),
    all(data$abund %% 1L == 0),
    !all(colSums(data$abund) == 100)
  )
}

#' @title Normalise read counts to 100, i.e. in percent relative abundance per sample
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#'
#' @return A modifed ampvis2 object
#' @keywords internal
normaliseTo100 <- function(data) {
  ### Data must be in ampvis2 format
  is_ampvis2(data)

  if (!abundAreCounts(data)) {
    warning("The data has already been normalised. Setting normalise = TRUE (the default) will normalise the data again and the relative abundance information about the original data of which the provided data is a subset will be lost.", call. = FALSE)
  }
  # normalise each sample to sample totals, skip samples with 0 sum to avoid NaN's
  tmp <- data$abund[, which(colSums(data$abund) != 0), drop = FALSE]
  if (nrow(tmp) == 1L) {
    # apply returns a vector and drops rownames if only 1 row, therefore set to 100 instead
    tmp[1L, ] <- 100L
  } else if (nrow(tmp) > 1L) {
    tmp <- as.data.frame(apply(tmp, 2, function(x) {
      x / sum(x) * 100
    }))
  }
  data$abund[, which(colSums(data$abund) != 0)] <- tmp
  attributes(data)$normalised <- TRUE
  return(data)
}

#' @title Filter species by a threshold in percent
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param filter_species Remove low abundant OTU's across all samples below this threshold in percent. (\emph{default}: \code{0})
#'
#' @importFrom ape drop.tip
#' @return An ampvis2 object
#' @keywords internal
filter_species <- function(data, filter_species = 0) {
  ### Data must be in ampvis2 format
  is_ampvis2(data)

  if (is.numeric(filter_species)) {
    if (filter_species > 0) {
      # For printing removed OTUs
      nOTUsbefore <- nrow(data$abund)

      # First transform to percentages
      abund_pct <- data$abund
      abund_pct[
        ,
        colSums(abund_pct) != 0
      ] <- as.data.frame(
        apply(
          abund_pct[, colSums(abund_pct) != 0, drop = FALSE],
          2,
          function(x) x / sum(x) * 100
        )
      )
      rownames(abund_pct) <- rownames(data$abund) # keep rownames

      # Then filter low abundant OTU's where ALL samples have
      # below the threshold set with filter_species in percent
      abund_subset <- abund_pct[
        !apply(
          abund_pct,
          1,
          function(row) all(row <= filter_species)
        ), ,
        drop = FALSE
      ] # remove low abundant OTU's

      data$abund <- data$abund[rownames(data$abund) %in% rownames(abund_subset), , drop = FALSE]
      rownames(data$tax) <- data$tax$OTU

      # also filter taxonomy, tree, and sequences
      data$tax <- data$tax[rownames(data$tax) %in% rownames(abund_subset), , drop = FALSE]

      if (!is.null(data$tree)) {
        data$tree <- ape::drop.tip(
          phy = data$tree,
          tip = data$tree$tip.label[!data$tree$tip.label %in% data$tax$OTU]
        )
      }
      if (!is.null(data$refseq)) {
        if (!is.null(names(data$refseq))) {
          # sometimes there is taxonomy alongside the OTU ID's. Anything after a ";" will be ignored
          names_stripped <- stringr::str_split(names(data$refseq), ";", simplify = TRUE)[, 1]
          data$refseq <- data$refseq[names_stripped %in% rownames(data$abund)]
        } else if (is.null(names(data$refseq))) {
          warning("DNA sequences have not been subsetted, could not find the names of the sequences in data$refseq.", call. = FALSE)
        }
      }
      nOTUsafter <- nrow(data$abund)
      if (nOTUsbefore == nOTUsafter) {
        message("0 OTU's have been filtered.")
      } else {
        message(
          paste0(
            nOTUsbefore - nOTUsafter,
            " OTUs not present in more than ",
            filter_species,
            "% relative abundance in any sample have been filtered \nBefore: ",
            nOTUsbefore,
            " OTUs\nAfter: ",
            nOTUsafter,
            " OTUs"
          )
        )
      }
    }
  }
  return(data)
}

#' @title Check if data has class "ampvis2"
#' @description Checks if the object is of class "ampvis2".
#'
#' @param data Object to check
#'
#' @return Returns error if \code{!inherits(data, "ampvis2")}, otherwise \code{invisible(TRUE)}.
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @keywords internal
is_ampvis2 <- function(data) {
  if (!inherits(data, "ampvis2")) {
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis2 functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)", call. = FALSE)
  }
  invisible(TRUE)
}
#' @title Valid taxonomic levels
#'
#' @description Valid taxonomic levels
#'
#' @name tax_levels
#' @docType data
#'
#' @format A charactor vector with the taxonomic levels Kingdom->OTU
#'
#' @keywords data internal
#'
tax_levels <- c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species",
  "OTU"
)

#' @title Check for installed package
#' @description Returns error if a required package is not installed. Mostly used for checking whether packages listed under the Suggests field in the DESCRIPTION file is installed.
#'
#' @param pkg The package to check for, a character of length 1.
#' @param msg Optionally additional text appended (with \code{paste0}) to the default error message.
#'
#' @return Returns error and message if not installed, otherwise \code{invisible(TRUE)}
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @keywords internal
checkReqPkg <- function(pkg, msg = "") {
  stopifnot(is.character(pkg), length(pkg) == 1L, nchar(pkg) > 0)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stopifnot(is.character(msg))
    stop(
      paste0("Package '", pkg, "' is required but not installed.", msg),
      call. = FALSE
    )
  }
  require(pkg, quietly = TRUE, character.only = TRUE)
}

#' @title Coerce DNAbin object to data.table
#' @description S3 generic function to create a data table from a DNAbin class object
#'
#' @param x (\emph{required}) A DNAbin class object
#' @param ... Not used
#' 
#' @importFrom data.table data.table
#'
#' @return a data.table
as.data.table.DNAbin <- function(x, ...) {
  dt <- data.table(
    name = names(x),
    seq = sapply(as.character(x), paste, collapse = "")
  )
  dt
}

#' @title Rename OTU's by sequence matching with a FASTA file
#' @description Match and rename OTU's in an ampvis2 object by sequence to a FASTA file
#'
#' @param data data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param fasta Path to a FASTA file or a \code{DNAbin} class object with sequences whose names will be used as OTU names by exact matches (i.e. same length, 100\% sequence identity). (\emph{default:} \code{NULL})
#' @param unmatched_prefix Prefix used to name any unmatched sequences when \code{refseq_names} is provided. An integer counting from 1 will be appended to this prefix, so for example the 123th unmatched sequence will be named \code{unmatched123}, and so on. (\emph{default:} \code{"unmatched"})
#' @param rename_unmatched Whether to rename any unmatched sequences or not when \code{refseq_names} is provided. (\emph{default:} \code{TRUE})
#'
#' @return An ampvis2 class object
matchOTUs <- function(
  data,
  fasta,
  unmatched_prefix = "unmatched",
  rename_unmatched = TRUE
) {
  #assure data is an ampvis2 class object
  is_ampvis2(data)
  
  if (!inherits(data$refseq, c("DNAbin", "AAbin"))) {
    stop("No DNA sequences available in the ampvis2 object", call. = FALSE)
  }
  
  #if fasta is a character expect a file path
  if (is.character(fasta) &
      length(fasta) == 1 &
      is.null(dim(fasta))) {
    fasta_dt <- as.data.table(read.FASTA(fasta))
  } else if (!inherits(fasta, c("DNAbin", "AAbin"))) {
    stop("refseq_names must be of class \"DNAbin\" or \"AAbin\" as loaded with the ape::read.FASTA() function.", call. = FALSE)
  } else if(inherits(fasta, c("DNAbin", "AAbin"))) {
    fasta_dt <- as.data.table(fasta)
  }
  
  refseq_dt <- as.data.table(data$refseq)
  
  #inner join (i.e. keep only rows in d_seqs_dt and order rows according to d_seqs_dt)
  merged_seqs <- fasta_dt[refseq_dt, on = "seq"]
  
  #generate new names for those with no match in ASV DB
  if (isTRUE(rename_unmatched)) {
    merged_seqs[is.na(name), name := paste0(unmatched_prefix, 1:nrow(.SD))]
  } else if (isFALSE(rename_unmatched)) {
    merged_seqs[is.na(name), name := i.name]
  }
  
  #rename everywhere in ampvis2 object
  rownames(data$abund) <- merged_seqs[["name"]]
  if(any(colnames(data$abund) %in% "OTU")) {
    data$abund[["OTU"]] <- merged_seqs[["name"]]
  }
  rownames(data$tax) <- merged_seqs[["name"]]
  data$tax[["OTU"]] <- merged_seqs[["name"]]
  names(data$refseq) <- merged_seqs[["name"]]
  
  return(data)
}

#' Replacement for ":::" to suppress R CMD CHECK warnings
# `%:::%` <- function(pkg, fun)
#  get(fun, envir = asNamespace(pkg), inherits = FALSE)
