#' @title Rarefy ampvis2 object (internal function)
#' @description This is just a wrapper of \code{\link[vegan]{rrarefy}} with convenient error messages and adjusted to work with ampvis2 objects.
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param rarefy (\emph{required}) Passed directly to \code{\link[vegan]{rrarefy}}.
#'
#' @return An ampvis2 object with rarefied OTU abundances.
#'
#' @importFrom vegan rrarefy
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
      data$abund <- suppressWarnings(vegan::rrarefy(t(data$abund), sample = rarefy)) %>% t() %>% as.data.frame()
      warning("The chosen rarefy size (", as.character(rarefy), ") is smaller than the smallest amount of reads in any sample (", as.character(minreads), ").", call. = FALSE)
    } else {
      data$abund <- suppressWarnings(vegan::rrarefy(t(data$abund), sample = rarefy)) %>% t() %>% as.data.frame()
      if (minreads < rarefy) {
        message("The following sample(s) have not been rarefied (less than ", as.character(rarefy), " reads):\n", paste(rownames(data$metadata[which(reads < rarefy), ]), collapse = ", "))
      }
    }
  } else if (!is.numeric(rarefy)) {
    stop("Argument rarefy must be numerical.", call. = FALSE)
  }
  return(data)
}

#' Tidy up taxonomy (internal function)
#'
#' Used internally in other ampvis functions.
#'
#' @usage amp_rename(data)
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
#'
#' @importFrom dplyr mutate
#'
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
amp_rename <- function(data, tax_class = NULL, tax_empty = "best", tax_level = "Genus") {
  tax <- data[["tax"]]

  ## First make sure that all entries are strings
  for (i in 1:ncol(tax)) {
    tax[, i] <- as.character(tax[, i])
  }

  ## Change a specific phylum to class level
  if (!is.null(tax_class)) {
    for (i in 1:nrow(tax)) {
      if (!is.na(tax$Phylum[i]) & tax$Phylum[i] %in% tax_class) {
        tax$Phylum[i] <- tax$Class[i]
      }
    }
  }

  ## Remove the underscore classifier from the data
  tax$Kingdom <- gsub("^k_*", "", tax$Kingdom)
  tax$Phylum <- gsub("^p_*", "", tax$Phylum)
  tax$Phylum <- gsub("^c_*", "", tax$Phylum)
  tax$Class <- gsub("^c_*", "", tax$Class)
  tax$Order <- gsub("^o_*", "", tax$Order)
  tax$Family <- gsub("^f_*", "", tax$Family)
  tax$Genus <- gsub("^g_*", "", tax$Genus)
  tax$Kingdom <- gsub("uncultured", "", tax$Kingdom)
  tax$Phylum <- gsub("uncultured", "", tax$Phylum)
  tax$Class <- gsub("uncultured", "", tax$Class)
  tax$Order <- gsub("uncultured", "", tax$Order)
  tax$Family <- gsub("uncultured", "", tax$Family)
  tax$Genus <- gsub("uncultured", "", tax$Genus)

  ## Check if there is a species level otherwise add it for consistency
  if (!is.null(tax$Species)) {
    tax$Species <- gsub("^s_*", "", tax$Species)
  } else {
    tax$Species <- ""
  }

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

#' Functional information tool
#'
#' Makes raw MiDAS function data compatible with ampvis format. Internal function, not exported.
#'
#' @usage amp_cleanMiF(data)
#'
#' @param data (required) A data frame with MiDAS functions.
#' @importFrom dplyr mutate
#' @return A data frame.
#'
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
amp_cleanMiF <- function(data) {
  MiF <- mutate(data,
    MiDAS = "POS",
    FIL = paste(Filamentous.Other, Filamentous.In.situ),
    FIL = ifelse(FIL %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", FIL),
    FIL = ifelse(FIL %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", FIL),
    FIL = ifelse(FIL %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", FIL),
    FIL = ifelse(FIL == "NT NT", "NT", FIL),

    AOB = paste(AOB.Other, AOB.In.situ),
    AOB = ifelse(AOB %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", AOB),
    AOB = ifelse(AOB %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", AOB),
    AOB = ifelse(AOB %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", AOB),
    AOB = ifelse(AOB == "NT NT", "NT", AOB),

    NOB = paste(NOB.Other, NOB.In.situ),
    NOB = ifelse(NOB %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", NOB),
    NOB = ifelse(NOB %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", NOB),
    NOB = ifelse(NOB %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", NOB),
    NOB = ifelse(NOB == "NT NT", "NT", NOB),

    Anammox = paste(Anammox.Other, Anammox.In.situ),
    Anammox = ifelse(Anammox %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", Anammox),
    Anammox = ifelse(Anammox %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", Anammox),
    Anammox = ifelse(Anammox %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", Anammox),
    Anammox = ifelse(Anammox == "NT NT", "NT", Anammox),

    AU.MIX = paste(Autotroph.Mixotroph.Other, Autotroph.Mixotroph.In.situ),
    AU.MIX = ifelse(AU.MIX %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", AU.MIX),
    AU.MIX = ifelse(AU.MIX %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", AU.MIX),
    AU.MIX = ifelse(AU.MIX %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", AU.MIX),
    AU.MIX = ifelse(AU.MIX == "NT NT", "NT", AU.MIX),

    PAO = paste(PAO.Other, PAO.In.situ),
    PAO = ifelse(PAO %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", PAO),
    PAO = ifelse(PAO %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", PAO),
    PAO = ifelse(PAO %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", PAO),
    PAO = ifelse(PAO == "NT NT", "NT", PAO),

    GAO = paste(GAO.Other, GAO.In.situ),
    GAO = ifelse(GAO %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", GAO),
    GAO = ifelse(GAO %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", GAO),
    GAO = ifelse(GAO %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", GAO),
    GAO = ifelse(GAO == "NT NT", "NT", GAO),

    HET = paste(Aerobic.heterotroph.Other, Aerobic.heterotroph.In.situ),
    HET = ifelse(HET %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", HET),
    HET = ifelse(HET %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", HET),
    HET = ifelse(HET %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", HET),
    HET = ifelse(HET == "NT NT", "NT", HET),

    DN = paste(Nitrite.reduction.Other, Nitrite.reduction.In.situ),
    DN = ifelse(DN %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", DN),
    DN = ifelse(DN %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", DN),
    DN = ifelse(DN %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", DN),
    DN = ifelse(DN == "NT NT", "NT", DN),

    FER = paste(Fermentation.Other, Fermentation.In.situ),
    FER = ifelse(FER %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", FER),
    FER = ifelse(FER %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", FER),
    FER = ifelse(FER %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", FER),
    FER = ifelse(FER == "NT NT", "NT", FER),

    SUL = paste(Sulphate.reduction.Other, Sulphate.reduction.In.situ),
    SUL = ifelse(SUL %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", SUL),
    SUL = ifelse(SUL %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", SUL),
    SUL = ifelse(SUL %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", SUL),
    SUL = ifelse(SUL == "NT NT", "NT", SUL),

    ACE = paste(Acetogen.Other, Acetogen.In.situ),
    ACE = ifelse(ACE %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", ACE),
    ACE = ifelse(ACE %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", ACE),
    ACE = ifelse(ACE %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", ACE),
    ACE = ifelse(ACE == "NT NT", "NT", ACE),

    MET = paste(Methanogen.Other, Methanogen.In.situ),
    MET = ifelse(MET %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", MET),
    MET = ifelse(MET %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", MET),
    MET = ifelse(MET %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", MET),
    MET = ifelse(MET == "NT NT", "NT", MET),

    FA = paste(Fatty.acids.Other, Fatty.acids.In.situ),
    FA = ifelse(FA %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", FA),
    FA = ifelse(FA %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", FA),
    FA = ifelse(FA %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", FA),
    FA = ifelse(FA == "NT NT", "NT", FA),

    SUG = paste(Sugars.Other, Sugars.In.situ),
    SUG = ifelse(SUG %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", SUG),
    SUG = ifelse(SUG %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", SUG),
    SUG = ifelse(SUG %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", SUG),
    SUG = ifelse(SUG == "NT NT", "NT", SUG),

    PRO = paste(Proteins.Amino.acids.Other, Proteins.Amino.acids.In.situ),
    PRO = ifelse(PRO %in% c("POS POS", "NEG POS", "NT POS", "POS NT", "POS VAR", "VAR POS"), "POS", PRO),
    PRO = ifelse(PRO %in% c("NEG NEG", "POS NEG", "NT NEG", "NEG NT", "VAR NEG"), "NEG", PRO),
    PRO = ifelse(PRO %in% c("VAR NT", "NT VAR", "NEG VAR", "VAR VAR"), "VAR", PRO),
    PRO = ifelse(PRO == "NT NT", "NT", PRO)
  )
  return(MiF)
}

#' Calculate weighted or unweighted UniFrac distances. Adopted from fastUniFrac() from phyloseq
#'
#' @param abund Abundance table with OTU counts, in \code{ampvis2} objects it is available with simply data$abund
#' @param tree Phylogenetic tree (rooted and with branch lengths) as loaded with \code{\link[ape]{read.tree}}.
#' @param weighted Calculate weighted or unweighted UniFrac distances.
#' @param normalise Should the output be normalised such that values range from 0 to 1 independent of branch length values? Note that unweighted UniFrac is always normalised. (\emph{default:} \code{TRUE})
#' @param num_threads The number of threads to be used for calculating UniFrac distances. If set to more than \code{1} then this is set by using \code{\link[doParallel]{registerDoParallel}} (\emph{default:} \code{1})
#' @importFrom ape is.rooted node.depth node.depth.edgelength reorder.phylo
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach registerDoSEQ %dopar%
#' @return A distance matrix of class \code{dist}
unifrac <- function(abund,
                    tree,
                    weighted = FALSE,
                    normalise = TRUE,
                    num_threads = 1L) {
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
    stop("OTU table and phylogenetic tree do not match. (Note: This may be the result of subsetting if the provided data is a subset of larger data as phylogenetic trees are not subsetted)", call. = FALSE)
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
      numerator <- sum({
        tree$edge.length * wUF_branchweight
      }, na.rm = TRUE)
      if (isFALSE(normalise)) {
        return(numerator)
      } else {
        denominator <- sum({
          tipAges * (OTU[, A] / AT + OTU[, B] / BT)
        }, na.rm = TRUE)
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

#' Find lowest taxonomic level
#'
#' @description Finds the lowest taxonomic level of two given levels in \code{tax}. The hierarchic order of taxonomic levels is simply taken from the order of column names in \code{tax}
#'
#' @param tax_aggregate A character vector of one or more taxonomic levels (exactly as in the column names of \code{ampvis2obj$tax}), fx Genus or Species. (\emph{default:} \code{NULL})
#' @param tax_add A second character vector similar to \code{tax_aggregate}. (\emph{default:} \code{NULL})
#' @param tax The taxonomy table from an ampvis2 object (\code{ampvis2obj$tax})
#'
#' @return A length one character vector with the lowest taxonomic level
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

#' Aggregate OTUs to a specific taxonomic level
#'
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
#' @return A data.table
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
  newabund <- abund[rownames(newtax), ]
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
    value.name = "abundance",
    variable.factor = FALSE
  )
  if (nrow(newtax) != nrow(tax) & isTRUE(calcSums)) {
    abundAggr[,
      Sum := sum(abundance),
      keyby = .(Display, Sample)
    ]
  }
  if (format == "long") {
    out <- abundAggr
  } else if (format == "abund") {
    out <- as.data.frame(
      data.table::dcast(abundAggr,
        Display ~ Sample,
        value.var = "abundance",
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
