#' Load data for ampvis2 functions
#'
#' This function reads an OTU-table and corresponding sample metadata, and returns a list for use in all ampvis2 functions. It is therefore required to load data with \code{\link{amp_load}} before any other ampvis2 functions can be used.
#'
#' @param otutable (\emph{required}) OTU-table with the read counts of all OTU's. Rows are OTU's, columns are samples, otherwise you must transpose. The taxonomy of the OTU's can be placed anywhere in the table and will be extracted by name (Kingdom/Domain -> Species). Can be a data frame, matrix, or path to a delimited text file or excel file which will be read using either \code{\link[data.table]{fread}} or \code{\link[readxl]{read_excel}}, respectively.
#' @param metadata (\emph{recommended}) Sample metadata with any information about the samples. The first column must contain sample ID's matching those in the otutable. If none provided, dummy metadata will be created. Can be a data frame, matrix, or path to a delimited text file or excel file which will be read using either \code{\link[data.table]{fread}} or \code{\link[readxl]{read_excel}}, respectively. (\emph{default:} \code{NULL})
#' @param taxonomy (\emph{recommended}) Taxonomy table where rows are OTU's and columns are up to 7 levels of taxonomy named Kingdom/Domain->Species. If taxonomy is also present in otutable, it will be discarded and only this will be used. Can be a data frame, matrix, or path to a delimited text file or excel file which will be read using either \code{\link[data.table]{fread}} or \code{\link[readxl]{read_excel}}, respectively. (\emph{default:} \code{NULL})
#' @param fasta (\emph{optional}) Path to a FASTA file with reference sequences for all OTU's in the OTU-table. (\emph{default:} \code{NULL})
#' @param tree (\emph{optional}) Path to a phylogenetic tree file which will be read using \code{\link[ape]{read.tree}}, or an object of class \code{"phylo"}. (\emph{default:} \code{NULL})
#' @param pruneSingletons (\emph{logical}) Remove OTU's only observed once in all samples. (\emph{default:} \code{FALSE})
#' @param ... (\emph{optional}) Additional arguments are passed on to any of the file reader functions used.
#'
#' @return A list of class \code{"ampvis2"} with 3 to 5 elements.
#'
#' @importFrom ape read.FASTA
#' @importFrom stringr str_replace_all str_to_title
#' @importFrom dplyr intersect mutate_at
#' @importFrom data.table fread
#' @importFrom readxl read_excel
#' @importFrom tools file_ext
#'
#' @export
#'
#' @details The \code{\link{amp_load}} function validates and corrects the provided data frames in different ways to make it suitable for the rest of the ampvis2 functions. It is important that the provided data frames match the requirements as described in the following sections to work properly.
#'
#' @section The OTU-table:
#' The OTU-table contains information about the OTUs, their read counts in each sample, and optionally their assigned taxonomy. The provided OTU-table must be a data frame with the following requirements:
#'
#' \itemize{
#'   \item The rows are OTU IDs and the columns are samples.
#'   \item The last 7 columns are optionally the corresponding taxonomy assigned to the OTUs, named \code{"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"}.
#'   \item The OTU ID's are expected to be in either the row names of the data frame or in a column called "OTU", "ASV", or "#OTU ID". Otherwise the function will stop with a message.
#'   \item The column names of the data frame are the sample IDs, exactly matching those in the metadata, (and taxonomy columns named \code{Kingdom} -> \code{Species} if present, of course).
#'   \item Generally avoid special characters and spaces in row- and column names.
#' }
#'
#' A minimal example is available with \code{data("example_otutable")}.
#'
#' @section The metadata:
#' The metadata contains additional information about the samples, for example where each sample was taken, date, pH, treatment etc, which is used to compare and group the samples during analysis. The amount of information in the metadata is unlimited, it can contain any number of columns (variables), however there are a few requirements:
#'
#' \itemize{
#'   \item The sample IDs must be in the first column. These sample IDs must match exactly to those in the OTU-table.
#'   \item Column classes matter, categorical variables should be loaded either \code{as.character()} or \code{as.factor()}, and continuous variables \code{as.numeric()}. See below.
#'   \item Generally avoid special characters and spaces in row- and column names.
#' }
#'
#' If for example a column is named "Year" and the entries are simply entered as numbers (2011, 2012, 2013 etc), then R will automatically consider these as numerical values (\code{as.numeric()}) and therefore the column as a continuous variable, while it is a categorical variable and should be loaded \code{as.factor()} or \code{as.character()} instead. This has consequences for the analysis as R treats them differently. Therefore either use the \code{colClasses = } argument when loading a csv file or \code{col_types = } when loading an excel file, or manually adjust the column classes afterwards with fx \code{metadata$Year <- as.character(metadata$Year)}.
#'
#' The \code{\link{amp_load}} function will automatically use the sample IDs in the first column as row names, but it is important to also have an actual column with sample IDs, so it's possible to fx group by that column during analysis. Any unmatched samples between the otutable and metadata will be removed with a warning.
#'
#' A minimal example is available with \code{data("example_metadata")}.
#'
#' @seealso
#' \code{\link{amp_load}}, \code{\link{amp_subset_samples}}, \code{\link{amp_subset_taxa}}
#'
#' @examples
#'
#' library(ampvis2)
#' \dontrun{
#' # Load data by either giving file paths or by passing already loaded R objects
#' ### example load with file paths
#' d <- amp_load(
#'   otutable = "path/to/otutable.tsv",
#'   metadata = "path/to/metadata.xlsx",
#'   taxonomy = "path/to/taxonomy.txt"
#' )
#'
#' ### example load with R objects
#' # Read the OTU-table as a data frame. It is important to set check.names = FALSE
#' myotutable <- read.delim("data/otutable.txt", check.names = FALSE)
#'
#' # Read the metadata, probably an excel sheet
#' mymetadata <- read_excel("data/metadata.xlsx", col_names = TRUE)
#'
#' # Read the taxonomy
#' mytaxonomy <- read.csv("data/taxonomy.csv", check.names = FALSE)
#'
#' # Combine the data with amp_load()
#' d <- amp_load(
#'   otutable = myotutable,
#'   metadata = mymetadata,
#'   taxonomy = mytaxonomy,
#'   pruneSingletons = FALSE,
#'   fasta = "path/to/fastafile.fa", # optional
#'   tree = "path/to/tree.tree" # optional
#' )
#'
#' ### Show a short summary about the data by simply typing the name of the object in the console
#' d
#' }
#'
#' ### Minimal example metadata:
#' data("example_metadata")
#' example_metadata
#'
#' ### Minimal example otutable:
#' data("example_otutable")
#' example_otutable
#'
#' ### Minimal example taxonomy:
#' data("example_taxonomy")
#' example_taxonomy
#'
#' # load example data
#' d <- amp_load(
#'   otutable = example_otutable,
#'   metadata = example_metadata,
#'   taxonomy = example_taxonomy
#' )
#' 
#' # show a summary of the data
#' d
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#'
amp_load <- function(otutable,
                     metadata = NULL,
                     taxonomy = NULL,
                     fasta = NULL,
                     tree = NULL,
                     pruneSingletons = FALSE,
                     ...) {
  ### the following functions are only useful in the context of amp_load()
  # function to check if provided object looks like a file path
  # and if so try to read the file, otherwise expect a data.frame
  import <- function(x, ...) {
    if (is.character(x) &
      length(x) == 1 &
      is.null(dim(x))) {
      ext <- tolower(tools::file_ext(x))
      if (ext %in% c("csv", "txt", "tsv", "")) {
        DF <- data.table::fread(x, fill = TRUE, data.table = FALSE, ...)
      } else if (ext %in% c("xls", "xlsx")) {
        DF <- readxl::read_excel(x, ...)
      } else {
        stop(paste("Unsupported file type \"", ext, "\""), call. = FALSE)
      }
    } else {
      DF <- x
    }
    DF <- as.data.frame(DF, check.names = FALSE)
    if (sum(dim(DF)) == 0L) {
      stop(paste(deparse(substitute(x)), "is empty"), call. = FALSE)
    }
    return(DF)
  }

  # function to detect OTU/ASV column
  findOTUcol <- function(x) {
    DF <- as.data.frame(x)
    otucol <- tolower(colnames(DF)) %in% c("otu", "asv", "#otu id")
    if (sum(otucol) > 1L) {
      stop(
        paste(
          "More than one column in",
          deparse(substitute(x)),
          "is named OTU/ASV, don't know which one to use."
        ),
        call. = FALSE
      )
    }
    if (any(otucol)) {
      rownames(DF) <- as.character(DF[, which(otucol)])
      colnames(DF)[which(otucol)] <- "OTU"
    } else if (!any(otucol)) {
      warning(
        paste0(
          "Could not find a column named OTU/ASV in ",
          deparse(substitute(x)),
          ", using rownames as sample ID's"
        ),
        call. = FALSE
      )
      DF[["OTU"]] <- rownames(DF)
    }
    return(DF)
  }

  # function to extract and parse taxonomy correctly
  tax.levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")
  parseTaxonomy <- function(x) {
    # create empty dummy taxonomy data frame to fill into
    tax <- as.data.frame(
      matrix(
        nrow = nrow(x),
        ncol = length(tax.levels),
        dimnames = list(rownames(x), tax.levels)
      )
    )
    tax[["OTU"]] <- rownames(tax)

    # rename all taxonomy columns except OTU column
    taxcols <- tolower(colnames(x)) %in% c("domain", tolower(tax.levels[-8]))
    colnames(x)[taxcols] <- stringr::str_to_title(colnames(x)[taxcols])

    # identify which columns contain taxonomy
    taxcolnames <- dplyr::intersect(tax.levels, colnames(x))

    # allow both Kingdom or Domain column, but not both at once
    if (all(c("Domain", "Kingdom") %in% taxcolnames)) {
      stop("Cannot have both Domain and Kingdom columns at the same time in taxonomy.", call. = FALSE)
    } else if (any("Domain" %in% taxcolnames)) {
      tax.levels[1] <- "Domain" -> colnames(tax)[1]
    }

    # fill into dummy taxonomy by merging by OTU
    tax <- merge(
      tax[, !colnames(tax) %in% taxcolnames[!taxcolnames %in% "OTU"], drop = FALSE],
      x[, taxcolnames, drop = FALSE],
      by = "OTU",
      sort = FALSE,
      all.x = TRUE,
      all.y = FALSE
    )
    rownames(tax) <- tax[["OTU"]]

    if (length(taxcolnames) == 1L & all(taxcolnames %in% "OTU")) {
      warning("Could not find or parse taxonomy, creating a dummy taxonomy table with only OTUs.", call. = FALSE)
    }

    # select and sort columns correctly in Kingdom/Domain -> Species order
    tax <- tax[, tax.levels, drop = FALSE]

    # remove whitespaces at the either side
    tax[] <- lapply(tax, stringr::str_replace_all, pattern = "^\\s+|\\s+$", replacement = "")

    # only character columns allowed and convert any NA to empty string
    tax[] <- lapply(tax, as.character)
    tax[is.na(tax)] <- ""

    return(tax)
  }

  ### import and check otutable (with or without taxonomy)
  otutable <- import(otutable)
  otutable <- findOTUcol(otutable)

  ### extract read abundances from otutable (same if taxonomy present or not)
  taxcols <- tolower(colnames(otutable)) %in% c("domain", tolower(tax.levels))
  abund <- otutable[, !taxcols, drop = FALSE]
  abund[is.na(abund)] <- 0L

  ### extract taxonomy from otutable or separate table if provided
  if (is.null(taxonomy)) {
    tax <- parseTaxonomy(otutable)
  } else if (!is.null(taxonomy)) {
    taxonomy <- import(taxonomy)
    taxonomy <- findOTUcol(taxonomy)
    tax <- parseTaxonomy(taxonomy)

    # check if provided otutable and taxonomy match, message with missing/excess OTUs
    if (!setequal(rownames(tax), rownames(abund))) {
      sharedOTUs <- dplyr::intersect(rownames(tax), rownames(abund))
      if (length(sharedOTUs) == 0L) {
        stop("No OTU's match between otutable and taxonomy", call. = FALSE)
      } else {
        warningMsg <- "The OTU's between otutable and taxonomy do not match exactly."
        if (all(sharedOTUs %in% rownames(abund)) && nrow(tax) > length(sharedOTUs)) {
          nOTUremoved <- nrow(tax) - length(sharedOTUs)
          warningMsg <- paste0(
            warningMsg,
            " ",
            nOTUremoved,
            " OTU",
            if (nOTUremoved > 1L) {
              "'s in taxonomy not present in otutable have "
            } else if (nOTUremoved == 1L) {
              " in taxonomy not present in otutable has "
            },
            "been removed from taxonomy."
          )
        }
        if (!all(rownames(abund) %in% sharedOTUs)) {
          OTUmissing <- setdiff(rownames(abund), sharedOTUs)
          warningMsg <- paste0(
            warningMsg,
            " ",
            length(OTUmissing),
            " OTU",
            if (length(OTUmissing) == 1L) {
              paste0(" (", OTUmissing, ") is ")
            } else {
              "'s are "
            },
            "missing in taxonomy compared to otutable",
            if (length(OTUmissing) == 1L) {
              "."
            } else {
              paste0(", some of which are:\n\"", paste0(head(OTUmissing), collapse = "\", \""), "\"")
            }
          )
        }
      }
      warning(warningMsg, call. = FALSE)
    }
  }

  # filter and order tax based on abund
  tax <- tax[rownames(tax) %in% rownames(abund), , drop = FALSE]

  ### prune singletons, check if abundances are whole numbers
  if (isTRUE(pruneSingletons)) {
    if (any(abund %% 1 != 0L)) {
      stop("The data contains non-integer values, could not identify singletons.", call. = FALSE)
    } else {
      singletons <- as.integer(rowSums(abund)) == 1L
      if (!any(singletons)) {
        message("No singletons found in the data, none removed")
      } else if (any(singletons)) {
        nOTUsBefore <- nrow(abund)
        abund <- abund[-which(singletons), , drop = FALSE]
        tax <- tax[rownames(tax) %in% rownames(abund), , drop = FALSE]
        nSingletons <- nOTUsBefore - nrow(abund)
        message(
          paste0(
            nSingletons,
            if (nSingletons == 1L) {
              " singleton has "
            } else {
              " singletons have "
            },
            "been removed"
          )
        )
      }
    }
  }

  ### create dummy metadata if none provided
  if (!is.null(metadata)) {
    metadata <- import(metadata)
  } else if (is.null(metadata)) {
    metadata <- data.frame(
      "SampleID" = colnames(abund),
      "DummyVariable" = "All samples",
      check.names = FALSE
    )
    warning("No sample metadata provided, creating dummy metadata.\n", call. = FALSE)
  }

  # sample ID's must be character
  metadata[[1]] <- as.character(metadata[[1]])
  rownames(metadata) <- metadata[[1]]

  # only alphanumeric characters in metadata column names, otherwise replace with "_"
  colnames(metadata) <- stringr::str_replace_all(colnames(metadata), "[^[:alnum:]]", "_")

  abund0 <- abund
  metadata0 <- metadata

  ### check if metadata and otutable match
  if (!setequal(rownames(metadata), colnames(abund))) {
    if (!any(colnames(abund) %in% rownames(metadata))) {
      # No samples match
      stop("No sample names match between metadata and otutable. Check that you have loaded matching files and that they meet the requirements described in ?amp_load(). Remember to use check.names = FALSE when loading the files.", call. = FALSE)
    } else {
      # Find matching samples between otutable and metadata
      sharedSamples <- dplyr::intersect(colnames(abund), rownames(metadata))

      # subset both based on shared samples
      abund0 <- abund[, match(sharedSamples, colnames(abund)), drop = FALSE]
      abund0 <- abund0[rowSums(abund0) > 0, , drop = FALSE] # after subset, rows with all 0's may be introduced, remove
      metadata0 <- metadata[match(sharedSamples, rownames(metadata)), , drop = FALSE]

      # Vectors with unique sample names in either
      metadataUniques <- metadata[-which(rownames(metadata) %in% rownames(metadata0)), , drop = FALSE] %>% rownames()
      abundUniques <- abund[, -which(colnames(abund) %in% colnames(abund0)), drop = FALSE] %>% colnames()

      # print the removed samples
      warning(
        "Only ",
        ncol(abund0),
        " of ",
        length(unique(c(rownames(metadata), colnames(abund)))),
        " unique sample names match between metadata and otutable. The following unmatched samples have been removed:",
        ifelse(
          length(metadataUniques) > 0,
          paste0("\nmetadata (", length(metadataUniques), "): \n\t\"", paste(metadataUniques, collapse = "\", \""), "\""),
          ""
        ),
        ifelse(
          length(abundUniques) > 0,
          paste0("\notutable (", length(abundUniques), "): \n\t\"", paste(abundUniques, collapse = "\", \""), "\""),
          ""
        ),
        call. = FALSE
      )
    }
  }

  ### data: return the data in a combined list w or w/o refseq. The rows of tax are ordered by
  # the rownames of abund, and the columns of abund are ordered by the metadata rownames
  data <- list(
    abund = abund0[, rownames(metadata0), drop = FALSE],
    tax = tax[rownames(tax) %in% rownames(abund0), , drop = FALSE],
    metadata = metadata0
  )

  ### append refseq if provided
  # check tree
  if (!is.null(fasta)) {
    if (is.character(fasta) &
      length(fasta) == 1 &
      is.null(dim(fasta))) {
      refseq <- ape::read.FASTA(fasta, ...)[rownames(abund0)]
    } else if (!inherits(fasta, c("DNAbin", "AAbin"))) {
      stop("fasta must be of class \"DNAbin\" or \"AAbin\" as loaded with the ape::read.FASTA() function.", call. = FALSE)
    }
    data$refseq <- refseq
  }

  ### append phylogenetic tree if provided
  # check tree
  if (!is.null(tree)) {
    if (is.character(tree) &
      length(tree) == 1 &
      is.null(dim(tree))) {
      tree <- ape::read.tree(tree, ...)
    } else if (!inherits(tree, "phylo")) {
      stop("The provided phylogenetic tree must be of class \"phylo\" as loaded with the ape::read.tree() function.", call. = FALSE)
    }
    tree <- ape::drop.tip(
      phy = tree,
      tip = tree$tip.label[!tree$tip.label %in% data$tax$OTU]
    )
    data$tree <- tree
  }

  return(
    structure(
      data,
      class = "ampvis2",
      normalised = FALSE
    )
  )
}
