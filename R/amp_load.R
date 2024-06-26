### the following functions are only useful in the context of amp_load()
#' @title Extract/parse taxonomy from a data.frame
#' 
#' @param x Data frame or file path
#'
#' @return A data frame
#' @keywords internal
parseTaxonomy <- function(x) {
  #assert object is a data frame here
  
  # create empty dummy taxonomy data frame to fill into
  tax <- as.data.frame(
    matrix(
      nrow = nrow(x),
      ncol = length(tax_levels),
      dimnames = list(rownames(x), tax_levels)
    )
  )
  tax[["OTU"]] <- rownames(tax)
  
  # rename all taxonomy columns except OTU column
  taxcols <- tolower(colnames(x)) %in% c("domain", tolower(tax_levels[-8]))
  colnames(x)[taxcols] <- stringr::str_to_title(colnames(x)[taxcols])
  
  # identify which columns contain taxonomy
  taxcolnames <- dplyr::intersect(tax_levels, colnames(x))
  
  # allow both Kingdom or Domain column, but not both at once
  if (all(c("Domain", "Kingdom") %in% taxcolnames)) {
    stop("Cannot have both Domain and Kingdom columns at the same time in taxonomy.", call. = FALSE)
  } else if (any("Domain" %in% taxcolnames)) {
    tax_levels[1] <- "Domain" -> colnames(tax)[1]
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
  tax <- tax[, tax_levels, drop = FALSE]
  
  # remove whitespaces at the either side
  tax[] <- lapply(tax, stringr::str_replace_all, pattern = "^\\s+|\\s+$", replacement = "")
  
  # only character columns allowed and convert any NA to empty string
  tax[] <- lapply(tax, as.character)
  tax[is.na(tax)] <- ""
  
  return(tax)
}

#' @title Detect OTU/ASV column in a data frame
#'
#' @param x A data frame
#' @param OTUcolname Character vector with columns names to search for OTUs/ASVs
#'
#' @return A data frame
#' @keywords internal
findOTUcol <- function(
  x,
  OTUcolname = c("OTU", "ASV", "#OTU ID")
) {
  DF <- as.data.frame(x)
  otucol <- tolower(colnames(DF)) %in% tolower(OTUcolname)
  if (sum(otucol) > 1L) {
    stop(
      paste(
        "More than one column in",
        deparse(substitute(x)),
        "contains OTUs/ASVs. I don't know which one to use."
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
        "Could not find a column named ",
        if(length(OTUcolname) > 1L) "one of \"",
        paste(OTUcolname, collapse = ", "),
        "\" in ",
        deparse(substitute(x)),
        ". Using row names as OTU ID's",
        if(identical(rownames(DF), as.character(1:nrow(DF))))
          ", however they seem to be just numbers and are likely made by R"
        else
          "."
      ),
      call. = FALSE
    )
    DF[["OTU"]] <- rownames(DF)
  }
  return(DF)
}

#' @title Import any file format
#' @description
#' Excel, txt, CSV, TSV, BIOM, sintax. If already a data frame, will return as-is.
#'
#' @param x File path or data frame.
#' @param ... Additional arguments passed on to read_excel() and fread().
#'
#' @return A data frame
#' @keywords internal
import <- function(x, ...) {
  #enlist additional arguments to be able to parse correctly
  add_args <- list(...)
  
  # check if x is a length 1 string (i.e. a file path)
  if (is.character(x) &
      length(x) == 1 &
      is.null(dim(x))) {
    ext <- tolower(tools::file_ext(x))
    
    # use fread() if it's a delimited text file
    if (ext %in% c("csv", "txt", "tsv", "gz", "zip", "bz2", "sintax", "")) {
      fread_args <- add_args[names(add_args) %chin% names(formals(fread))]
      
      # Currently, sintax format taxonomy CANNOT be loaded if gzip or bz2, only zip
      x <- unzip_file(x)
      ext <- tools::file_ext(x)
      
      # if ext is .sintax expect sintax format and parse accordingly
      if (ext %in% "sintax") {
        DF <- do.call(
          fread,
          c(
            input = x,
            sep = "\t",
            fill = TRUE,
            header = FALSE,
            data.table = TRUE,
            fread_args
          )
        )[, c(1, 4)]
        colnames(DF) <- c("OTU", "tax")
        
        # Separate each taxonomic level into individual columns.
        # This has to be done in separate lines as taxonomic 
        # levels can be blank in between two other levels.
        DF[, "Kingdom" := gsub("[dk]:", "k__", stringr::str_extract(tax, "[dk]:[^,]*"))]
        DF[, "Phylum" := gsub("p:", "p__", stringr::str_extract(tax, "p:[^,]*"))]
        DF[, "Class" := gsub("c:", "c__", stringr::str_extract(tax, "c:[^,]*"))]
        DF[, "Order" := gsub("o:", "o__", stringr::str_extract(tax, "o:[^,]*"))]
        DF[, "Family" := gsub("f:", "f__", stringr::str_extract(tax, "f:[^,]*"))]
        DF[, "Genus" := gsub("g:", "g__", stringr::str_extract(tax, "g:[^,]*"))]
        DF[, "Species" := gsub("s:", "s__", stringr::str_extract(tax, "s:[^,]*"))]
        DF <- DF[, -2]
        
        # coerce to data.frame
        data.table::setDF(DF, rownames = DF[, OTU])
      } else {
        DF <- do.call(
          fread,
          c(
            input = x,
            data.table = FALSE,
            fill = TRUE,
            fread_args
          ))
      }
      # use read_excel() if xls, xlsx
    } else if (ext %in% c("xls", "xlsx")) {
      checkReqPkg("readxl")
      DF <- readxl::read_excel(x, ...)
      # if ext is .biom expect BIOM format and parse correctly
    } else if (ext %in% "biom") {
      checkReqPkg(
        "biomformat",
        " Please install with:\n  install.packages(\"BiocManager\"); BiocManager::install(\"biomformat\")"
      )
      
      biom <- biomformat::read_biom(x)
      
      # extract read counts
      abund <- biomformat::biom_data(biom) %>%
        as.matrix(check.names = FALSE) %>%
        as.data.frame(check.names = FALSE)
      
      # if only one sample biom_data() drops the name of it for some reason
      if (ncol(abund) == 1L) {
        colnames(abund) <- biom$columns[[1]][["id"]]
      }
      
      # extract the taxonomy
      taxlist <- lapply(biom$rows, function(x) {
        x$metadata$taxonomy
      })
      
      # check if taxonomy is empty for all OTU's, use ID's as OTU's if so
      if (all(is.null(unlist(taxlist, use.names = FALSE)))) {
        warning("Could not find the taxonomy of one or more OTU's in the provided .biom file", call. = FALSE)
        DF <- abund
      } else {
        # extract OTU names
        names(taxlist) <- lapply(biom$rows, function(x) {
          x$id
        })
        
        # coerce to data frame
        tax <- as.data.frame(t(as.data.frame(taxlist, check.names = FALSE, stringsAsFactors = FALSE)))
        
        # rename taxonomic levels
        colnames(tax) <- tax_levels[1:ncol(tax)]
        
        if (ncol(tax) < 7L) {
          warning("Taxonomy had less than 7 levels for all OTU's (Kingdom->Species), filling with NA from Species level and up.", call. = FALSE)
        }
        
        # combine abundances and taxonomy and return
        DF <- cbind(abund, tax) # no need for merge()
      }
      
      # extract the sample metadata
      metadatalist <- lapply(biom$columns, function(x) {
        x$metadata
      })
      
      # check whether metadata is empty
      if (all(is.null(unlist(metadatalist, use.names = FALSE)))) {
        warning("Could not find any sample metadata in the provided .biom file", call. = FALSE)
        metadata <- NULL
      } else {
        # extract sample names
        names(metadatalist) <- lapply(biom$columns, function(x) {
          x$id
        })
        
        # coerce to data table, then data frame
        metadata <- setDF(rbindlist(metadatalist, idcol = "SampleID"))
        
        # append metadata to DF as an attribute
        # (best solution for now)
        attr(DF, "metadata") <- metadata
      }
    } else {
      stop(paste0("Unsupported file type \".", ext, "\""), call. = FALSE)
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

#' Load data for ampvis2 functions
#'
#' This function reads an OTU-table and corresponding sample metadata, and returns a list for use in all ampvis2 functions. It is therefore required to load data with \code{\link{amp_load}} before any other ampvis2 functions can be used.
#'
#' @param otutable (\emph{required}) File path, data frame, or a phyloseq-class object. OTU-table with the read counts of all OTU's. Rows are OTU's, columns are samples, otherwise you must transpose. The taxonomy of the OTU's can be placed anywhere in the table and will be extracted by name (Kingdom/Domain -> Species). If a file path is provided it will be attempted being read by either \code{\link[data.table]{fread}} or \code{\link[readxl]{read_excel}}, respectively. Compressed files (zip, bzip2, gzip) are supported if not an excel file (bzip2 and gzip requires \code{data.table} 1.14.3 or later). Can also be a path to a BIOM file, which will then be parsed using the \href{https://github.com/joey711/biomformat}{biomformat} package, so both the JSON and HDF5 versions of the BIOM format are supported.
#' @param metadata (\emph{recommended}) File path or a data frame. Sample metadata with any information about the samples. The first column must contain sample ID's matching those in the otutable. If none provided, dummy metadata will be created. Can be a data frame, matrix, or path to a delimited text file or excel file which will be read using either \code{\link[data.table]{fread}} or \code{\link[readxl]{read_excel}}, respectively. Compressed files (zip, bzip2, gzip) are supported if not an excel file (bzip2 and gzip requires \code{data.table} 1.14.3 or later). If \code{otutable} is a BIOM file and contains sample metadata, \code{metadata} will take precedence if provided. (\emph{default:} \code{NULL})
#' @param taxonomy (\emph{recommended}) File path or a data frame. Taxonomy table where rows are OTU's and columns are up to 7 levels of taxonomy named Kingdom/Domain->Species. If taxonomy is also present in otutable, it will be discarded and only this will be used. Can be a data frame, matrix, or path to a delimited text file or excel file which will be read using either \code{\link[data.table]{fread}} or \code{\link[readxl]{read_excel}}, respectively. Compressed files (zip, bzip2, gzip) are supported if not an excel file (bzip2 and gzip requires \code{data.table} 1.14.3 or later). Can also be a path to a .sintax taxonomy table from a \href{http://www.drive5.com/usearch/}{USEARCH} analysis \href{http://www.drive5.com/usearch/manual/ex_miseq.html}{pipeline}, file extension must be \code{.sintax}. bzip2 or gzip compression is currently NOT supported if sintax format. (\emph{default:} \code{NULL})
#' @param fasta (\emph{optional}) Path to a FASTA file with reference sequences for all OTU's in the OTU-table. (\emph{default:} \code{NULL})
#' @param tree (\emph{optional}) Path to a phylogenetic tree file which will be read using \code{\link[ape]{read.tree}}, or an object of class \code{"phylo"}. (\emph{default:} \code{NULL})
#' @param pruneSingletons (\emph{logical}) Remove OTU's only observed once in all samples. (\emph{default:} \code{FALSE})
#' @param removeAbsentOTUs (\emph{logical}) Remove OTU's with 0 abundance in all samples. Absent OTUs are rarely present in the input data itself, but can occur when some samples are removed because of a mismatch between samples in the OTU-table and sample metadata. (\emph{default:} \code{TRUE})
#' @param otutable_OTUcolname Character vector with the name(s) of the column in the otutable that contains the OTUs/ASVs.
#' @param taxonomy_OTUcolname Character vector with the name(s) of the column in the taxonomy that contains the OTUs/ASVs.
#' @param ... (\emph{optional}) Additional arguments are passed on to any of the file reader functions used.
#'
#' @return A list of class \code{"ampvis2"} with 3 to 5 elements.
#'
#' @importFrom ape read.FASTA
#' @importFrom stringr str_replace_all str_to_title
#' @importFrom dplyr intersect mutate_at
#' @importFrom data.table fread setDF rbindlist
#' @importFrom tools file_ext
#' @importFrom utils head
#'
#' @export
#'
#' @details The \code{\link{amp_load}} function validates and corrects the provided data frames in different ways to make it suitable for the rest of the ampvis2 functions. It is important that the provided data frames match the requirements as described in the following sections to work properly. If a \code{phyloseq}-class object is provided the metadata, taxonomy, fasta, and tree arguments are ignored as they are expected to be provided in the \code{phyloseq} object.
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
#' \code{\link{amp_load}}, \code{\link{amp_filter_samples}}, \code{\link{amp_filter_taxa}}
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
#' ### Load a phyloseq object
#' d <- amp_load(physeq_object)
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
                     removeAbsentOTUs = TRUE,
                     otutable_OTUcolname = c("OTU", "ASV", "#OTU ID"),
                     taxonomy_OTUcolname = c("OTU", "ASV", "#OTU ID"),
                     ...) {
  ### Check whether otutable is a phyloseq object
  if(inherits(otutable, "phyloseq")) {
    message("Phyloseq object provided. Ignoring anything provided for the metadata, taxonomy, fasta, or tree arguments, using those from the phyloseq object instead.")
    physeq <- otutable
    checkReqPkg(
      "phyloseq",
      " Please install with:\n  install.packages(\"BiocManager\"); BiocManager::install(\"phyloseq\")"
    )
    #ampvis2 requires taxonomy and abundance table, phyloseq checks for the latter
    if(is.null(otutable@tax_table))
      stop("No taxonomy found in the phyloseq object and is required for ampvis2", call. = FALSE)
    
    #OTUs must be in rows, not columns
    if(taxa_are_rows(physeq)) {
      otutable <- as.data.frame(otu_table(physeq)@.Data)
    } else {
      otutable <- as.data.frame(t(otu_table(physeq)@.Data))
    }
    
    #tax_table is assumed to have OTUs in rows too
    taxonomy <- tax_table(physeq)@.Data
    
    #extract sample_data (metadata)
    if(!is.null(physeq@sam_data)) {
      metadata <- data.frame(
        sample_data(physeq),
        row.names = sample_names(physeq), 
        stringsAsFactors = FALSE, 
        check.names = FALSE
      )
      
      #check if any columns match exactly with rownames
      #if none matched assume row names are sample identifiers
      samplesCol <- unlist(lapply(metadata, function(x) {
        identical(x, rownames(metadata))}))
      
      if(any(samplesCol)) {
        #error if a column matched and it's not the first
        if(!samplesCol[[1]])
          stop("Sample ID's must be in the first column in the sample metadata, please reorder", call. = FALSE)
      } else {
        #assume rownames are sample identifiers, merge at the end with name "SampleID"
        if(any(colnames(metadata) %in% "SampleID"))
          stop("A column in the sample metadata is already named \"SampleID\" but does not seem to contain sample ID's", call. = FALSE)
        metadata$SampleID <- rownames(metadata)
        
        #reorder columns so SampleID is the first
        metadata <- metadata[, c(which(colnames(metadata) %in% "SampleID"), 1:(ncol(metadata)-1L)), drop = FALSE]
      }
    } else
      metadata <- NULL
    
    #extract phylogenetic tree, assumed to be of class "phylo"
    if(!is.null(physeq@phy_tree)) {
      tree <- phy_tree(physeq)
    } else
      tree <- NULL
    
    #extract OTU DNA sequences, assumed to be of class "XStringSet"
    if(!is.null(physeq@refseq)) {
      checkReqPkg(
        "Biostrings",
        " Please install with:\n  install.packages(\"BiocManager\"); BiocManager::install(\"Biostrings\")"
      )
      #convert XStringSet to DNAbin using a temporary file (easiest)
      fasta <- tempfile(pattern = "ampvis2_", fileext = ".fa")
      writeXStringSet(physeq@refseq, filepath = fasta)
    } else
      fasta <- NULL
    
  }
  
  ### import and check otutable (with or without taxonomy)
  otutable <- import(otutable, ...)
  otutable <- findOTUcol(otutable, OTUcolname = otutable_OTUcolname)

  ### extract read abundances from otutable (same if taxonomy present or not)
  taxcols <- tolower(colnames(otutable)) %in% c("domain", tolower(tax_levels))
  abund <- otutable[, !taxcols, drop = FALSE]
  abund[is.na(abund)] <- 0L

  if(!all(sapply(abund, is.numeric)))
    stop("Abundance table contains non-numeric columns", call. = FALSE)

  # warn if any empty sample(s)
  if (any(colSums(abund) == 0L)) {
    warning("One or more sample(s) have 0 reads", call. = FALSE)
  }

  ### extract taxonomy from otutable or separate table if provided
  if (is.null(taxonomy)) {
    tax <- parseTaxonomy(otutable)
  } else if (!is.null(taxonomy)) {
    taxonomy <- import(taxonomy)
    taxonomy <- findOTUcol(taxonomy, OTUcolname = taxonomy_OTUcolname)
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

  ### create dummy metadata if not provided or passed on as an attribute of otutable
  if (!is.null(metadata)) {
    metadata <- import(metadata)
  } else if (!is.null(attr(otutable, "metadata"))) {
    # if otutable is loaded from a biom file, sample metadata is also there
    metadata <- attr(otutable, "metadata")
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

  abund0 <- abund
  nOTUsbefore <- nrow(abund) #used if removeAbsentOTUs = TRUE
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
  
  if (isTRUE(removeAbsentOTUs)) {
    abund0 <- abund0[rowSums(abund0) > 0, , drop = FALSE]
    nOTUsremoved <- nOTUsbefore - nrow(abund0)
    if (nOTUsremoved > 0) {
      message(
        nOTUsremoved,
        if (nOTUsremoved == 1L) {
          " OTU"
        } else {
          " OTUs"
        },
        " with 0 total abundance across all samples",
        if (nOTUsremoved == 1L) {
          " has"
        } else {
          " have"
        },
        " been removed."
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
  # check file
  if (!is.null(fasta)) {
    if (is.character(fasta) &
      length(fasta) == 1 &
      is.null(dim(fasta))) {
      refseq <- ape::read.FASTA(unzip_file(fasta), ...)[data$tax$OTU]
    } else if (!inherits(fasta, c("DNAbin", "AAbin"))) {
      stop("fasta must be of class \"DNAbin\" or \"AAbin\" as loaded with the ape::read.FASTA() function.", call. = FALSE)
    } else if(inherits(fasta, c("DNAbin", "AAbin"))) {
      refseq <- fasta[data$tax$OTU]
    }
    if (all(sapply(refseq, is.null))) {
      stop("No sequences match any OTU's", call. = FALSE)
    }
    data$refseq <- refseq
  }

  ### append phylogenetic tree if provided
  # check tree
  if (!is.null(tree)) {
    if (is.character(tree) &
      length(tree) == 1 &
      is.null(dim(tree))) {
      tree <- ape::read.tree(unzip_file(tree), ...)
    } else if (!inherits(tree, "phylo")) {
      stop("The provided phylogenetic tree must be of class \"phylo\" as loaded with the ape::read.tree() function.", call. = FALSE)
    }
    tree <- ape::drop.tip(
      phy = tree,
      tip = tree$tip.label[!tree$tip.label %in% data$tax$OTU]
    )
    if (is.null(tree)) {
      stop("No tree tip labels match any OTU's", call. = FALSE)
    }
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
