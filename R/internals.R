#' @title Rarefy ampvis2 object (internal function)
#' @description This is just a wrapper of \code{\link[vegan]{rrarefy}} with convenient error messages and adjusted to work with ampvis2 objects.
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param rarefy (\emph{required}) Passed directly to \code{\link[vegan]{rrarefy}}.
#' 
#' @return An ampvis2 object with rarefied OTU abundances.
#' 
#' @importFrom magrittr %>%
#' @importFrom vegan rrarefy
amp_rarefy <- function(data, rarefy) {
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis2 functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)", call. = FALSE)
  
  set.seed(0) #for reproducibility
  
  reads <- colSums(data$abund)
  maxreads <- max(reads)
  minreads <- min(reads)
  if(is.numeric(rarefy)){
    if(rarefy > maxreads ) {
      stop("The chosen rarefy size is larger than the largest amount of reads in any sample (", as.character(maxreads), ").", call. = FALSE)
    } else if (rarefy < minreads) {
      data$abund <- suppressWarnings(vegan::rrarefy(t(data$abund), sample = rarefy)) %>% t() %>% as.data.frame() 
      warning("The chosen rarefy size (", as.character(rarefy), ") is smaller than the smallest amount of reads in any sample (", as.character(minreads), ").", call. = FALSE)
    } else {
      data$abund <- suppressWarnings(vegan::rrarefy(t(data$abund), sample = rarefy)) %>% t() %>% as.data.frame()
      if (minreads < rarefy) {
        message("The following sample(s) have not been rarefied (less than ", as.character(rarefy), " reads):\n", paste(rownames(data$metadata[which(reads < rarefy),]), collapse = ", "))
      }
    }
  } else if(!is.numeric(rarefy)) {
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
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
amp_rename <- function(data, tax_class = NULL, tax_empty = "best", tax_level = "Genus") {
  
  tax = data[["tax"]]
  
  ## First make sure that all entries are strings
  for ( i in 1:ncol(tax) ){
    tax[,i] <- as.character(tax[,i])  
  }
  
  ## Change a specific phylum to class level
  if(!is.null(tax_class)){
    for (i in 1:nrow(tax)){
      if (!is.na(tax$Phylum[i]) & tax$Phylum[i] %in% tax_class){
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
  if (!is.null(tax$Species)){
    tax$Species <- gsub("^s_*", "", tax$Species)
  } else {
    tax$Species <- ""
  }
  
  tax[is.na(tax)] <- ""
  
  ## How to handle empty taxonomic assignments
  if (tax_empty == "OTU"){
    for (i in 1:nrow(tax)) {
      if (tax[i,"Species"] == "") {tax[i,"Species"] <- rownames(tax)[i]}
      if (tax[i,"Genus"] == "") {tax[i,"Genus"] <- rownames(tax)[i]}
      if (tax[i,"Family"] == "") {tax[i,"Family"] <- rownames(tax)[i]}
      if (tax[i,"Order"] == "") {tax[i,"Order"] <- rownames(tax)[i]}
      if (tax[i,"Class"] == "") {tax[i,"Class"] <- rownames(tax)[i]}
      if (tax[i,"Phylum"] == "") {tax[i,"Phylum"] <- rownames(tax)[i]}
    }
  }
  
  ## Handle empty taxonomic strings
  rn <- rownames(tax) #damn rownames are silently dropped by mutate()
  if(tax_empty == "best"){
    tax <- mutate(tax, Kingdom, Kingdom = ifelse(Kingdom == "", "Unclassified", Kingdom)) %>%
      mutate(Phylum, Phylum = ifelse(Phylum == "", paste("k__", Kingdom, "_", rownames(tax), sep = ""), Phylum)) %>%
      mutate(Class, Class = ifelse(Class == "", ifelse(grepl("__", Phylum), Phylum, paste("p__", Phylum, "_", rownames(tax), sep = "")), Class)) %>%
      mutate(Order, Order = ifelse(Order == "", ifelse(grepl("__", Class), Class, paste("c__", Class, "_", rownames(tax), sep = "")), Order)) %>%
      mutate(Family, Family = ifelse(Family == "", ifelse(grepl("__", Order), Order, paste("o__", Order, "_", rownames(tax), sep = "")), Family)) %>%
      mutate(Genus, Genus = ifelse(Genus == "", ifelse(grepl("__", Family), Family, paste("f__", Family, "_", rownames(tax), sep = "")), Genus)) %>%
      mutate(Species, Species = ifelse(Species == "", ifelse(grepl("__", Genus), Genus, paste("g__", Genus, "_", rownames(tax), sep = "")), Species))
  }
  rownames(tax) <- rn
  
  if(tax_empty == "remove"){
    abund <- data[["abund"]]
    tax <- subset(tax, tax[,tax_level] != "")
    abund <- subset(abund, rownames(abund) %in% rownames(tax))
    data[["abund"]] <- abund
  }
  data[["tax"]] <- tax
  rownames(data[["tax"]]) <- rownames(tax)
  
  return(data)
}
#' Prints heatmap+functions plot (internal function)
#'
#' @param hmfunplot plot
#'
#' @export
#' @importFrom cowplot plot_grid
#' 
print.hmfunplot <- function(hmfunplot) {
  print(cowplot::plot_grid(hmfunplot$heatmap,
                     hmfunplot$functions,
                     ncol = 2,
                     rel_widths = attributes(hmfunplot)[["rel_widths"]],
                     align = "h",
                     axis = "tb"))
}
#' Prints ampvis2 object summary (internal function)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#'
#' @importFrom crayon underline
#' @export
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' 
print.ampvis2 <- function(data) {
  ###calculate basic statistics and useful information about the data, print it
  if(!isTRUE(attributes(data)$normalised)) {
    #calculate basic stats and store in attributes for use in print.ampvis2
    readstats <- attributes(data)$readstats <- list(
      "Total#Reads" = as.character(sum(data$abund)),
      "Min#Reads" = as.character(min(colSums(data$abund))),
      "Max#Reads" = as.character(max(colSums(data$abund))),
      "Median#Reads" = as.character(median(colSums(data$abund))),
      "Avg#Reads" = as.character(round(mean(colSums(data$abund)), digits = 2))
    )
  } else if(isTRUE(attributes(data)$normalised)) {
    readstats <- attributes(data)$readstats
  }
  cat(class(data), "object with", length(data),"elements.", crayon::underline("\nSummary of OTU table:\n"))
  print.table(c("Samples" = as.character(ncol(data$abund)),
                "OTUs" = as.character(nrow(data$abund)),
                readstats), 
              justify = "right")
  if(isTRUE(attributes(data)$normalised))
    cat("(The read counts have been normalised)\n")
  cat(crayon::underline("\nAssigned taxonomy:\n"))
  print.table(c("Kingdom" = paste0(sum(nchar(data$tax$Kingdom) > 3), "(", round(sum(nchar(data$tax$Kingdom) > 3) / nrow(data$abund), digits = 2) * 100, "%)"),
                "Phylum" = paste0(sum(nchar(data$tax$Phylum) > 3), "(", round(sum(nchar(data$tax$Phylum) > 3) / nrow(data$abund) * 100, digits = 2), "%)"),
                "Class" = paste0(sum(nchar(data$tax$Class) > 3), "(", round(sum(nchar(data$tax$Class) > 3) / nrow(data$abund) * 100, digits = 2), "%)"),
                "Order" = paste0(sum(nchar(data$tax$Order) > 3), "(", round(sum(nchar(data$tax$Order) > 3) / nrow(data$abund) * 100, digits = 2), "%)"),
                "Family" = paste0(sum(nchar(data$tax$Family) > 3), "(", round(sum(nchar(data$tax$Family) > 3) / nrow(data$abund) * 100, digits = 2), "%)"),
                "Genus" = paste0(sum(nchar(data$tax$Genus) > 3), "(", round(sum(nchar(data$tax$Genus) > 3) / nrow(data$abund) * 100, digits = 2), "%)"),
                "Species" = paste0(sum(nchar(data$tax$Species) > 3), "(", round(sum(nchar(data$tax$Species) > 3) / nrow(data$abund) * 100, digits = 2), "%)")),
              justify = "right")
  cat(crayon::underline("\nMetadata variables:"), as.character(ncol(data$metadata)), "\n", paste(as.character(colnames(data$metadata)), collapse = ", "))
}