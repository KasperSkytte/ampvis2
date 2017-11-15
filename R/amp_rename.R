#' Tidy up taxonomy
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
#' @import dplyr
#' @export
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_rename <- function(data, tax_class = NULL, tax_empty = "best", tax_level = "Genus"){
  
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
      mutate(Class, Class = ifelse(Class == "", ifelse(grepl("__", Phylum), Phylum, paste("c__", Phylum, "_", rownames(tax), sep = "")), Class)) %>%
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
