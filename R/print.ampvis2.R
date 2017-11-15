#' Prints ampvis2 object summary. Internal function.
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#'
#' @return Text output in console
#' @export
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' 
print.ampvis2 <- function(data) {
  ###calculate basic statistics and useful information about the data, print it
  cat(class(data), "object with", length(data),"elements.\nSummary of OTU table:\n")
  print.table(c("Samples" = as.character(ncol(data$abund)),
    "OTUs" = as.character(nrow(data$abund)),
    "Total#Reads" = as.character(sum(data$abund)),
    "Min#Reads" = as.character(min(colSums(data$abund))),
    "Max#Reads" = as.character(max(colSums(data$abund))),
    "Median#Reads" = as.character(median(colSums(data$abund))),
    "Avg#Reads" = as.character(round(mean(colSums(data$abund)), digits = 2))), 
    justify = "right")
  cat("\nAssigned taxonomy:\n")
  print.table(c("Kingdom" = paste(sum(nchar(data$tax$Kingdom) > 3), "(", round(sum(nchar(data$tax$Kingdom) > 3) / nrow(data$abund), digits = 2) * 100, "%)", sep = ""),
    "Phylum" = paste(sum(nchar(data$tax$Phylum) > 3), "(", round(sum(nchar(data$tax$Phylum) > 3) / nrow(data$abund) * 100, digits = 2), "%)", sep = ""),
    "Class" = paste(sum(nchar(data$tax$Class) > 3), "(", round(sum(nchar(data$tax$Class) > 3) / nrow(data$abund) * 100, digits = 2), "%)", sep = ""),
    "Order" = paste(sum(nchar(data$tax$Order) > 3), "(", round(sum(nchar(data$tax$Order) > 3) / nrow(data$abund) * 100, digits = 2), "%)", sep = ""),
    "Family" = paste(sum(nchar(data$tax$Family) > 3), "(", round(sum(nchar(data$tax$Family) > 3) / nrow(data$abund) * 100, digits = 2), "%)", sep = ""),
    "Genus" = paste(sum(nchar(data$tax$Genus) > 3), "(", round(sum(nchar(data$tax$Genus) > 3) / nrow(data$abund) * 100, digits = 2), "%)", sep = ""),
    "Species" = paste(sum(nchar(data$tax$Species) > 3), "(", round(sum(nchar(data$tax$Species) > 3) / nrow(data$abund) * 100, digits = 2), "%)", sep = "")),
    justify = "right")
  cat("\nMetadata variables:", as.character(ncol(data$metadata)), "\n", paste(as.character(colnames(data$metadata)), collapse = ", "))
}