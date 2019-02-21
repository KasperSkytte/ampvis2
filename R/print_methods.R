#' Prints heatmap+functions plot (internal function)
#'
#' @param hmfunplot plot
#'
#' @export
#' @importFrom cowplot plot_grid
print.hmfunplot <- function(hmfunplot) {
  print(cowplot::plot_grid(hmfunplot$heatmap,
    hmfunplot$functions,
    ncol = 2,
    rel_widths = attributes(hmfunplot)[["rel_widths"]],
    align = "h",
    axis = "tb"
  ))
}

#' Prints ampvis2 object summary (internal function)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#'
#' @importFrom crayon underline
#' @export
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
print.ampvis2 <- function(data) {
  ### calculate basic statistics and useful information about the data, print it
  if (!isTRUE(attributes(data)$normalised)) {
    # calculate basic stats and store in attributes for use in print.ampvis2
    readstats <- attributes(data)$readstats <- list(
      "Total#Reads" = as.character(sum(data$abund)),
      "Min#Reads" = as.character(min(colSums(data$abund))),
      "Max#Reads" = as.character(max(colSums(data$abund))),
      "Median#Reads" = as.character(median(colSums(data$abund))),
      "Avg#Reads" = as.character(round(mean(colSums(data$abund)), digits = 2))
    )
  } else if (isTRUE(attributes(data)$normalised)) {
    readstats <- attributes(data)$readstats
  }
  cat(class(data), "object with", length(data), "elements.", crayon::underline("\nSummary of OTU table:\n"))
  print.table(c(
    "Samples" = as.character(ncol(data$abund)),
    "OTUs" = as.character(nrow(data$abund)),
    readstats
  ),
  justify = "right"
  )
  if (isTRUE(attributes(data)$normalised)) {
    cat("(The read counts have been normalised)\n")
  }
  cat(crayon::underline("\nAssigned taxonomy:\n"))
  print.table(c(
    "Kingdom" = paste0(sum(nchar(data$tax$Kingdom) > 3), "(", round(sum(nchar(data$tax$Kingdom) > 3) / nrow(data$abund), digits = 2) * 100, "%)"),
    "Phylum" = paste0(sum(nchar(data$tax$Phylum) > 3), "(", round(sum(nchar(data$tax$Phylum) > 3) / nrow(data$abund) * 100, digits = 2), "%)"),
    "Class" = paste0(sum(nchar(data$tax$Class) > 3), "(", round(sum(nchar(data$tax$Class) > 3) / nrow(data$abund) * 100, digits = 2), "%)"),
    "Order" = paste0(sum(nchar(data$tax$Order) > 3), "(", round(sum(nchar(data$tax$Order) > 3) / nrow(data$abund) * 100, digits = 2), "%)"),
    "Family" = paste0(sum(nchar(data$tax$Family) > 3), "(", round(sum(nchar(data$tax$Family) > 3) / nrow(data$abund) * 100, digits = 2), "%)"),
    "Genus" = paste0(sum(nchar(data$tax$Genus) > 3), "(", round(sum(nchar(data$tax$Genus) > 3) / nrow(data$abund) * 100, digits = 2), "%)"),
    "Species" = paste0(sum(nchar(data$tax$Species) > 3), "(", round(sum(nchar(data$tax$Species) > 3) / nrow(data$abund) * 100, digits = 2), "%)")
  ),
  justify = "right"
  )
  cat(crayon::underline("\nMetadata variables:"), as.character(ncol(data$metadata)), "\n", paste(as.character(colnames(data$metadata)), collapse = ", "))
}

#' Print method for figure caption created by amp_ordinate
#'
#' @param captionwithrefs Character vector with the caption
#'
#' @importFrom cli cat_line rule
#' @importFrom crayon italic
#' @export
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
print.figcaption <- function(captionwithrefs) {
  cli::cat_line(cli::rule("Auto-generated figure caption (start)"))
  captionwithrefs %>%
    strwrap() %>%
    crayon::italic() %>%
    cli::cat_line()
  cli::cat_line(cli::rule("Auto-generated figure caption (end)"))
}
