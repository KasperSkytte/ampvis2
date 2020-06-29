#' Prints heatmap+functions plot (internal function)
#'
#' @param x plot
#' @param ... not used
#'
#' @export
#' @importFrom patchwork plot_layout
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
print.hmfunplot <- function(x, ...) {
  p <- x$heatmap +
    x$functions +
    patchwork::plot_layout(
      ncol = 2,
      widths = attributes(x)[["rel_widths"]],
      guides = "collect"
    )
  print(p)
}
}

#' Prints ampvis2 object summary (internal function)
#'
#' @param x (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param ... not used
#'
#' @importFrom crayon underline
#' @export
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
print.ampvis2 <- function(x, ...) {
  ### calculate basic statistics and useful information about the data, print it
  abundAreCounts <- abundAreCounts(x)
  if (isTRUE(abundAreCounts)) {
    # calculate basic stats of reads
    readstats <- list(
      "Total#Reads" = as.character(sum(x$abund)),
      "Min#Reads" = as.character(min(colSums(x$abund))),
      "Max#Reads" = as.character(max(colSums(x$abund))),
      "Median#Reads" = as.character(median(colSums(x$abund))),
      "Avg#Reads" = as.character(round(mean(colSums(x$abund)), digits = 2))
    )
  } else if (isFALSE(abundAreCounts)) {
    readstats <- NULL
  }
  cat(class(x), "object with", length(x), "elements.", crayon::underline("\nSummary of OTU table:\n"))
  print.table(c(
    "Samples" = as.character(ncol(x$abund)),
    "OTUs" = as.character(nrow(x$abund)),
    readstats
  ),
  justify = "right"
  )
  if (isFALSE(abundAreCounts)) {
    cat("(The read counts have been normalised)\n")
  }
  cat(crayon::underline("\nAssigned taxonomy:\n"))
  print.table(c(
    "Kingdom" = paste0(sum(nchar(x$tax$Kingdom) > 3), "(", round(sum(nchar(x$tax$Kingdom) > 3) / nrow(x$abund), digits = 2) * 100, "%)"),
    "Phylum" = paste0(sum(nchar(x$tax$Phylum) > 3), "(", round(sum(nchar(x$tax$Phylum) > 3) / nrow(x$abund) * 100, digits = 2), "%)"),
    "Class" = paste0(sum(nchar(x$tax$Class) > 3), "(", round(sum(nchar(x$tax$Class) > 3) / nrow(x$abund) * 100, digits = 2), "%)"),
    "Order" = paste0(sum(nchar(x$tax$Order) > 3), "(", round(sum(nchar(x$tax$Order) > 3) / nrow(x$abund) * 100, digits = 2), "%)"),
    "Family" = paste0(sum(nchar(x$tax$Family) > 3), "(", round(sum(nchar(x$tax$Family) > 3) / nrow(x$abund) * 100, digits = 2), "%)"),
    "Genus" = paste0(sum(nchar(x$tax$Genus) > 3), "(", round(sum(nchar(x$tax$Genus) > 3) / nrow(x$abund) * 100, digits = 2), "%)"),
    "Species" = paste0(sum(nchar(x$tax$Species) > 3), "(", round(sum(nchar(x$tax$Species) > 3) / nrow(x$abund) * 100, digits = 2), "%)")
  ),
  justify = "right"
  )
  cat(crayon::underline("\nMetadata variables:"), as.character(ncol(x$metadata)), "\n", paste(as.character(colnames(x$metadata)), collapse = ", "))
  cat("\n")
}

#' Print method for figure caption created by amp_ordinate
#'
#' @param x Character vector with the caption
#' @param ... not used
#'
#' @importFrom cli cat_line rule
#' @importFrom crayon italic
#' @export
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
print.figcaption <- function(x, ...) {
  cli::cat_line(cli::rule("Auto-generated figure caption (start)"))
  x %>%
    strwrap() %>%
    crayon::italic() %>%
    cli::cat_line()
  cli::cat_line(cli::rule("Auto-generated figure caption (end)"))
}
