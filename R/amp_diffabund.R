#' Differential abundance test
#'
#' Tests if there is a significant difference in abundances between samples or groups hereof based on selected conditions. Returns a list containing test results as well as two different plots; an MA-plot and an abundance plot with taxa with the most significant p-value (below the threshold).
#'
#' @usage amp_diffabund(data, group)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param group (required) A categorical variable in the metadata that defines the sample groups to test. 
#' @param test The name of the test to use, either \code{"Wald} or \code{"LRT}. See \code{\link[DESeq2]{DESeq}}. (\emph{default:} \code{"Wald"})
#' @param fitType The type of fitting of dispersions to the mean intensity, either \code{"parametric"}, \code{"local"}, or \code{"mean"}. (\emph{default:} \code{"parametric"})
#' @param verbose (\emph{Logical}) Whether to print status messages during the test calculations. (\emph{Default: } \code{TRUE}) 
#' @param num_threads The number of threads to use for parallelization by the \href{https://www.bioconductor.org/packages/devel/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.pdf}{\code{BiocParallel}} backend. Parallelization is not supported on windows machines. (\emph{default:} \code{1})
#' @param signif_plot_type Either \code{"boxplot"} or \code{"point"}. (\emph{default:} \code{"point"})
#' @param signif_thrh Significance threshold. (\emph{default:} \code{0.01})
#' @param fold Log2fold filter for displaying significant results. (\emph{default:} \code{0})
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"Phylum"})
#' @param tax_add Additional taxonomic level(s) to display, e.g. \code{"Phylum"}. (\emph{default:} \code{"none"})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible. 
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param plot_nshow The amount of the most significant results to display in the most-significant plot. (\emph{default:} \code{10})
#' @param plot_point_size The size of the plotted points. (\emph{default:} \code{2})
#' @param adjust_zero Keep 0 abundances in ggplot2 median calculations by adding a small constant to these.
#'
#' @return A list with multiple elements:
#'   \itemize{
#'     \item \code{"DESeq2_results"}: The raw output result from \code{\link[DESeq2]{DESeq}}.
#'     \item \code{"DESeq2_results_signif"}: The raw output result from \code{\link[DESeq2]{DESeq}}, but subset to only taxa with p-value below the threshold set by \code{signif_thrh}.
#'     \item \code{"signif_plotdata"}: The data used to generate the ggplots, but subset to only taxa with p-value below the threshold set by \code{signif_thrh}.
#'     \item \code{"Clean_results"}: A simpler version of DESeq2_results_signif only with adjusted p-values, log2FoldChange, and average abundance of each taxa per \code{group}.
#'     \item \code{"plot_MA"}: MA-plot
#'     \item \code{"plot_MA_plotly"}: Interactive \code{plotly} plot of \code{MA-plot} with custom hover information.
#'     \item \code{"plot_signif"}: Abundance plot with taxa with the n most significant p-value (below the threshold), where n is set by \code{plot_nshow}.
#'     \item \code{"plot_signif_plotly"}: Interactive \code{plotly} plot of \code{plot_signif} with custom hover information. 
#'   }
#' 
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix results
#' @importFrom dplyr filter arrange group_by mutate select summarise
#' @importFrom data.table as.data.table setkey
#' @importFrom stringr str_split
#' @importFrom purrr imap
#' @importFrom plotly ggplotly
#' @importFrom tidyr gather spread unite
#' 
#' @export
#' @examples
#' #Load example data
#' data("AalborgWWTPs")
#'
#' #Save the results in an object
#' results <- amp_diffabund(AalborgWWTPs, group = "Plant", num_threads = 4)
#'
#' #Show plots
#' results$plot_signif
#' results$plot_MA
#'
#' #Or show raw results
#' results$Clean_results
#'
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
amp_diffabund <- function(data,
                          group,
                          test = "Wald",
                          fitType = "parametric",
                          num_threads = 1L,
                          signif_thrh = 0.01,
                          fold = 0,
                          verbose = TRUE,
                          signif_plot_type = "point",
                          plot_nshow = 10,
                          plot_point_size = 2,
                          tax_aggregate = "OTU",
                          tax_add = NULL,
                          tax_class = NULL,
                          tax_empty = "best",
                          adjust_zero = NULL) {
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)", call. = FALSE)
  
  ## Clean up the taxonomy
  data <- ampvis2:::amp_rename(data = data,
                               tax_class = tax_class,
                               tax_empty = tax_empty,
                               tax_level = tax_aggregate)
  
  #tax_add and tax_aggregate can't be the same
  if(!is.null(tax_aggregate) & !is.null(tax_add)) {
    if(identical(tax_aggregate, tax_add)) {
      stop("tax_aggregate and tax_add cannot be the same", call. = FALSE)
    }
  }
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  metadata <- data[["metadata"]]
  
  # fix group factors to be syntactically valid
  metadata[,group] %<>% stringr::str_replace_all("[^[:alnum:]_.]", "_") %>% as.factor()
  
  ## Make a name variable that can be used instead of tax_aggregate to display multiple levels
  suppressWarnings(
    if (!is.null(tax_add)){
      if (tax_add != tax_aggregate) {
        tax <- data.frame(tax, Display = apply(tax[,c(tax_add,tax_aggregate)], 1, paste, collapse="; "))
      }
    } else {
      tax <- data.frame(tax, Display = tax[,tax_aggregate])
    }
  )
  
  # Aggregate to a specific taxonomic level
  abund3 <- cbind.data.frame(Display = tax[,"Display"], abund) %>%
    tidyr::gather(key = Sample, value = Abundance, -Display) %>% as.data.table()
  
  abund3 <- abund3[, "sum":=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    as.data.frame() %>%
    select(-Abundance) %>%
    unique()
  
  ##### Convert to DESeq2 format and test for significant differential abundance #####
  abund4 <- tidyr::spread(data = abund3, key = Sample, value = sum) %>%
    column_to_rownames("Display")
  abund4 <- abund4[,metadata[[1]]]
  
  if(isTRUE(verbose))
    message("Running DESeq2 differential abundance test. This may take a while depending on the size of the data. \n---------------------------------")
  data_deseq <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = abund4, 
                                                                colData = metadata, 
                                                                design = as.formula(paste("~", group, sep=""))))
  
  data_deseq_test = DESeq2::DESeq(data_deseq, 
                                  test = test,
                                  fitType = fitType, 
                                  quiet = if(!isTRUE(verbose)) TRUE else FALSE,
                                  parallel = if(num_threads > 1L) TRUE else FALSE,
                                  BPPARAM = BiocParallel::MulticoreParam(num_threads))
  
  ## Extract the results
  res = DESeq2::results(data_deseq_test, 
                        cooksCutoff = FALSE
                        #,parallel = if(num_threads > 1L) TRUE else FALSE
                        #,BPPARAM = BiocParallel::MulticoreParam(num_threads)
                        )
  
  res_tax = data.frame(as.data.frame(res), Tax = rownames(res))
  
  res_tax_sig = filter(res_tax, padj < signif_thrh & fold < abs(log2FoldChange)) %>%
    arrange(padj)
  
  if(isTRUE(verbose))
    message("---------------------------------\nDone. Generating plots.")
  
  ##### MA plot #####
  res_tax$Significant <- ifelse(rownames(res_tax) %in% res_tax_sig$Tax ,
                                "Significant",
                                "Not significant")
  res_tax$Significant[is.na(res_tax$Significant)] <- "Not significant"
  
  MAplot <- ggplot(data = res_tax, 
                   aes(x = baseMean, 
                       y = log2FoldChange,
                       color = Significant)) +
    scale_x_log10() +
    scale_color_manual(values=c("black", "red")) +
    labs(x = "BaseMean read abundance", y = "Log2 fold change") +
    theme_classic() +
    theme(panel.grid.major.x = element_line(color = "grey90"),
          panel.grid.major.y = element_line(color = "grey90"),
          legend.title = element_blank())
  
  #find the lowest taxonomic level and generate a character vector with taxonomic
  #information for each taxa (only with taxonomy from the lowest taxonomic level and up)
  taxlevels <- factor(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
                      c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  if(!is.null(tax_add)) {
    lowestlevel <- as.character(taxlevels[max(as.numeric(c(taxlevels[which(taxlevels %in% tax_aggregate)], 
                                                           taxlevels[which(taxlevels %in% tax_add)])))])
  } else
    lowestlevel <- tax_aggregate
  
  data_plotly <- data$tax %>% 
    .[which(.[,which(colnames(.) == lowestlevel)] %in% unlist(stringr::str_split(res_tax$Tax, "; "))),] %>%
    .[,1:which(colnames(.) == lowestlevel)] %>%
    purrr::imap(~paste(.y, .x, sep = ": ")) %>%
    as.data.frame() %>% 
    tidyr::unite("test", sep = "<br>") %>%
    unlist(use.names = FALSE) %>%
    unique()
  
  ##### data for significance plot #####
  abund5 <- mutate(abund4, Tax = rownames(abund4)) %>%
    tidyr::gather(key=Sample, value=Count, -Tax) %>%
    group_by(Sample) %>%
    mutate(Abundance = Count / sum(Count)*100)
  
  abund6 <- suppressWarnings(dplyr::inner_join(abund5, res_tax, by = "Tax")) %>%
    filter(padj < signif_thrh & fold < abs(log2FoldChange)) %>%
    arrange(padj)
  
  if(nrow(abund6) == 0){stop("No significant differences found.", call. = FALSE)}
  
  if(!is.null(adjust_zero)){
    abund6$Abundance[abund6$Abundance==0] <- adjust_zero
  }
  
  colnames(metadata)[1] <- "Sample"
  metadata <- metadata[c("Sample",group)]
  colnames(metadata)[2] <- "Group"
  
  point_df <- dplyr::inner_join(x = abund6, y = metadata, by = "Sample") %>%
    group_by(Sample) %>%
    arrange(padj)
  
  colnames(point_df)[12] <- group
  
  if(!is.null(plot_nshow)){
    if(plot_nshow < nrow(abund6)){plot_nshow <- nrow(abund6)}
    point_df <- filter(point_df, Tax %in% as.character(unique(point_df$Tax))[1:plot_nshow])
  }
  
  point_df$Tax <- factor(point_df$Tax, levels = rev(as.character(unique(point_df$Tax))[1:plot_nshow]))
  
  ##### significance plot #####
  signifplot <- ggplot(data = point_df, aes_string(x = "Tax", y = "Abundance", color = group)) +
    labs(x = "", y = "Read Abundance (%)") +
    coord_flip() +
    theme_classic() +
    theme(panel.grid.major.x = element_line(color = "grey90"),
          panel.grid.major.y = element_line(color = "grey90"))
  
  if(signif_plot_type == "point") {
    signifplot <- signifplot + geom_jitter(position = position_jitter(width = .05), size = plot_point_size)
  } else if(signif_plot_type == "boxplot") {
    signifplot <- signifplot + geom_boxplot(outlier.size=1)
  }
  
  ##### return results #####
  clean_res0 <- suppressWarnings(dplyr::inner_join(abund5, res_tax, by = "Tax")) %>%
    dplyr::inner_join(y = metadata, by = "Sample") %>%
    group_by(Sample) %>%
    arrange(padj)
  
  colnames(clean_res0)[12] <- "Group"
  
  cr <- mutate(clean_res0,
               padj = signif(padj, 2),
               Log2FC = signif(log2FoldChange, 2),
               Taxonomy = Tax) %>%
    group_by(Group, Taxonomy, padj, Log2FC) %>%
    summarise(Avg = round(mean(Abundance), 3)) %>%
    tidyr::spread(key = Group, value = Avg) %>%
    arrange(padj)
  
  out <- list(DESeq2_results = res, 
              DESeq2_results_signif = res_tax_sig, 
              signif_plotdata = point_df, 
              Clean_results = cr,
              plot_MA = MAplot + 
                geom_point(size = plot_point_size),
              plot_MA_plotly = plotly::ggplotly(MAplot + 
                                                  suppressWarnings(geom_point(size = plot_point_size-1,
                                                             aes(text = data_plotly))),
                                                tooltip = "text"),
              plot_signif = signifplot,
              plot_signif_plotly = plotly::ggplotly(signifplot))
  return(out)
}
