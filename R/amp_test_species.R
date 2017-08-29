#' Test species abundances
#'
#' Tests if there is a significant difference in abundances between samples or groups hereof based on selected conditions. 
#'
#' @usage amp_test_species(data, group)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param group (required) The group to test.
#' @param signif_thrh Significance treshold. (\emph{default:} \code{0.01})
#' @param fold Log2fold filter for displaying significant results. (\emph{default:} \code{0})
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"OTU"})
#' @param tax_add Additional taxonomic level(s) to display, e.g. \code{"Phylum"}. (\emph{default:} \code{"none"})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible. 
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param plot_type Either \code{"boxplot"} or \code{"point"}. (\emph{default:} \code{"point"})
#' @param plot_nshow The amount of the most significant results to display. (\emph{default:} \code{10})
#' @param plot_point_size The size of the plotted points. (\emph{default:} \code{2})
#' @param adjust_zero Keep 0 abundances in ggplot2 median calculations by adding a small constant to these.
#' @param plotly Returns an interactive plot instead. (\emph{default:} \code{FALSE})
#' 
#' @return A list with multiple elements. 
#' @import dplyr
#' @import DESeq2
#' @import ggplot2
#' @import tidyr
#' @import plotly
#' @import data.table
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_test_species <- function(data,
                             group,
                             test = "Wald",
                             fitType = "parametric",
                             signif_thrh = 0.01,
                             fold = 0,
                             label = FALSE,
                             plot_type = "point",
                             plot_nshow = 10,
                             plot_point_size = 2,
                             tax_aggregate = "OTU",
                             tax_add = NULL,
                             tax_class = NULL,
                             tax_empty = "best",
                             adjust_zero = NULL,
                             plotly = FALSE){
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax_class = tax_class, tax_empty = tax_empty, tax_level = tax_aggregate)
  
  #tax_add and tax_aggregate can't be the same
  if(!is.null(tax_aggregate) & !is.null(tax_add)) {
    if(tax_aggregate == tax_add) {
      stop("tax_aggregate and tax_add cannot be the same")
    }
  }
  
  if (is.null(group)) {
    stop("Argument 'group' must be provided.")
  }
  
  ## Extract the data into seperate objects for readability
  abund <- data[["abund"]]  
  tax <- data[["tax"]]
  metadata <- data[["metadata"]]

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
    gather(key = Sample, value = Abundance, -Display) %>% as.data.table()
  
  abund3 <- abund3[, "sum":=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    unique() %>% 
    select(-Abundance)
  
  ## Convert to DESeq2 format
  abund4 <- spread(data = abund3, key = Sample, value = sum)
  rownames(abund4) <- abund4$Display
  abund4 <- abund4[,-1]
  
  groupF <- as.formula(paste("~", group, sep=""))
  
  
  data_deseq <- DESeqDataSetFromMatrix(countData = abund4,
                                       colData = metadata,
                                       design = groupF)
  
  ## Test for significant differential abundance
  data_deseq_test = DESeq(data_deseq, test=test, fitType=fitType)
  
  ## Extract the results
  res = results(data_deseq_test, cooksCutoff = FALSE)  
  res_tax = data.frame(as.data.frame(res), Tax = rownames(res))
  
  res_tax_sig = subset(res_tax, padj < signif_thrh & fold < abs(log2FoldChange)) %>%
    arrange(padj)
  
  ## Plot the data
  ### MA plot
  res_tax$Significant <- ifelse(rownames(res_tax) %in% res_tax_sig$Tax , "Yes", "No")
  res_tax$Significant[is.na(res_tax$Significant)] <- "No"
  
  
  
  
  p1 <- ggplot(data = res_tax, aes(x = baseMean, y = log2FoldChange, color = Significant)) + 
    scale_x_log10() +
    scale_color_manual(values=c("black", "red")) +
    labs(x = "BaseMean read abundance", y = "Log2 fold change") +
    theme_classic() +
    theme(panel.grid.major.x = element_line(color = "grey90"),
          panel.grid.major.y = element_line(color = "grey90"))
  
  if(plotly == F){
    p1 <- p1 + geom_point(size = plot_point_size)
  }
  
  if(plotly == T){
    data_plotly <- paste("Kingdom: ", data$tax[,1],"<br>",
                         "Phylum: ", data$tax[,2],"<br>",
                         "Class: ", data$tax[,3],"<br>",
                         "Order: ", data$tax[,4],"<br>",
                         "Family: ", data$tax[,5],"<br>",
                         "Genus: ", data$tax[,6],"<br>",
                         "Species: ", data$tax[,7],"<br>",
                         "OTU: ", data$tax[,8],sep = "")
    
    p1 <- p1 + geom_point(size = plot_point_size-1, aes(text = data_plotly))
  }
  
  
  
  ### Points plot of significant differential abundant entries
  abund5 <- mutate(abund4, Tax = rownames(abund4)) %>%
    gather(key=Sample, value=Count, -Tax) %>%
    group_by(Sample) %>%
    mutate(Abundance = Count / sum(Count)*100)
  
  
  abund6 <- merge(abund5, res_tax, by = "Tax") %>%
    filter(padj < signif_thrh & fold < abs(log2FoldChange)) %>%
    arrange(padj)
  
  if(nrow(abund6) == 0){stop("No significant differences found.")}
  
  if(!is.null(adjust_zero)){
    abund6$Abundance[abund6$Abundance==0] <- adjust_zero
  }
  
  
  colnames(metadata)[1] <- "Sample"
  metadata <- metadata[c("Sample",group)]
  colnames(metadata)[2] <- "Group"
  
  point_df <- merge(x = abund6, y = metadata, by = "Sample") %>%
    group_by(Sample) %>%
    arrange(padj)
  
  colnames(point_df)[12] <- group
  
  if(!is.null(plot_nshow)){
    if(plot_nshow < nrow(abund6)){plot_nshow <- nrow(abund6)}
    point_df <- subset(point_df, Tax %in% as.character(unique(point_df$Tax))[1:plot_nshow])
  }
  
  point_df$Tax <- factor(point_df$Tax, levels = rev(as.character(unique(point_df$Tax))[1:plot_nshow]))
  
  p2 <- ggplot(data = point_df, aes_string(x = "Tax", y = "Abundance", color = group)) +
          labs(x = "", y = "Read Abundance (%)") +
          coord_flip() +
          theme_classic() +
          theme(panel.grid.major.x = element_line(color = "grey90"),
                panel.grid.major.y = element_line(color = "grey90"))
  
  if (plot_type == "point"){
    p2 <- p2 + geom_jitter(position = position_jitter(width = .05), size = plot_point_size)
  } else{
    p2 <- p2 + geom_boxplot(outlier.size=1)
  }
  
  
  
  clean_res0 <- merge(abund5, res_tax, by = "Tax") %>% 
                merge(y = metadata, by = "Sample") %>%
                group_by(Sample) %>%
                arrange(padj)
  
  colnames(clean_res0)[12] <- "group"

    cr <- mutate(clean_res0, 
                      padj = signif(padj, 2), 
                      Log2FC = signif(log2FoldChange, 2),
                      Taxonomy = Tax) %>%
    group_by(group, Taxonomy, padj, Log2FC) %>%
    summarise(Avg = round(mean(Abundance), 3)) %>%
    spread(key = group, value = Avg) %>%
    arrange(padj)
  
  out <- list(DESeq2_results = res, plot_MA = p1, DESeq2_results_significant = res_tax_sig, plot_sig = p2 , sig_res_plot_data = point_df, Clean_results = cr)
  
  if(plotly == T){ggplotly(p1, tooltip = "text")} 
  else{ return(out) }
}
