#' Tests if there is a significant difference in abundance between selected conditions
#'
#' Tests if there is a significant difference in abundance between selected conditions
#'
#' @usage amp_test_species(data, design)
#'
#' @param data (required) Data list as loaded with `amp_load()`.
#' @param group (required) The group to test against.
#' @param sig Significance treshold (default: 0.01).
#' @param fold Log2fold filter default for displaying significant results (default: 0)
#' @param tax.aggregate Group data at specific taxonomic level (default: "OTU").
#' @param tax.class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param tax.empty Either "remove" OTUs without taxonomic information at X level, with "best" classification or add the "OTU" name (default: best).
#' @param tax.display Display additional taxonomic levels in the plot output e.g. "Genus".
#' @param label Label the significant entries with tax.display (default:F).
#' @param plot.type Either "boxplot" or "point" (default: point)
#' @param plot.show Display the X most significant results (default: 10).
#' @param plot.point.size The size of the plotted points.
#' @param adjust.zero Keep 0 abundances in ggplot2 median calculations by adding a small constant to these.
#' @param plotly Returns an interactive plot instead (default: F).
#' 
#' @return A p-value for each comparison.
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_test_species <- function(data, group, tax.aggregate = "OTU", tax.add = NULL, test = "Wald", fitType = "parametric", sig = 0.01, fold = 0, tax.class = NULL, tax.empty = "best", label = F, plot.type = "point", plot.show = 10, plot.point.size = 2, adjust.zero = NULL, plotly = F){
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax.class = tax.class, tax.empty = tax.empty, tax.level = tax.aggregate)
  
  ## Extract the data into seperate objects for readability
  abund <- data[["abund"]]  
  tax <- data[["tax"]]
  sample <- data[["metadata"]]

  ## Make a name variable that can be used instead of tax.aggregate to display multiple levels 
  suppressWarnings(
    if (!is.null(tax.add)){
      if (tax.add != tax.aggregate) {
        tax <- data.frame(tax, Display = apply(tax[,c(tax.add,tax.aggregate)], 1, paste, collapse="; "))
      }
    } else {
      tax <- data.frame(tax, Display = tax[,tax.aggregate])
    }
  )  
  
  # Aggregate to a specific taxonomic level
  abund3 <- cbind.data.frame(Display = tax[,"Display"], abund) %>%
    gather(key = Sample, value = Abundance, -Display)
  
  abund3 <- data.table(abund3)[, sum:=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    unique() %>% 
    as.data.frame() %>%
    select(-Abundance)
  
  ## Convert to DESeq2 format
  abund4 <- spread(data = abund3, key = Sample, value = sum)
  rownames(abund4) <- abund4$Display
  abund4 <- abund4[,-1]
  
  groupF <- as.formula(paste("~", group, sep=""))
  
  
  data_deseq <- DESeqDataSetFromMatrix(countData = abund4,
                                       colData = sample,
                                       design = groupF)
  
  ## Test for significant differential abundance
  data_deseq_test = DESeq(data_deseq, test=test, fitType=fitType)
  
  ## Extract the results
  res = results(data_deseq_test, cooksCutoff = FALSE)  
  res_tax = cbind(as.data.frame(res), Tax = rownames(res))
  
  res_tax_sig = subset(res_tax, padj < sig & fold < abs(log2FoldChange)) %>%
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
    p1 <- p1 + geom_point(size = plot.point.size)
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
    
    p1 <- p1 + geom_point(size = plot.point.size-1, aes(text = data_plotly))
  }
  
  
  
  ### Points plot of significant differential abundant entries
  abund5 <- mutate(abund4, Tax = rownames(abund4)) %>%
    gather(key=Sample, value=Count, -Tax) %>%
    group_by(Sample) %>%
    mutate(Abundance = Count / sum(Count)*100)
  
  
  abund6 <- merge(abund5, res_tax, by = "Tax") %>%
    filter(padj < sig & fold < abs(log2FoldChange)) %>%
    arrange(padj)
  
  if(nrow(abund6) == 0){stop("No significant differences found.")}
  
  if(!is.null(adjust.zero)){
    abund6$Abundance[abund6$Abundance==0] <- adjust.zero
  }
  
  
  colnames(sample)[1] <- "Sample"
  sample <- sample[c("Sample",group)]
  colnames(sample)[2] <- "Group"
  
  point_df <- merge(x = abund6, y = sample, by = "Sample") %>%
    group_by(Sample) %>%
    arrange(padj)
  
  colnames(point_df)[12] <- group
  
  if(!is.null(plot.show)){
    if(plot.show < nrow(abund6)){plot.show <- nrow(abund6)}
    point_df <- subset(point_df, Tax %in% as.character(unique(point_df$Tax))[1:plot.show])
  }
  
  point_df$Tax <- factor(point_df$Tax, levels = rev(as.character(unique(point_df$Tax))[1:plot.show]))
  
  p2 <- ggplot(data = point_df, aes_string(x = "Tax", y = "Abundance", color = group)) +
          labs(x = "", y = "Read Abundance (%)") +
          coord_flip() +
          theme_classic() +
          theme(panel.grid.major.x = element_line(color = "grey90"),
                panel.grid.major.y = element_line(color = "grey90"))
  
  if (plot.type == "point"){
    p2 <- p2 + geom_jitter(position = position_jitter(width = .05), size = plot.point.size)
  } else{
    p2 <- p2 + geom_boxplot(outlier.size=1)
  }
  
  
  
  clean_res0 <- merge(abund5, res_tax, by = "Tax") %>% 
                merge(y = sample, by = "Sample") %>%
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
