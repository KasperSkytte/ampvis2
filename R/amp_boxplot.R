#' Generates a ggplot2 style boxplot of the most abundant taxa
#'
#' Generates a ggplot2 style boxplot of the most abundant taxa
#'
#' @usage amp_boxplot(data)
#'
#' @param data (required) Data list as loaded with `amp_load()`.
#' @param group Group the data based on a sample variable.
#' @param order_group A vector defining the order of groups.
#' @param order_y A vector to order the y-axis by.
#' @param tax_show The number of taxa to show or a vector of taxa names (default: 50).
#' @param tax_clean Replace the phylum Proteobacteria with the respective Classes instead (default: T).
#' @param tax_aggregate The taxonomic level that the data should be aggregated to (defualt: "Genus")
#' @param tax_add Additional taxonomic levels to display for each entry (default: "Phylum") 
#' @param tax_empty Either "remove" OTUs without taxonomic information, "rename" with best classification or add the "OTU" name (default: rename).
#' @param tax_class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param plot_flip Flip the axis of the plot (default: F).
#' @param plot_log Log10 scale the data (default: F)
#' @param adjust_zero Keep 0 abundances in ggplot2 median calculations by adding a small constant to these.
#' @param point_size Size of points (default: 1).
#' @param output Either "plot" or "complete" (default: "plot").
#' @param sort_by Sort the boxplot by either "median", "mean" or "total" (default = "median").
#' @param plot_type Show data using boxplot or points (default: boxplot).
#' @param raw Display raw input instead of converting to percentages (default: F).
#' 
#' @return A ggplot2 object
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_boxplot <- function(data, group = "Sample", order_group = NULL, tax_show = 50, tax_clean = T, plot_log = F, output = "plot", tax_add = NULL, tax_aggregate = "Genus", tax_empty = "best", tax_class = NULL, point_size = 1, plot_flip = F, sort_by = "median", adjust_zero = NULL, order_y = NULL, raw = F, plot_type = "boxplot"){
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax_class = tax_class, tax_empty = tax_empty, tax_level = tax_aggregate)
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  sample <- data[["metadata"]]
  
  if (raw == F){
    abund <- as.data.frame(sapply(abund, function(x) x/sum(x)*100))
  }
  
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
    gather(key = Sample, value = Abundance, -Display)
  
  abund3 <- data.table(abund3)[, Abundance:=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    unique() %>% 
    as.data.frame()
  
  ## Add group information
  suppressWarnings(
    if (group != "Sample"){
      if (length(group) > 1){
        grp <- data.frame(Sample = rownames(sample), Group = apply(sample[,group], 1, paste, collapse = " ")) 
      } else{
        grp <- data.frame(Sample = rownames(sample), Group = sample[,group]) 
      }
      abund3$Group <- grp$Group[match(abund3$Sample, grp$Sample)]
      abund5 <- abund3
    } else{ abund5 <- data.frame(abund3, Group = abund3$Sample)}
  ) 
  
    ## Find the x most abundant levels and sort
    TotalCounts <- group_by(abund5, Display) %>%
      summarise(Median = median(Abundance), Total = sum(Abundance), Mean = mean(Abundance))
    if(sort_by == "median"){TotalCounts %<>% arrange(desc(Median)) %>% as.data.frame()}
    if(sort_by == "mean"){TotalCounts %<>% arrange(desc(Mean)) %>% as.data.frame()}
    if(sort_by == "total"){TotalCounts %<>% arrange(desc(Total)) %>% as.data.frame()}
    
    abund5$Display <- factor(abund5$Display, levels = rev(TotalCounts$Display))
    
    ## Subset to the x most abundant levels
    if (is.numeric(tax_show)){
      if (tax_show > nrow(TotalCounts)){  
        tax_show <- nrow(TotalCounts)
      }
      abund7 <- subset(abund5, abund5$Display %in% TotalCounts[1:tax_show,"Display"])  
    }
    ## Subset to a list of level names
    if (!is.numeric(tax_show)){
      if (length(tax_show) > 1){
        abund7 <- subset(abund5, as.character(abund5$Display) %in% tax_show)
      }
      if ((length(tax_show) == 1) && (tax_show != "all")){
        abund7 <- subset(abund5, as.character(abund5$Display) %in% tax_show)
      }
      ### Or just show all  
      if ((length(tax_show) == 1) && (tax_show == "all")){
        tax_show <- nrow(TotalCounts)  
        abund7 <- subset(abund5, abund5$Display %in% TotalCounts[1:tax_show,"Display"])  
      }
    }
    
    ## Add a small constant to handle ggplot2 removal of 0 values in log scaled plots
    if(!is.null(adjust_zero)){
      abund7$Abundance[abund7$Abundance==0] <- adjust_zero
    }
    
    ## Order y based on a vector
    if (length(order_y) > 1){
      abund7$Display <- factor(abund7$Display, levels = order_y)
      abund7 <- subset(abund7, !is.na(Display))
    }
    
    ## plot the data
    if (group == "Sample"){
      p <-ggplot(abund7, aes(x = Display, y = Abundance))   
    }
    if (group != "Sample"){
      if(!is.null(order_group)){
        abund7$Group <- factor(abund7$Group, levels = rev(order_group))
      }
      p <-ggplot(abund7, aes(x = Display, y = Abundance, color = Group))   
    }
    
    p <- p +  ylab("Read Abundance (%)") + 
      guides(col = guide_legend(reverse = TRUE)) + 
      xlab("") +
      theme_classic() +
      theme(panel.grid.major.x = element_line(color = "grey90"),
            panel.grid.major.y = element_line(color = "grey90"))
    
    if (plot_flip == F){ p <- p + coord_flip() } else{
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }
    
    if (plot_type == "point"){ p <- p + geom_point(size = point_size) }
    if (plot_type == "boxplot"){p <- p + geom_boxplot(outlier.size = point_size)}
    if (plot_log ==T){ p <- p + scale_y_log10()}
    
    outlist <- list(plot = p, data = abund7)
  
  if(output == "complete"){ return(outlist) }
  if(output == "plot"){ return(p) }
}
