#' Generates a ggplot2 style rank abundance plots
#'
#' A nice long description
#'
#' @usage amp_rankabundance(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param group Group the data based on a sample variable.
#' @param order.group A vector defining the order of groups.
#' @param tax.clean Replace the phylum Proteobacteria with the respective Classes instead (default: T).
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (defualt: "Genus")
#' @param tax.add Additional taxonomic levels to display for each entry (default: "Phylum") 
#' @param tax.empty Either "remove" OTUs without taxonomic information, "rename" with best classification or add the "OTU" name (default: rename).
#' @param tax.class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param plot.log Log10 scale the data (default: F)
#' @param output Either "plot" or "complete" (default: "plot").
#' @param raw Display raw input instead of converting to percentages (default: F).
#' 
#' @return A ggplot2 object
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_rankabundance <- function(data, group = "Sample", tax.clean = T, plot.log = F, output = "plot", tax.add = NULL, tax.aggregate = "Genus", tax.empty = "best", tax.class = NULL, raw = F, order.group = NULL){
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax.class = tax.class, tax.empty = tax.empty, tax.level = tax.aggregate)
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  sample <- data[["metadata"]]
  
  if (raw == F){
    abund <- as.data.frame(sapply(abund, function(x) x/sum(x)*100))
  }
  
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
  
  
    temp3 <- group_by(abund5, Display, Group) %>%
      summarise(Mean = mean(Abundance))
    
    TotalCounts <- temp3[with(temp3, order(-Mean)),] %>%
      group_by(Group) %>%
      mutate(dummy = 1) %>%
      mutate(Cumsum = cumsum(Mean), Rank = cumsum(dummy)) %>%
      as.data.frame()
    
    if(!is.null(order.group)){
      TotalCounts$Group <- factor(TotalCounts$Group, levels = rev(order.group))
    }
    
    p <- ggplot(data = TotalCounts, aes(x = Rank, y = Cumsum, color = Group)) +
      geom_line(size = 1) +
      ylim(0,100) +
      xlab("Rank abundance") +
      ylab("Cumulative read abundance (%)") +
      theme_classic()
    
    if (plot.log ==T){
      p <- p + scale_x_log10() 
    } 
    
    outlist <- list(plot = p, data = TotalCounts)
  
  if(output == "complete"){ return(outlist) }
  if(output == "plot"){ return(p) }
}
