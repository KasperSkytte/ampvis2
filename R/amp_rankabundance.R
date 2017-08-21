#' Rank abundance plot
#'
#' Generates a rank abundance curve (rank abundance vs cumulative read abundance) for each sample. 
#'
#' @usage amp_rankabundance(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param group_by Group the samples by a variable in the metadata.
#' @param order_group A vector to order the groups by.
#' @param tax_clean (\emph{logical}) Replace the phylum Proteobacteria with the respective Classes instead. (\emph{default:} \code{TRUE})
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"Genus"})
#' @param tax_add Additional taxonomic level(s) to display, e.g. \code{"Phylum"}. (\emph{default:} \code{"Phylum"})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible. 
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param plot_log (\emph{logical}) Log10-scale the plot. (\emph{default:} \code{FALSE})
#' @param raw (\emph{logical}) Display raw input instead of converting to percentages. (\emph{default:} \code{FALSE}) 
#' @param detailed_output (\emph{logical}) Return additional details or not. If \code{TRUE}, it is recommended to save to an object and then access the additional data by \code{View(object$data)}. (\emph{default:} \code{FALSE})
#' 
#' @return A ggplot2 object. If \code{detailed_output = TRUE} a list with a ggplot2 object and additional data.
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_rankabundance <- function(data,
                              group_by = "Sample",
                              tax_clean = TRUE,
                              plot_log = FALSE,
                              detailed_output = FALSE,
                              tax_add = NULL,
                              tax_aggregate = "Genus",
                              tax_empty = "best",
                              tax_class = NULL,
                              raw = FALSE,
                              order_group = NULL){
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax_class = tax_class, tax_empty = tax_empty, tax_level = tax_aggregate)
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  metadata <- data[["metadata"]]
  
  if (raw == FALSE){
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
    if (group_by != "Sample"){
      if (length(group_by) > 1){
        grp <- data.frame(Sample = rownames(metadata), Group = apply(metadata[,group_by], 1, paste, collapse = " ")) 
      } else{
        grp <- data.frame(Sample = rownames(metadata), Group = metadata[,group_by]) 
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
    
    if(!is.null(order_group)){
      TotalCounts$Group <- factor(TotalCounts$Group, levels = rev(order_group))
    }
    
    p <- ggplot(data = TotalCounts, aes(x = Rank, y = Cumsum, color = Group)) +
      geom_line(size = 1) +
      ylim(0,100) +
      xlab("Rank abundance") +
      ylab("Cumulative read abundance (%)") +
      theme_classic()
    
    if (plot_log ==T){
      p <- p + scale_x_log10() 
    } 
    
    ## Define the output 
    if (detailed_output){
      outlist <- list(plot = p, data = TotalCounts)
      return(outlist)  
    }
    if (!detailed_output)
      return(p)
}
