#' Boxplot
#'
#' Generates boxplots of the most abundant taxa.
#'
#' @usage amp_boxplot(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param group_by Group the samples by a variable in the metadata.
#' @param order_group A vector to order the groups by.
#' @param order_y A vector to order the y-axis by.
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"Genus"})
#' @param tax_add Additional taxonomic level(s) to display, e.g. \code{"Phylum"}. (\emph{default:} \code{"none"})
#' @param tax_show The number of taxa to show, or a vector of taxa names. (\emph{default:} \code{20})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible. 
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param plot_flip (\emph{logical}) Flip the axes of the plot axis. (\emph{default:} \code{FALSE})
#' @param plot_log (\emph{logical}) Log10-scale the plot. (\emph{default:} \code{FALSE})
#' @param adjust_zero Keep abundances of 0 in the calculation of medians by adding this value. (\emph{default:} \code{NULL})
#' @param point_size The size of points. (\emph{default:} \code{1})
#' @param sort_by Sort the boxplots by \code{"median"}, \code{"mean"} or \code{"total"}. (\emph{default:} \code{"median"})
#' @param plot_type Plot type. \code{"boxplot"} or \code{"point"}. (\emph{default:} \code{"boxplot"})
#' @param raw (\emph{logical}) Display raw input instead of converting to percentages. (\emph{default:} \code{FALSE})
#' @param detailed_output (\emph{logical}) Return additional details or not. If \code{TRUE}, it is recommended to save to an object and then access the additional data by \code{View(object$data)}. (\emph{default:} \code{FALSE})
#' 
#' @return A ggplot2 object. If \code{detailed_output = TRUE} a list with a ggplot2 object and additional data.
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import data.table
#' @export
#' 
#' @examples 
#' #Load example data
#' data("AalborgWWTPs")
#' 
#' #10 boxplots grouped by WWTP with phylum name added
#' amp_boxplot(AalborgWWTPs, 
#'             group_by = "Plant",
#'             tax_show = 10,
#'             tax_add = "Phylum")
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_boxplot <- function(data, 
                        group_by = "Sample",
                        sort_by = "median",
                        plot_type = "boxplot",
                        point_size = 1,
                        tax_aggregate = "Genus",
                        tax_add = NULL, 
                        tax_show = 20, 
                        tax_empty = "best", 
                        tax_class = NULL,
                        order_group = NULL,
                        order_y = NULL,
                        plot_flip = FALSE, 
                        plot_log = FALSE, 
                        adjust_zero = NULL,
                        raw = FALSE,
                        detailed_output = FALSE){
  
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
    tidyr::gather(key = Sample, value = Abundance, -Display) %>% as.data.table()
  
  abund3 <- abund3[,"Abundance" := sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    unique() 
  
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
    if (group_by == "Sample"){
      p <-ggplot(abund7, aes(x = Display, y = Abundance))   
    }
    if (group_by != "Sample"){
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
    if (plot_log == TRUE){ p <- p + scale_y_log10()}
  
    if (detailed_output) {
      return(list(plot = p, data = abund7))
    } else 
      return(p)
}
