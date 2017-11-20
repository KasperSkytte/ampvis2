#' Core community analysis
#'
#' Generates a core community plot.
#'
#' @usage amp_core(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param group_by Group the samples by a variable in the metadata. (\emph{default:} \code{"Sample"})
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"OTU"})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible. 
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param abund_thrh Threshold in percent for defining "abundant"/"core" taxa. (\emph{default:} \code{0.1})
#' @param plotly (\emph{logical}) Returns an interactive plot instead. (\emph{default:} \code{FALSE})
#' @param raw (\emph{logical}) Display raw input instead of converting to percentages. (\emph{default:} \code{FALSE})
#' @param detailed_output (\emph{logical}) Return additional details or not. If \code{TRUE}, it is recommended to save to an object and then access the additional data by \code{View(object$data)}. (\emph{default:} \code{FALSE})
#' 
#' @return A ggplot2 object. If \code{detailed_output = TRUE} a list with a ggplot2 object and additional data.
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @import data.table
#' @export
#' 
#' @details Plots the number of samples a given OTU is observed in (first axis) against the number of samples it is abundant in (second axis), as defined by the threshold. By setting the \code{group_by} argument to a variable in the metadata, then the axes will show their frequency in this group instead of per sample. 
#' 
#' The OTU points can be aggregated by the taxonomy by setting fx \code{tax_aggregate = "Phylum"}. To see the corresponding taxonomy of the OTUs use \code{plotly = TRUE} for an interactive plot.
#' 
#' @examples 
#' #Load example data
#' data("AalborgWWTPs")
#' 
#' #Core plot
#' amp_core(AalborgWWTPs)
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_core <- function(data, 
                     group_by = "Sample", 
                     abund_thrh = 0.1, 
                     tax_aggregate = "OTU", 
                     tax_class = NULL,
                     tax_empty = "best", 
                     plotly = FALSE, 
                     raw = FALSE,
                     detailed_output = FALSE){
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax_class = tax_class, tax_empty = tax_empty, tax_level = tax_aggregate)
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  metadata <- data[["metadata"]]
  
  if (raw == F){
    abund <- as.data.frame(sapply(abund, function(x) x/sum(x)*100))
  }
  
  # Aggregate to a specific taxonomic level
  abund1 <- cbind.data.frame(Display = tax[,tax_aggregate], abund) %>%
    tidyr::gather(key = Sample, value = Abundance, -Display) %>% as.data.table()
  
  abund1 <- abund1[, "sum":=sum(Abundance), by=list(Display, Sample)] %>%
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
      abund1$Group <- grp$Group[match(abund1$Sample, grp$Sample)]
      abund2 <- abund1
    } else{ abund2 <- data.frame(abund1, Group = abund1$Sample)}
  )
  
  ## Take the average to group level
  abund3 <- data.table(abund2)[, Abundance:=mean(sum), by=list(Display, Group)] %>%
    setkey(Display, Group) %>%
    unique() %>%
    filter(Abundance > 0) %>%  
    mutate(freq = 1)
  

    abund3$HA <- ifelse(abund3$Abundance > abund_thrh, 1, 0)
    temp3 <- group_by(abund3, Display) %>%
      summarise(Frequency = sum(freq), freq_A= sum(HA), Abundance = round(mean(Abundance),2)) %>%
      as.data.frame()

    p <- ggplot(data = temp3, aes(x = Frequency, y = freq_A)) +
      ylab(paste("Abundant in N ", group_by, "s (>", abund_thrh, "%)" , sep="")) +
      xlab(paste("Observed in N ", group_by, "s", sep="")) +
      theme_classic() +
      theme(panel.grid.major.x = element_line(color = "grey90"),
            panel.grid.major.y = element_line(color = "grey90"))
    
    if(plotly == F){
      p <- p + geom_jitter(size = 3, alpha = 0.5)
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
      
      p <- p + geom_jitter(size = 2, alpha = 0.5, aes(text = data_plotly))
      }

  if(tax_aggregate == "OTU"){
    colnames(temp3)[1] <- "OTU"
    core <- merge(x = temp3, y = tax, by = "OTU") 
    temp3 <- core
  }
  
  if(plotly){
    ggplotly(p, tooltip = "text") %>% 
      layout(showlegend = FALSE)
  } else if (!plotly) {
    if (detailed_output) {
      return(list(data = temp3, plot = p, abund = abund, tax = tax, metadata = metadata))
    } else if (!detailed_output)
      return(p)
  }
}
