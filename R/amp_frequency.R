#' Frequency plot
#'
#' Generates a barplot with frequency vs read abundance.
#'
#' @usage amp_frequency(data)
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
#' @param weight (\emph{logical}) Weight the frequency by abundance. (\emph{default:} \code{TRUE})
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
#' @examples 
#' #Load example data
#' data("AalborgWWTPs")
#' 
#' #Frequency plot
#' amp_frequency(AalborgWWTPs)
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_frequency <- function(data,
                          group_by = "Sample", 
                          tax_class = NULL,
                          tax_empty = "best",
                          tax_aggregate = "OTU",
                          weight = TRUE, 
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
  
  if (raw == FALSE){
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
  
  ## Make a nice frequency plot
  
    temp3 <- group_by(abund3, Display) %>%
      summarise(Frequency = sum(freq), Total = sum(sum)) %>%
      mutate(Percent = Total / length(unique(abund3$Group))) %>%
      as.data.frame()
    
    if(weight == T){
      p <- ggplot(data = temp3, aes(x = Frequency, weight = Percent)) +
        ylab("Read abundance (%)") +
        xlab(paste("Frequency (Observed in N ", group_by, "s)", sep = ""))
      if(max(temp3$Frequency) > 30){ p <- p + geom_bar()}
      if(max(temp3$Frequency) <= 30){ p <- p + geom_bar(binwidth = 1)}
    }
    
    if(weight == F){
      p <- ggplot(data = temp3, aes(x = Frequency)) +
        ylab("Count") +
        xlab(paste("Frequency (Observed in N ", group_by, "s)", sep=""))
      if(max(temp3$Frequency) > 30){ p <- p + geom_bar()}
      if(max(temp3$Frequency) <= 30){ p <- p + geom_bar(binwidth = 1)}
    } 

    p <- p +  
      theme_classic() +
      theme(panel.grid.major.x = element_line(color = "grey90"),
            panel.grid.major.y = element_line(color = "grey90"))
      
  if(tax_aggregate == "OTU"){
    colnames(temp3)[1] <- "OTU"
    core <- merge(x = temp3, y = tax, by = "OTU") 
    temp3 <- core
  }
  
    if (detailed_output) {
      return(list(data = temp3, plot = p, abund = abund, tax = tax, metadata = metadata))
    } else if (!detailed_output)
      return(p)
}
