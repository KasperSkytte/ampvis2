#' Timeseries
#'
#' Generates a timeseries plot showing relative read abundances over time.
#' 
#' @usage amp_timeseries(data, time_variable = "")
#' 
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param time_variable (required) The name of the column in the metadata containing the time variables, e.g. \code{"Date"}.
#' @param group_by Group the samples by a variable in the metadata.
#' @param split Split the plot into subplots of each taxa. (\emph{default:} \code{FALSE}) 
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"OTU"})
#' @param tax_add Additional taxonomic level(s) to display, e.g. \code{"Phylum"}. (\emph{default:} \code{"none"})
#' @param tax_show The number of taxa to show, or a vector of taxa names. (\emph{default:} \code{5})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible. 
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param raw (\emph{logical}) Display raw input instead of converting to percentages. (\emph{default:} \code{FALSE}) 
#' 
#' @keywords timeseries
#' @import dplyr
#' @import ggplot2
#' @import data.table
#' @return A ggplot2 object.
#' 
#' @export 
#' 
#' @author Julie Klessner Thun Pedersen \email{julieklessnerthun@@gmail.com}

amp_timeseries <- function(data,
                           time_variable = NULL, 
                           group_by = "Sample", 
                           split = FALSE, 
                           tax_show = 5, 
                           tax_aggregate="OTU", 
                           tax_add=NULL, 
                           tax_class=NULL,
                           tax_empty="best", 
                           raw = FALSE){
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  # Required arguments
  if(is.null(time_variable))
    stop("Argument \"time_variable\" is required.")
  # Clean and rename taxonomy ---------------------------------------------------------
  
  data <- amp_rename(data = data, 
                     tax_class = tax_class, 
                     tax_empty = tax_empty, 
                     tax_level = tax_aggregate)
  
  # Divide data to seperate data frames ----------------------------------------------
  
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  metadata <- data[["metadata"]]
  
  # Convert to percentages -----------------------------------------------------------
  
  if (raw == F){
    abund <- as.data.frame(sapply(abund, function(x) x/sum(x)*100))
  }
  
  # Display multiple levels using tax_add argument ------------------------------------
  
  suppressWarnings(
    if (!is.null(tax_add)){
      
      if (tax_add != tax_aggregate){
        tax <- data.frame(tax, 
                          Display = apply(tax[,c(tax_add,tax_aggregate)], 1, 
                                          paste, collapse="; "))
      }
    } else {
      tax <- data.frame(tax, Display = tax[, tax_aggregate])
    }
  )  
  
  # Aggregate to a specific taxonomic level using tax_aggregate argument--------------
  
  abund3 <- cbind(Display = tax[,"Display"], abund) %>%
    melt(id.var = "Display", 
         value.name = "Abundance", 
         variable.name = "Sample")
  
  abund3 <- data.table(abund3)[, sum:=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    unique() %>% 
    as.data.frame()
  
  # Reshaping data  ---------------------------------------------------------
  
  abund4 <- abund3[, c("Display","Sample", "sum")]
  
  ## Add group information
  suppressWarnings(
    if (group_by != "Sample"){
      if (length(group_by) > 1){
        grp <- data.frame(Sample = metadata[,1], Group = apply(metadata[,group_by], 1, paste, collapse = " ")) 
        oldGroup <- unique(cbind.data.frame(metadata[,group_by], Group = grp$Group))
      } else{
        grp <- data.frame(Sample = metadata[,1], Group = metadata[,group_by]) 
      }
      abund5 <- merge(abund4, grp)
    } else{ abund5 <- data.frame(abund4, Group = "Sample")}
  )
  
  ## Find the x most abundant levels and sort
  TotalCounts <- group_by(abund5, Display) %>%
    summarise(Median = median(sum), Total = sum(sum), Mean = mean(sum)) %>% 
    arrange(desc(Mean))
  
  ## Subset to the x most abundant levels
  if (is.numeric(tax_show)){
    if (tax_show > nrow(TotalCounts)){  
      tax_show <- nrow(TotalCounts)
    }
    abund5$Display <- as.character(abund5$Display)
    abund7 <- subset(abund5, Display %in% as.character(unlist(TotalCounts[1:tax_show,"Display"])))
  }
  ## Subset to a list of level names
  if (!is.numeric(tax_show)){
    if (length(tax_show) > 1){
      abund7 <- subset(abund5, Display %in% tax_show)
    }
    if ((length(tax_show) == 1) && (tax_show != "all")){
      abund7 <- subset(abund5, Display %in% tax_show)
    }
    ### Or just show all  
    if ((length(tax_show) == 1) && (tax_show == "all")){
      tax_show <- nrow(TotalCounts)  
      abund7 <- subset(abund5, Display %in% as.character(unlist(TotalCounts[1:tax_show,"Display"])))  
    }
  }
  
  abund8 <- merge(abund7, metadata, by.x = "Sample", by.y = colnames(metadata)[1])
  
  abund8[, time_variable] <- as.Date(abund8[, time_variable])

  abund9 <- mutate(abund8, DG = paste(Display, Group))
    
  
  # Base plot, user may choose labels ------------------------------------------------
    
  
  if(length(levels(abund9$Group)) > 1) {
    p <- ggplot(abund9, aes_string(x=time_variable, y="sum", col = "Display", group = "DG", linetype = group_by))
  } else{
    p <- ggplot(abund9, aes_string(x=time_variable, y="sum", col = "Display", group = "DG"))
  }
  
  p <-  p +
         geom_line()+
         geom_point()+
         scale_color_discrete(name = "") +
         scale_linetype_discrete() +
         xlab("Date") +
         ylab("Read abundance (%)") +
         theme_classic() +
         theme(axis.text.x = element_text(size = 10, vjust = 0.3, angle = 90),
               panel.grid.major.x = element_line(color = "grey90"),
               panel.grid.major.y = element_line(color = "grey90"))
  
  if(split == T){
    p <- p + facet_wrap(~Display) +
      theme(strip.background = element_rect(colour=NA, fill="grey95"),
            panel.grid.major.x = element_line(color = "grey90"),
            panel.grid.major.y = element_line(color = "grey90"),
            legend.position = "bottom") +
      scale_color_discrete(guide = F)
  }
    
    
    # Output----------------------------------------------------------------------------
    return(p)
}
