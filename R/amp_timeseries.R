#' Timeseries
#'
#' This function allows you plot relative abundance reads over time in a timeseries plot.
#' 
#' @usage amp_timeseries(data)
#' 
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param time (required) The name of the column containing the time variables, e.g. "Date".
#' @param group A variable from the associated sample data to group/split samples by (default: "Samples").
#' @param split Split the plot into subplot of each taxa (default: F).
#' @param tax.show The number of taxa to show or a vector of taxa names (default: 5).
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (default: "OTU")
#' @param tax.add Additional taxonomic levels to display for each entry, e.g. "Phylum" (default: "none").
#' @param tax.class Converts a specific phyla to class level instead, e.g. "p__Proteobacteria" (default: "none) .
#' @param tax.empty Either "remove" OTUs without taxonomic information, add "best" classification or add the "OTU" name (default: "best").
#' @param raw Display raw input instead of converting to percentages (default: F)  
#' 
#' @keywords timeseries
#' 
#' @return A ggplot2 object or a list with the ggplot2 object and associated dataframes.
#' 
#' @export 
#' 
#' @author Julie Klessner Thun Pedersen \email{julieklessnerthun@@gmail.com}

amp_timeseries <- function(data, time, group = "Sample", split = F, tax.show = 5, tax.aggregate="OTU", tax.add=NULL, tax.class=NULL, tax.empty="best", layout="dygraph", raw = F){
  
  # Clean and rename taxonomy ---------------------------------------------------------
  
  data <- amp_rename(data = data, 
                     tax.class = tax.class, 
                     tax.empty = tax.empty, 
                     tax.level = tax.aggregate)
  
  # Divide data to seperate data frames ----------------------------------------------
  
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  sample <- data[["metadata"]]
  
  # Convert to percentages -----------------------------------------------------------
  
  if (raw == F){
    abund <- as.data.frame(sapply(abund, function(x) x/sum(x)*100))
  }
  
  # Display multiple levels using tax.add argument ------------------------------------
  
  suppressWarnings(
    if (!is.null(tax.add)){
      
      if (tax.add != tax.aggregate){
        tax <- data.frame(tax, 
                          Display = apply(tax[,c(tax.add,tax.aggregate)], 1, 
                                          paste, collapse="; "))
      }
    } else {
      tax <- data.frame(tax, Display = tax[, tax.aggregate])
    }
  )  
  
  # Aggregate to a specific taxonomic level using tax.aggregate argument--------------
  
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
    if (group != "Sample"){
      if (length(group) > 1){
        grp <- data.frame(Sample = sample$SeqID, Group = apply(sample[,group], 1, paste, collapse = " ")) 
        oldGroup <- unique(cbind.data.frame(sample[,group], Group = grp$Group))
      } else{
        grp <- data.frame(Sample = sample$SeqID, Group = sample[,group]) 
      }
      abund5 <- merge(abund4, grp)
    } else{ abund5 <- data.frame(abund4, Group = "Sample")}
  )
  
  ## Find the x most abundant levels and sort
  TotalCounts <- group_by(abund5, Display) %>%
    summarise(Median = median(sum), Total = sum(sum), Mean = mean(sum)) %>% 
    arrange(desc(Mean))
  
  ## Subset to the x most abundant levels
  if (is.numeric(tax.show)){
    if (tax.show > nrow(TotalCounts)){  
      tax.show <- nrow(TotalCounts)
    }
    abund5$Display <- as.character(abund5$Display)
    abund7 <- subset(abund5, Display %in% as.character(unlist(TotalCounts[1:tax.show,"Display"])))
  }
  ## Subset to a list of level names
  if (!is.numeric(tax.show)){
    if (length(tax.show) > 1){
      abund7 <- subset(abund5, Display %in% tax.show)
    }
    if ((length(tax.show) == 1) && (tax.show != "all")){
      abund7 <- subset(abund5, Display %in% tax.show)
    }
    ### Or just show all  
    if ((length(tax.show) == 1) && (tax.show == "all")){
      tax.show <- nrow(TotalCounts)  
      abund7 <- subset(abund5, Display %in% as.character(unlist(TotalCounts[1:tax.show,"Display"])))  
    }
  }
  
  abund8 <- merge(abund7, sample, by.x = "Sample", by.y = colnames(sample)[1])
  
  abund8[, time] <- as.Date(abund8[, time])

  abund9 <- mutate(abund8, DG = paste(Display, Group))
    
  
  # Base plot, user may choose labels ------------------------------------------------
    
  
  if(length(levels(abund9$Group)) > 1) {
    p <- ggplot(abund9, aes_string(x=time, y="sum", col = "Display", group = "DG", linetype = group))
  } else{
    p <- ggplot(abund9, aes_string(x=time, y="sum", col = "Display", group = "DG"))
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
