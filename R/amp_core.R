#' Generates a ggplot2 style core community plot
#'
#' Generates a ggplot2 style core community plot
#'
#' @usage amp_core(data)
#'
#' @param data (required) Data list as loaded with `amp_load()`.
#' @param group Group the data based on a sample variable (default: "Sample").
#' @param tax.empty Either "remove" OTUs without taxonomic information, add "best" classification or add the "OTU" name (default: best).
#' @param tax.class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (default: OTU).
#' @param abund.treshold Treshold for considering something abundant in percent (default: 0.1).
#' @param plotly Returns an interactive plot instead (default: F).
#' @param output Either plot or complete (default: "plot").
#' @param raw Display raw input instead of converting to percentages (default: F).
#' 
#' @return A ggplot2 object
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_core <- function(data, group = "Sample", abund.treshold = 0.1, tax.aggregate = "OTU", tax.class = NULL, tax.empty = "best", plotly = F, raw = F, output = "plot"){
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax.class = tax.class, tax.empty = tax.empty, tax.level = tax.aggregate)
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  sample <- data[["metadata"]]
  
  if (raw == F){
    abund <- as.data.frame(sapply(abund, function(x) x/sum(x)*100))
  }
  
  # Aggregate to a specific taxonomic level
  abund1 <- cbind.data.frame(Display = tax[,tax.aggregate], abund) %>%
    gather(key = Sample, value = Abundance, -Display)
  
  abund1 <- data.table(abund1)[, sum:=sum(Abundance), by=list(Display, Sample)] %>%
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
  

    abund3$HA <- ifelse(abund3$Abundance > abund.treshold, 1, 0)
    temp3 <- group_by(abund3, Display) %>%
      summarise(Frequency = sum(freq), freq_A= sum(HA), Abundance = round(mean(Abundance),2)) %>%
      as.data.frame()

    p <- ggplot(data = temp3, aes(x = Frequency, y = freq_A)) +
      ylab(paste("Abundant in N ", group, "s (>", abund.treshold, "%)" , sep="")) +
      xlab(paste("Observed in N ", group, "s", sep="")) +
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

  if(tax.aggregate == "OTU"){
    colnames(temp3)[1] <- "OTU"
    core <- merge(x = temp3, y = tax, by = "OTU") 
    temp3 <- core
  }
  
  if(plotly == T){
    ggplotly(p, tooltip = "text") %>% 
      layout(showlegend = FALSE)
  }
  else if(output == "complete"){ return(list(data = temp3, plot = p, abund = abund, tax = tax, sample = sample))}
  else if(output == "plot"){ return(p) }
}
