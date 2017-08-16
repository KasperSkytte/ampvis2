#' Generates a ggplot2 style frequency plot
#'
#' Generates a ggplot2 style frequency plot
#'
#' @usage amp_frequency(data)
#'
#' @param data (required) Data list as loaded with `amp_load()`.
#' @param group Group the data based on a sample variable (default: "Sample").
#' @param tax.empty Either "remove" OTUs without taxonomic information, add "best" classification or add the "OTU" name (default: best).
#' @param tax.class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (default: OTU).
#' @param weight Weight the frequency by abundance (default: T).
#' @param output Either plot or complete (default: "plot").
#' @param raw Display raw input instead of converting to percentages (default: F).
#' 
#' @return A ggplot2 object
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_frequency <- function(data, group = "Sample", tax.class = NULL, tax.empty = "best", output = "plot",  tax.aggregate = "OTU", weight = T, raw = F){
  
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
  
  ## Make a nice frequency plot
  
    temp3 <- group_by(abund3, Display) %>%
      summarise(Frequency = sum(freq), Total = sum(sum)) %>%
      mutate(Percent = Total / length(unique(abund3$Group))) %>%
      as.data.frame()
    
    if(weight == T){
      p <- ggplot(data = temp3, aes(x = Frequency, weight = Percent)) +
        ylab("Read abundance (%)") +
        xlab(paste("Frequency (Observed in N ", group, "s)", sep = ""))
      if(max(temp3$Frequency) > 30){ p <- p + geom_bar()}
      if(max(temp3$Frequency) <= 30){ p <- p + geom_bar(binwidth = 1)}
    }
    
    if(weight == F){
      p <- ggplot(data = temp3, aes(x = Frequency)) +
        ylab("Count") +
        xlab(paste("Frequency (Observed in N ", group, "s)", sep=""))
      if(max(temp3$Frequency) > 30){ p <- p + geom_bar()}
      if(max(temp3$Frequency) <= 30){ p <- p + geom_bar(binwidth = 1)}
    } 

    p <- p +  
      theme_classic() +
      theme(panel.grid.major.x = element_line(color = "grey90"),
            panel.grid.major.y = element_line(color = "grey90"))
      
  if(tax.aggregate == "OTU"){
    colnames(temp3)[1] <- "OTU"
    core <- merge(x = temp3, y = tax, by = "OTU") 
    temp3 <- core
  }
  
  if(output == "complete"){ return(list(data = temp3, plot = p, abund = abund, tax = tax, sample = sample)) }
  if(output == "plot"){ return(p) }
}
