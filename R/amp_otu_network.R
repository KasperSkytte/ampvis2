#' Generate a taxa/sample network
#'
#' Generate a taxa/sample network
#'
#' @usage amp_otu_network(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param color A metadata variable to color the samples by.
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (default: "Phylum")
#' @param tax.add Additional taxonomic levels to display for each entry e.g. "Phylum" (default: "none") 
#' @param tax.show The number of taxa to show (default: "all").
#' @param tax.empty Either "remove" OTUs without taxonomic information, add "best" classification or add the "OTU" name (default: "best").
#' @param tax.class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param min.abundance Minimum taxa abundance pr. sample (default: 0)
#' @param output To output a plot or the complete data inclusive dataframes (default: "plot")
#' @param raw Display raw input instead of converting to percentages (default: F) 
#' 
#' @return A ggplot2 object or a list with the ggplot2 object and associated dataframes.
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_otu_network <- function(data, tax.aggregate = "Phylum", tax.add = NULL, tax.show = 10, tax.class = NULL, tax.empty = "best", output = "plot", raw = F, min.abundance = 0, color = NULL){
  
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
    gather(key = Sample, value = Abundance, -Display) %>%
    mutate(Display = paste("Taxa; ", Display))

    abund3 <- data.table(abund3)[, sum:=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    unique() %>% 
    as.data.frame()
  
  ## Add group information
    abund5 <- data.frame(abund3, Group = abund3$Sample)
  
  ## Take the average to group level
  
    abund6 <- data.table(abund5)[, Abundance:=mean(sum), by=list(Display, Group)] %>%
      setkey(Display, Group) %>%
      unique() %>% 
      as.data.frame() %>% 
      select(-sum)
  
  ## Find the X most abundant levels
    TotalCounts <- group_by(abund6, Display) %>%
      summarise(Abundance = sum(Abundance)) %>%
      arrange(desc(Abundance))

    if (tax.show > nrow(TotalCounts)){tax.show <- nrow(TotalCounts)}
    abund7 <- filter(abund6, Display %in% TotalCounts$Display[1:tax.show])
    
  ## Convert to network  
    netw <- data.frame(SeqID = as.character(abund7$Group), 
                       Taxa = abund7$Display, 
                       Abundance = abund7$Abundance, 
                       stringsAsFactors = F) %>%
      subset(Abundance > min.abundance) %>% # Subset to dominant species in each sample
      select(-Abundance) %>%
      network(directed = FALSE)
    
  ## Add data to nodes
    x = data.frame(SeqID = network.vertex.names(netw), stringsAsFactors = F)
    
    xsamples <- filter(x, !grepl("Taxa", SeqID)) %>%
      merge(metadata, all.x = T, by = "SeqID") 
    
    if (is.null(color)){
      xsamples$Description <- "Samples"
    } else{
      xsamples$Description <- as.character(xsamples[, color]) 
    }
    
    set.vertex.attribute(netw,"snames", c(rep("Sample", nrow(xsamples)), rep("Taxa", nrow(x)-nrow(xsamples))))
    set.vertex.attribute(netw,"stype", c(xsamples$Description, rep("Taxa", nrow(x)-nrow(xsamples))))
    set.vertex.attribute(netw,"nsize", c(rep(3, nrow(xsamples)), rep(1, nrow(x)-nrow(xsamples))))
    
  ## Make network plot
    
    p <- ggnet2(netw, label = F,  color = "stype", node.alpha = 0.7, edge.color = "grey80", node.size = "nsize") +
      scale_color_brewer(palette = "Set1", name = "") +
      scale_size_discrete(guide = F, range = c(3,6))

  ## Define the output 
  if (output == "complete"){
    outlist <- list(heatmap = p, data = abund7)
    return(outlist)  
  }
  if (output == "plot"){
    return(p)
  }
}
