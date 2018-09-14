#' Network plot
#'
#' Generates network plot of taxa and samples based on ggnet2.
#'
#' @usage amp_otu_network(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param color_by A metadata variable to color the samples by.
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"Phylum"})
#' @param tax_add Additional taxonomic level(s) to display, e.g. \code{"Phylum"}. (\emph{default:} \code{"none"})
#' @param tax_show The number of taxa to show, or a vector of taxa names. (\emph{default:} \code{10})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible. 
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param min_abundance Minimum taxa abundance pr. sample. (\emph{default:} \code{0})
#' @param normalise (\emph{logical}) Transform the OTU read counts to be in percent per sample. (\emph{default:} \code{TRUE})
#' @param detailed_output (\emph{logical}) Return additional details or not. If \code{TRUE}, it is recommended to save to an object and then access the additional data by \code{View(object$data)}. (\emph{default:} \code{FALSE})
#' 
#' @return A ggplot2 object. If \code{detailed_output = TRUE} a list with a ggplot2 object and additional data.
#' 
#' @import ggplot2
#' @importFrom dplyr filter arrange group_by mutate select summarise desc
#' @importFrom tidyr gather
#' @importFrom data.table as.data.table data.table setkey
#' @importFrom network network network.vertex.names set.vertex.attribute
#' @importFrom ggnet ggnet2
#' 
#' @export
#' 
#' @section Preserving relative abundances in a subset of larger data:
#' See \code{?\link{amp_subset_samples}} or the \href{https://madsalbertsen.github.io/ampvis2/articles/faq.html#preserving-relative-abundances-in-a-subset-of-larger-data}{ampvis2 FAQ}.
#' 
#' @details See \code{\link[GGally]{ggnet2}}
#' 
#' @seealso 
#' \code{\link{amp_load}}
#' 
#' @examples 
#' #Load example data
#' data("AalborgWWTPs")
#' 
#' #OTU network plot
#' amp_otu_network(AalborgWWTPs)
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_otu_network <- function(data,
                            min_abundance = 0,
                            color_by = NULL,
                            tax_aggregate = "Phylum",
                            tax_add = NULL,
                            tax_show = 10,
                            tax_class = NULL,
                            tax_empty = "best",
                            normalise = TRUE,
                            detailed_output = FALSE){
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis2 functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)", call. = FALSE)
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax_class = tax_class, tax_empty = tax_empty, tax_level = tax_aggregate)
  
  #tax_add and tax_aggregate can't be the same
  if(!is.null(tax_aggregate) & !is.null(tax_add)) {
    if(tax_aggregate == tax_add) {
      stop("tax_aggregate and tax_add cannot be the same", call. = FALSE)
    }
  }
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  metadata <- data[["metadata"]]
  
  ## SampleID column is used to merge data later, so it must be there!
  colnames(metadata)[1] <- "SampleID"
  
  if (isTRUE(normalise)){
    if(isTRUE(attributes(data)$normalised))
      warning("The data has already been normalised by either amp_subset_samples or amp_subset_taxa. Setting normalise = TRUE (the default) will normalise the data again and the relative abundance information about the original data of which the provided data is a subset will be lost.", call. = FALSE)
    #calculate sample percentages, skip columns with 0 sum to avoid NaN's
    abund[,which(colSums(abund) != 0)] <- as.data.frame(apply(abund[,which(colSums(abund) != 0), drop = FALSE], 2, function(x) x/sum(x)*100))
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
    tidyr::gather(key = Sample, value = Abundance, -Display) %>%
    mutate(Display = paste("Taxa; ", Display)) %>% as.data.table()

    abund3 <- abund3[, "sum":=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
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

    if (tax_show > nrow(TotalCounts)){tax_show <- nrow(TotalCounts)}
    abund7 <- filter(abund6, Display %in% TotalCounts$Display[1:tax_show])
    
  ## Convert to network  
    netw <- data.frame(SampleID = as.character(abund7$Group), 
                       Taxa = abund7$Display, 
                       Abundance = abund7$Abundance, 
                       stringsAsFactors = F) %>%
      subset(Abundance > min_abundance) %>% # Subset to dominant species in each sample
      select(-Abundance) %>%
      network::network(directed = FALSE)
    
  ## Add data to nodes
    x = data.frame(SampleID = network::network.vertex.names(netw), stringsAsFactors = F)
    
    xsamples <- filter(x, !grepl("Taxa", SampleID)) %>%
      merge(metadata, all.x = T, by = 1) 
    
    if (is.null(color_by)){
      xsamples$Description <- "Sample"
    } else{
      xsamples$Description <- as.character(xsamples[, color_by]) 
    }
    
    network::set.vertex.attribute(netw,"snames", c(rep("Sample", nrow(xsamples)), rep("Taxa", nrow(x)-nrow(xsamples))))
    network::set.vertex.attribute(netw,"stype", c(xsamples$Description, rep("Taxa", nrow(x)-nrow(xsamples))))
    network::set.vertex.attribute(netw,"nsize", c(rep(3, nrow(xsamples)), rep(1, nrow(x)-nrow(xsamples))))
    
  ## Make network plot
    
    p <- ggnet::ggnet2(netw, label = F,  color = "stype", node.alpha = 0.7, edge.color = "grey80", node.size = "nsize") +
      scale_color_brewer(palette = "Set1", name = "") +
      scale_size_discrete(guide = F, range = c(3,6))

  ## Define the output 
    if (detailed_output) {
      return(list(heatmap = p, data = abund7))
    } else 
      return(p)
}
