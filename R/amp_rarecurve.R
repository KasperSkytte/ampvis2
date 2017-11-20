#' Rarefaction curve
#'
#' Generates a rarefaction curve (number of reads vs number of observed OTUs) for each sample. 
#'
#' @usage amp_rarecurve(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param stepsize Step size for the curves. Lower is prettier but takes more time to generate. (\emph{default:} \code{1000})
#' @param color_by Color curves by a variable in the metadata. 
#' 
#' @export
#' @import dplyr
#' @import vegan
#' @import ggplot2
#' 
#' @return A ggplot2 object.
#' @examples 
#' #Load example data
#' data("AalborgWWTPs")
#' 
#' #Rarecurve
#' amp_rarecurve(AalborgWWTPs)
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}


amp_rarecurve <- function (data,
                           stepsize = 1000,
                           color_by = NULL){
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  maxreads <- max(colSums(data$abund))
  if(maxreads < stepsize) {
    stop("\"stepsize\" too high, maximum number of reads in any sample is: ", maxreads)
  }
  
  abund <- data[["abund"]] %>% as.matrix() %>% t()
  metadata <- data[["metadata"]]

  if (!identical(all.equal(abund, round(abund)), TRUE)) stop("Function accepts only integers (counts)")

  tot <- rowSums(abund)
  nr <- nrow(abund)
  out <- lapply(seq_len(nr), 
                function(i) {
                             n <- seq(1, tot[i], by = stepsize)
                             if (n[length(n)] != tot[i]) 
                             n <- c(n, tot[i])
                             drop(vegan::rarefy(abund[i, ], n))
                             }
                )
  
  df <- data.frame(Reads = as.numeric(), Species = as.numeric(), SampleID = as.character())
  
  for (i in 1:length(out)){
    tsample <- attributes(out[[i]])$Subsample[length(out[[i]])] %>% names()
    tspecies <- unlist(out[[i]])
    treads <- attributes(out[[i]])$Subsample
    tdf <- data.frame(Reads = treads, Species = tspecies, SampleID = tsample)
    df <- rbind.data.frame(df, tdf)
  }
  
  metadata_col1name <- colnames(metadata)[1]
  colnames(df)[which(colnames(df) == "SampleID")] <- metadata_col1name
  dfm <- merge(metadata, df, by = metadata_col1name)
  
  ## Plot the data
  p <- ggplot(dfm, aes_string(x = "Reads", y = "Species", group = metadata_col1name, color = color_by)) +
    geom_line() +
    theme_classic() +
    xlab("Sequencing depth (reads)") +
    ylab("Number of observed OTUs")
  
  return(p)
}