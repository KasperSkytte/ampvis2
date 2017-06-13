#' Calculate rarefaction curve for each sample.
#'
#' Calculate rarefaction curve for each sample using the vegan rarecurve function directly from an ampvis object.
#'
#' @usage amp_rarecurve(data)
#'
#' @param data (required) A ampvis object.
#' @param step Step size for sample sizes in rarefaction curves (default: 1000).
#' @param color A metadata variable to color by.
#' 
#' @export
#' @import vegan
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}


amp_rarecurve <- function (data, step = 1000, color = NULL){

  abund <- data[["abund"]] %>% as.matrix() %>% t()
  sample <- data[["metadata"]]

  if (!identical(all.equal(abund, round(abund)), TRUE)) stop("function accepts only integers (counts)")

  tot <- rowSums(abund)
  nr <- nrow(abund)
  out <- lapply(seq_len(nr), 
                function(i) {
                             n <- seq(1, tot[i], by = step)
                             if (n[length(n)] != tot[i]) 
                             n <- c(n, tot[i])
                             drop(rarefy(abund[i, ], n))
                             }
                )
  
  df <- data.frame(Reads = as.numeric(), Species = as.numeric(), SeqID = as.character())
  
  for (i in 1:length(out)){
    tsample <- attributes(out[[i]])$Subsample[length(out[[i]])] %>% names()
    tspecies <- unlist(out[[i]])
    treads <- attributes(out[[i]])$Subsample
    tdf <- data.frame(Reads = treads, Species = tspecies, SeqID = tsample)
    df <- rbind.data.frame(df, tdf)
  }
  
  dfm <- merge(sample, df, by = "SeqID") # Could add a check if this coloumn is correct
  
  ## Plot the data
  p <- ggplot(dfm, aes_string(x = "Reads", y = "Species", group = "SeqID", color = color)) +
    geom_line() +
    theme_classic() +
    xlab("Sequencing depth (reads)") +
    ylab("Number of observed OTUs")
  
  return(p)
}