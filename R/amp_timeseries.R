#' Abundance timeseries
#'
#' Generates a timeseries plot showing relative read abundances over time.
#' 
#' @usage amp_timeseries(data, time_variable = "")
#' 
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param time_variable (required) The name of the column in the metadata containing the time variables, e.g. \code{"Date"}. Must be directly compatible with \code{\link[lubridate]{as_date}} and preferably of the form \code{"yyyy-mm-dd"} or \code{"\%Y-\%m-\%d"}.
#' @param group_by Group the samples by a variable in the metadata.
#' @param split Split the plot into subplots of each taxa. (\emph{default:} \code{FALSE}) 
#' @param tax_aggregate The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"OTU"})
#' @param tax_add Additional taxonomic level(s) to display, e.g. \code{"Phylum"}. (\emph{default:} \code{"none"})
#' @param tax_show The number of taxa to show, or a vector of taxa names. (\emph{default:} \code{6})
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible. 
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param raw (\emph{logical}) Display raw input instead of converting to percentages. (\emph{default:} \code{FALSE}) 
#' @param plotly (\emph{logical}) Returns an interactive plot instead. (\emph{default:} \code{FALSE})
#' @param ... Additional arguments passed to \code{\link[lubridate]{as_date}} to make the time_variable compatible with the timeseries plot, fx the \code{format} or \code{tz} arguments, see \code{?as_date}.
#' 
#' @keywords timeseries
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @import lubridate
#' @import plotly
#' @import data.table
#' @return A ggplot2 object.
#' 
#' @export 
#' @examples 
#' 
#' #Load example data
#' data("AalborgWWTPs")
#' 
#' #Timeseries of the 5 most abundant OTUs based on the "Date" column
#' amp_timeseries(AalborgWWTPs,
#'                time_variable = "Date",
#'                tax_aggregate = "OTU")
#' 
#' #As the above warning suggests, there are more than one sample per date in the data,
#' #in this case one from Aalborg East and one from Aalborg West. The average of the 
#' #two samples is then shown per date. In such case it is then recommended to either 
#' #subset the data, or group the samples by setting group_by = "" and split by tax_aggregate
#' #by setting split = TRUE:
#' amp_timeseries(AalborgWWTPs,
#'                time_variable = "Date",
#'                group_by = "Plant", 
#'                split = TRUE, 
#'                tax_show = 6,
#'                tax_aggregate = "Genus",
#'                tax_add = "Phylum")
#' 
#' @author Julie Klessner Thun Pedersen \email{julieklessnerthun@@gmail.com}
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}

amp_timeseries <- function(data,
                           time_variable = NULL,
                           group_by = NULL,
                           tax_aggregate = "OTU",
                           tax_add = NULL,
                           tax_show = 5,
                           tax_class = NULL,
                           tax_empty = "best",
                           split = FALSE,
                           raw = FALSE,
                           plotly = FALSE,
                           ...) {
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data,
                     tax_class = tax_class,
                     tax_empty = tax_empty, 
                     tax_level = tax_aggregate)
  
  #tax_add and tax_aggregate can't be the same
  if(!is.null(tax_aggregate) & !is.null(tax_add)) {
    if(tax_aggregate == tax_add) {
      stop("tax_aggregate and tax_add cannot be the same")
    }
  }
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  metadata <- data[["metadata"]]
  
  ## try to find a date column
  if(is.null(time_variable)) {
    dateCols <- unlist(lapply(metadata, lubridate::is.Date))
    if(sum(dateCols) == 1) {
      time_variable <- colnames(metadata)[which(dateCols)]
      message("No \"time_variable\" provided, assuming the column \"", time_variable, "\" contains the dates.\n")
    } else {
      stop("Please provide a valid date column by the argument time_variable.")
    }
  }
  
  ## Coerce the group_by and facet_by variables to factor to always be considered categorical. Fx Year is automatically loaded as numeric by R, but it should be considered categorical. 
  if(!is.null(group_by)) {
    metadata[group_by] <- lapply(metadata[group_by], factor)
  }
  
  if (raw == FALSE){
    abund <- as.data.frame(sapply(abund, function(x) x/sum(x)*100))
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
    tidyr::gather(key = Sample, value = Abundance, -Display) %>% as.data.table()
  
  abund3 <- abund3[, "sum":=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    as.data.frame()
  
  suppressWarnings(
    if (!is.null(group_by)){
      if (length(group_by) > 1){
        grp <- data.frame(Sample = metadata[,1], Group = apply(metadata[,group_by], 1, paste, collapse = "; ")) 
        oldGroup <- unique(cbind.data.frame(metadata[,group_by], Group = grp$Group))
      } else{
        grp <- data.frame(Sample = metadata[,1], Group = metadata[,group_by]) 
      }
      abund3$Group <- grp$Group[match(abund3$Sample, grp$Sample)]
      abund5 <- abund3
    } else {
      abund5 <- data.frame(abund3, Group = abund3$Sample)
    }
  )
  
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
  abund7 <- as.data.frame(abund7)
  abund7$Display <- factor(abund7$Display, levels = TotalCounts$Display)
  
  if (length(group_by) > 1) {
    abund7 <- merge(abund7, oldGroup)
  }
  
  abund7 <- merge(abund7, metadata, by.x = "Sample", by.y = colnames(metadata)[1])
  abund7[,time_variable] <- lubridate::as_date(abund7[,time_variable], ...)
  abund7$DisplayGroup <- paste(abund7$Display, abund7$Group)
  colnames(abund7)[which(colnames(abund7) == "Display")] <- tax_aggregate
  colnames(abund7)[which(colnames(abund7) == "sum")] <- "Value"
  abund7 <- abund7[,-c(which(colnames(abund7) == "Abundance"))]
  abund7 <- unique(abund7)
  
  if(is.null(group_by)) {
    if(any(duplicated(metadata[,time_variable]))) {
      warning("Duplicate dates in column ", time_variable, ", displaying the average for each date.\n Consider grouping dates using the group_by argument or subset the data using amp_subset_samples.\n")
      abund7 %>% 
        dplyr::group_by_(time_variable, tax_aggregate) %>% 
        dplyr::summarise_at("Value", mean, na.rm = TRUE) -> abund7
    }
    if(isTRUE(split)) {
      p <- ggplot(abund7, aes_string(x=time_variable, y="Value"))
    } else if(!isTRUE(split)) {
      p <- ggplot(abund7, aes_string(x=time_variable, y="Value", color = tax_aggregate))
    }
  } else if(!is.null(group_by)) {
    if(isTRUE(split)) {
      p <- ggplot(abund7, aes_string(x=time_variable, y="Value", color = "Group", group = "DisplayGroup"))
    } else if(!isTRUE(split)) {
      p <- ggplot(abund7, aes_string(x=time_variable, y="Value", color = "Group", group = "DisplayGroup", linetype = tax_aggregate, shape = tax_aggregate))
    }
  }
  p <- p +
    geom_line() +
    geom_point() +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 10, vjust = 0.3, angle = 90),
          panel.grid.major.x = element_line(color = "grey90"),
          panel.grid.major.y = element_line(color = "grey90"))
  
  if(isTRUE(raw)) {
    p <- p + ylab("Read counts")
  } else if(!isTRUE(raw)) {
    p <- p + ylab("Read abundance (%)")
  }
  
  if(split == T){
    p <- p + facet_wrap(tax_aggregate) +
      theme(strip.background = element_rect(colour=NA, fill="grey95"),
            panel.grid.major.x = element_line(color = "grey90"),
            panel.grid.major.y = element_line(color = "grey90"),
            legend.position = "bottom",
            strip.text = element_text(size = 10),
            legend.title=element_blank())
  }
  
  if(isTRUE(plotly)) {
    return(ggplotly(p, tooltip = c("Group", tax_aggregate, "Date", "Value")))
  } else if(!isTRUE(plotly)) {
    return(p)
  }
}
