#' Generate a heatmap from amplicon data
#'
#' Generate a heatmap in ggplot2 format from amplicon data in ampvis format. Use sample metadata to aggregate samples and taxonomy to aggregate OTUs.
#'
#' @usage amp_heatmap(data)
#'
#' @param data (required) A ampvis object.
#' @param group A variable from the associated sample data to group samples by.
#' @param facet A variable to facet the plot by.
#' @param scale A variable from the associated sample data to scale the abundance by.
#' @param normalise A specific sample or group to normalise the counts to, or "relative".
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (default: "Phylum").
#' @param tax.add Additional taxonomic levels to display for each entry e.g. "Phylum" (default: "none"). 
#' @param tax.show The number of taxa to show or a vector of taxa names (default: 10).
#' @param tax.empty Either "remove" OTUs without taxonomic information, add "best" classification or add the "OTU" name (default: "best").
#' @param tax.class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param calc Calculate and display "mean", "max" or "median" across the groups (default: "mean").
#' @param sort.by Sort the heatmap by a specific value of the "group", e.g. "Treatment A".
#' @param order.x A taxonomy group or vector to order the x-axis by, alternatively "cluster".
#' @param order.y A sample or vector to order the y-axis by, alternatively "cluster".
#' @param plot.numbers Plot the values on the heatmap (default: T)
#' @param plot.breaks A vector of breaks for the abundance legend.
#' @param plot.colorscale Either "sqrt" or "log10" (default: "log10")
#' @param plot.na Whether to color missing values with the lowest color in the scale (default: T).
#' @param plot.text.size The size of the plotted text (default: 4).
#' @param min.abundance All values below are given the same color (default: 0.1).
#' @param max.abundance All values above are given the same color.
#' @param output To output a plot or the complete data inclusive dataframes (default: "plot")
#' @param color.vector Vector with colors for colorscale e.g. c("red","white") (default: NULL)
#' @param round Number of digits to plot (default: 1)
#' @param raw Display raw input instead of converting to percentages (default: F) 
#' 
#' @return A ggplot2 object or a list with the ggplot2 object and associated dataframes.
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_heatmap <- function(data, group = "Sample", facet = NULL, normalise = NULL, scale = NULL, tax.aggregate = "Phylum", tax.add = NULL, tax.show = 10, tax.class = NULL, tax.empty = "best", order.x = NULL, order.y = NULL, plot.numbers = T, plot.breaks = NULL, plot.colorscale = "log10", plot.na = T, output = "plot",plot.text.size = 4, plot.theme = "normal", calc = "mean", min.abundance = 0.1, max.abundance = NULL, sort.by = NULL, color.vector = NULL, round = 1, raw = F){
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax.class = tax.class, tax.empty = tax.empty, tax.level = tax.aggregate)
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  sample <- data[["metadata"]]
  
  ## Scale the data by a selected metadata sample variable
  if (!is.null(scale)){
    variable <- as.numeric(sample[,scale])
    abund <- t(t(abund)*variable)
  }
  
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
    gather(key = Sample, value = Abundance, -Display)
  
  abund3 <- data.table(abund3)[, sum:=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    unique() %>% 
    as.data.frame()
  
  ## Add group information
  
  if(!is.null(facet)){
    ogroup <- group
    group <- c(group, facet)
  }
  
  suppressWarnings(
    if (group != "Sample"){
      if (length(group) > 1){
        grp <- data.frame(Sample = sample$SeqID, Group = apply(sample[,group], 1, paste, collapse = " ")) 
        oldGroup <- unique(cbind.data.frame(sample[,group], Group = grp$Group))
      } else{
        grp <- data.frame(Sample = sample$SeqID, Group = sample[,group]) 
      }
      abund3$Group <- grp$Group[match(abund3$Sample, grp$Sample)]
      abund5 <- abund3
    } else{ abund5 <- data.frame(abund3, Group = abund3$Sample)}
  )
  
  ## Take the average to group level
  
  if (calc == "mean"){
    abund6 <- data.table(abund5)[, Abundance:=mean(sum), by=list(Display, Group)] %>%
      setkey(Display, Group) %>%
      unique() %>% 
      as.data.frame()
  }
  
  if (calc == "max"){
    abund6 <- data.table(abund5)[, Abundance:=max(sum), by=list(Display, Group)] %>%
      setkey(Display, Group) %>%
      unique() %>% 
      as.data.frame()
  }  
  
  if (calc == "median"){
    abund6 <- data.table(abund5)[, Abundance:=median(sum), by=list(Display, Group)] %>%
              setkey(Display, Group) %>%
              unique() %>% 
              as.data.frame()
  }
  
  
  ## Find the X most abundant levels
  if (calc == "mean"){
    TotalCounts <- group_by(abund6, Display) %>%
      summarise(Abundance = sum(Abundance)) %>%
      arrange(desc(Abundance))
  }
  
  if (calc == "max"){
    TotalCounts <- group_by(abund6, Display) %>%
      summarise(Abundance = max(Abundance)) %>%
      arrange(desc(Abundance))
  }
  
  if (calc == "median"){
    TotalCounts <- group_by(abund6, Display) %>%
      summarise(Abundance = median(Abundance)) %>%
      arrange(desc(Abundance))
  }
  
  if (!is.null(sort.by)){
    TotalCounts <- filter(abund6, Group == sort.by) %>%
      arrange(desc(Abundance))
  }
  
  ## Subset to X most abundant levels
  if (is.numeric(tax.show)){
    if (tax.show > nrow(TotalCounts)){  
      tax.show <- nrow(TotalCounts)
    }
    abund7 <- filter(abund6, Display %in% TotalCounts$Display[1:tax.show])
  }
  
  ## Subset to a list of level names
  if (!is.numeric(tax.show)){
    if (tax.show != "all"){
      abund7 <- filter(abund6, Display %in% tax.show)    
    }
    ### Or just show all  
    if (tax.show == "all"){
      tax.show <- nrow(TotalCounts)  
      abund7 <- filter(abund6, Display %in% TotalCounts$Display[1:tax.show]) 
    }
  }
  abund7 <- as.data.frame(abund7)
  
  ## Normalise to a specific group (The Abundance of the group is set as 1)  
  
  if(!is.null(normalise)){
    if (normalise != "relative"){
      #temp <- dcast(abunds7, Display~Group, value.var = "Abundance") # Is this working?
      temp <- spread(abund7, key = Group, value = Abundance) 
      temp1 <- cbind.data.frame(Display = temp$Display, temp[,-1]/temp[,normalise])   
      abund7 <- gather(temp1, key = Group, value = Abundance, -Display)
    }
  } 
  if(!is.null(normalise)){
    if (normalise == "relative"){
      #temp <- dcast(abund7, Display~Group, value.var = "Abundance") # Is this working?
      temp <- spread(abund7, key = Group, value = Abundance) 
      temp1 <- cbind.data.frame(Display = temp[,1], temp[,-1]/apply(as.matrix(temp[,-1]), 1, mean))    
      abund7 <- gather(temp1, key = Group, value = Abundance, -Display)
    }
  }
  
  ## Order.y
  if (is.null(order.y)){
    abund7$Display <- factor(abund7$Display, levels = rev(TotalCounts$Display))
  }
  if (!is.null(order.y)){
    if ((length(order.y) == 1) && (order.y != "cluster")){       
      temp1 <- filter(abund7, Group == order.y) %>%
        group_by(Display) %>%
        summarise(Mean = mean(Abundance)) %>%
        arrange(desc(Mean))
      
      abund7$Display <- factor(abund7$Display, levels = rev(temp1$Display))
    }
    if (length(order.y) > 1){
      abund7$Display <- factor(abund7$Display, levels = order.y)
    }
    if ((length(order.y) == 1) && (order.y == "cluster")){
      if (is.null(max.abundance)){max.abundance <- max(abund7$Abundance)}
      tdata <- mutate(abund7, 
                      Abundance = ifelse(Abundance < min.abundance, min.abundance, Abundance),
                      Abundance = ifelse(Abundance > max.abundance, max.abundance, Abundance))
      tdata <- dcast(tdata, Display~Group, value.var = "Abundance") #is this working?
      #tdata <- spread(tdata, key = Group, value = Abundance) 
      rownames(tdata) <- tdata$Display
      tdata2 <- tdata[,-1]
      tclust <- hclust(dist(tdata2))
      tnames <- levels(droplevels(tdata$Display))[tclust$order]
      abund7$Display <- factor(abund7$Display, levels = tnames)
    }
  }
  
  ## Order.x
  if (!is.null(order.x)){
    if ((length(order.x) == 1) && (order.x != "cluster")){
      temp1 <- filter(abund7, Display == order.x) %>%
        group_by(Group) %>%
        summarise(Mean = mean(Abundance)) %>%
        arrange(desc(Mean))
      abund7$Group <- factor(abund7$Group, levels = as.character(temp1$Group))
    }    
    if (length(order.x) > 1){
      abund7$Group <- factor(abund7$Group, levels = order.x)
    }
    if ((length(order.x) == 1) && (order.x == "cluster")){
      if (is.null(max.abundance)){max.abundance <- max(abund7$Abundance)}
      tdata <- mutate(abund7, 
                      Abundance = ifelse(Abundance < min.abundance, min.abundance, Abundance),
                      Abundance = ifelse(Abundance > max.abundance, max.abundance, Abundance))
      tdata <- dcast(tdata, Display~Group, value.var = "Abundance") #is this working?
      #tdata <- spread(tdata, key = Group, value = Abundance) 
      rownames(tdata) <- tdata$Display
      tdata2 <- tdata[,-1]
      tclust <- hclust(dist(t(tdata2)))
      tnames <- tclust$labels[tclust$order]
      abund7$Group <- factor(abund7$Group, levels = tnames) 
    }
  }
  
  ## Handle NA values
  if(plot.na == F){ plot.na <- "grey50" }else{ if(!is.null(color.vector)) {plot.na <-color.vector[1]} else {plot.na <-"#67A9CF"}}  
  
  ## Scale to percentages if not normalised and scaled
  
  if (length(group) > 1 ){ abund7 <- merge(abund7, oldGroup)}
  
  if (is.null(min.abundance)){
    min.abundance <- ifelse(min(abund7$Abundance) > 0.001, min(abund7$Abundance), 0.001)
    }
  if (is.null(max.abundance)){
    max.abundance <- max(abund7$Abundance)
    }
  
  ## Make a heatmap style plot
  p <- ggplot(abund7, aes_string(x = "Group", y = "Display", label = formatC("Abundance", format = "f", digits = 1))) +     
    geom_tile(aes(fill = Abundance), colour = "white", size = 0.5) +
    theme(axis.text.y = element_text(size = 12, color = "black", vjust = 0.4),
          axis.text.x = element_text(size = 10, color = "black", vjust = 0.5, angle = 90),
          axis.title = element_blank(),
          text = element_text(size = 8, color = "black"),
          axis.ticks.length = unit(1, "mm"),
          plot.margin = unit(c(1,1,1,1), "mm"),
          title = element_text(size = 8))
  
  ## Get colorpalette for colorscale or set default
  if (!is.null(color.vector)){
    color.pal = color.vector
  } else {
    color.pal = rev(brewer.pal(3, "RdBu"))
  }
  
  if (plot.numbers == T){
    abund8 <- abund7
    abund8$Abundance <- round(abund8$Abundance, round)
    p <- p + geom_text(data = abund8, size = plot.text.size, colour = "grey10", check_overlap = TRUE) +
      theme(legend.position = "none")
  }
  if (is.null(plot.breaks)){
    p <- p +scale_fill_gradientn(colours = color.pal, trans = plot.colorscale, na.value=plot.na, oob = squish, limits = c(min.abundance, max.abundance))
  }
  if (!is.null(plot.breaks)){
    p <- p +scale_fill_gradientn(colours = color.pal, trans = plot.colorscale, breaks=plot.breaks, na.value=plot.na , oob = squish, limits = c(min.abundance, max.abundance))
  }
  
  
  if (is.null(normalise)){
    p <- p + labs(x = "", y = "", fill = "% Read\nAbundance")  
  }
  if (!is.null(normalise)){
    p <- p + labs(x = "", y = "", fill = "Relative")  
  }
  
  if(!is.null(facet)){
    if(length(ogroup) > 1){
      p$data$Group <- apply(p$data[,ogroup], 1, paste, collapse = " ")  
    } else{
      p$data$Group <- p$data[,ogroup]
    }
    
    if(plot.numbers == T){
      if(length(ogroup) > 1){
        p$layers[[2]]$data$Group <- apply(p$layers[[2]]$data[,ogroup], 1, paste, collapse = " ")  
      } else{
        p$layers[[2]]$data$Group <- p$layers[[2]]$data[,ogroup]
      }
    }
    p <- p + facet_grid(reformulate(facet), scales = "free_x", space = "free")
  }
  
  
  ## Define the output 
  if (output == "complete"){
    outlist <- list(heatmap = p, data = abund7)
    return(outlist)  
  }
  if (output == "plot"){
    return(p)
  }
}
