#' Heatmap
#'
#' Generates a heatmap of amplicon data by using sample metadata to aggregate samples and taxonomy to aggregate OTUs.
#'
#' @usage amp_heatmap(data, group_by = "")
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param group_by (\emph{recommended}) Group the samples by a categorical variable in the metadata. If \code{NULL} then all samples are shown.
#' @param facet_by Facet the samples by a categorical variable in the metadata. 
#' @param scale_by Scale the abundances by a variable in the metadata.
#' @param normalise_by A variable or a specific sample in the metadata to normalise the counts by.
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
#' @param measure Calculate and display either \code{"mean"}, \code{"max"} or \code{"median"} across the groups. (\emph{default:} \code{"mean"})
#' @param sort_by Sort the heatmap by a specific value of the \code{"group_by"} argument, e.g. \code{"Treatment A"}.
#' @param order_x_by A sample or vector to order the y-axis by, or \code{"cluster"} for hierarchical clustering by \code{\link[stats]{hclust}}.
#' @param order_y_by A taxonomy group or vector to order the x-axis by, or \code{"cluster"} for hierarchical clustering by \code{\link[stats]{hclust}}.
#' @param plot_values (\emph{logical}) Plot the values on the heatmap or not. (\emph{default:} \code{TRUE})
#' @param plot_values_size The size of the plotted values. (\emph{default:} \code{4})
#' @param plot_legendbreaks A vector of breaks for the abundance legend, fx \code{c(1, 10, 20)}.
#' @param plot_colorscale The type of scale used for the coloring of abundances, either \code{"sqrt"} or \code{"log10"}. (\emph{default:} \code{"log10"})
#' @param plot_na (\emph{logical}) Whether to color missing values with the lowest color in the scale or not. (\emph{default:} \code{TRUE})
#' @param min_abundance All values below this value are given the same color. (\emph{default:} \code{0.1})
#' @param max_abundance All values above this value are given the same color.
#' @param color_vector Vector of colors for the colorscale, e.g. \code{c("white", "red")}.
#' @param round Number of digits to show with the values. (\emph{default:} \code{1})
#' @param normalise (\emph{logical}) Transform the OTU read counts to be in percent per sample. (\emph{default:} \code{TRUE})
#' @param textmap (\emph{logical}) Return a data frame to print as raw text instead of a ggplot2 object. (\emph{default:} \code{FALSE})
#' @param plot_functions Return a 2-column grid plot instead, showing known functional information about the Genus-level OTUs next to the heatmap. When using this feature, make sure that either \code{tax_aggregate} is set to "Genus" or that \code{tax_add} contains "Genus". (\emph{default:} \code{FALSE})
#' @param function_data A data frame with functional information about genus-level OTUs in each column. If \code{NULL} the \code{data("MiF")} dataset will be used. (\emph{default:} \code{NULL})
#' @param functions A vector with the functions to be displayed. (\emph{default:} \code{c("MiDAS","FIL", "AOB", "NOB", "PAO", "GAO")})
#' @param functions_point_size Size of the plotted points in the function grid. (\emph{default:} \code{5})
#' @param rel_widths A vector with the relative widths of the heatmap and function grid when \code{plot_functions = TRUE}. (\emph{default:} \code{c(0.75, 0.25)})
#' 
#' @return A ggplot2 object, or a data frame if \code{textmap = TRUE}.
#' 
#' @export
#' 
#' @section Preserving relative abundances in a subset of larger data:
#' By default the raw read counts in the abundance matrix are normalised (transformed to percentages) by some plotting functions automatically (for example \code{\link{amp_heatmap}}, \code{\link{amp_timeseries}}, and more). This means that the relative abundances shown will be calculated based on the remaining taxa after the subset, not including the removed taxa, if any. To circumvent this, set \code{normalise = TRUE} when subsetting with the \code{\link{amp_subset_taxa}} and \code{\link{amp_subset_samples}} functions, and then set \code{normalise = FALSE} in the plotting function. This will transform the OTU counts to relative abundances BEFORE the subset, and setting \code{normalise = FALSE} will skip the transformation in the plotting function, see the example below.
#' 
#' \preformatted{
#' data("MiDAS")
#' subsettedData <- amp_subset_samples(MiDAS,
#'                                     Plant \%in\% c("Aalborg West", "Aalborg East"),
#'                                     normalise = TRUE
#'                                     )
#' amp_heatmap(subsettedData,
#'             group_by = "Plant",
#'             tax_aggregate = "Phylum",
#'             tax_add = "Genus",
#'             normalise = FALSE
#'             )
#' }
#' 
#' @section Accessing detailed raw data:
#' The complete raw data used to generate any ggplot can always be accessed with \code{ggplot2_object$data} when the plot is saved as a ggplot2 object. Additionally, a "textmap" version of the generated heatmap can also be generated by setting \code{textmap = TRUE} to only extract the raw data as shown on the particular heatmap, see examples.
#' 
#' @examples 
#' #Load example data
#' data("AalborgWWTPs")
#' 
#' #Heatmap grouped by WWTP
#' amp_heatmap(AalborgWWTPs, group_by = "Plant")
#' 
#' #Heatmap of 20 most abundant Genera (by mean) grouped by WWTP, split by Year, 
#' #values not plotted for visibility, phylum name added and colorscale adjusted manually
#' amp_heatmap(AalborgWWTPs,
#'            group_by = "Plant", 
#'            facet_by = "Year",
#'            plot_values = FALSE,
#'            tax_show = 20,
#'            tax_aggregate = "Genus", 
#'            tax_add = "Phylum",
#'            color_vector = c("white", "red"), 
#'            plot_colorscale = "sqrt",
#'            plot_legendbreaks = c(1, 5, 10)
#'            )
#'            
#' #Heatmap with known functional information about the Genera shown to the right
#' amp_heatmap(AalborgWWTPs, 
#'             group_by = "Plant", 
#'             tax_aggregate = "Genus", 
#'             plot_functions = TRUE, 
#'             functions = c("PAO", "GAO", "AOB", "NOB")
#'             )
#'            
#' #A raw text version of the heatmap can be printed or saved as a data frame with textmap = TRUE:
#' textmap <- amp_heatmap(AalborgWWTPs, 
#'                        group_by = "Plant", 
#'                        tax_aggregate = "Genus",
#'                        plot_functions = TRUE,
#'                        functions = c("PAO", "GAO", "AOB", "NOB"),
#'                        textmap = TRUE
#'                        )
#' textmap
#' 
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr filter desc arrange group_by mutate summarise
#' @importFrom tidyr gather spread
#' @importFrom data.table as.data.table data.table setkey dcast melt
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales squish
#' @importFrom cowplot plot_grid
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_heatmap <- function(data,
                        group_by = NULL,
                        facet_by = NULL,
                        normalise = TRUE,
                        tax_aggregate = "Phylum",
                        tax_add = NULL,
                        tax_show = 10,
                        tax_class = NULL,
                        tax_empty = "best",
                        order_x_by = NULL,
                        order_y_by = NULL,
                        plot_values = TRUE,
                        plot_values_size = 4,
                        plot_legendbreaks = NULL,
                        plot_colorscale = "log10", 
                        plot_na = TRUE, 
                        measure = "mean",
                        min_abundance = 0.1,
                        max_abundance = NULL, 
                        sort_by = NULL, 
                        normalise_by = NULL,
                        scale_by = NULL,
                        color_vector = NULL,
                        round = 1,
                        textmap = FALSE,
                        plot_functions = FALSE,
                        function_data = NULL, 
                        functions = c("MiDAS","FIL", "AOB", "NOB", "PAO", "GAO"),
                        functions_point_size = 5,
                        rel_widths = c(0.75, 0.25)
                        ) {
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis2 functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)", call. = FALSE)
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data,
                     tax_class = tax_class,
                     tax_empty = tax_empty, 
                     tax_level = tax_aggregate)
  
  #tax_add and tax_aggregate can't be the same
  if(!is.null(tax_aggregate) & !is.null(tax_add)) {
    if(identical(tax_aggregate, tax_add)) {
      stop("tax_aggregate and tax_add cannot be the same", call. = FALSE)
    }
  }
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  metadata <- data[["metadata"]]
  
  #add functions check
  if(isTRUE(plot_functions)) {
    if(!any("Genus" %in% c(tax_add, tax_aggregate)))
      stop("One of the arguments tax_add or tax_aggregate must contain \"Genus\"", call. = FALSE)
  }
  # Retrieve the function data if not provided
  if (is.null(function_data)) {
    data(MiF, envir = environment())
    function_data <- MiF
  }
  
  ## Coerce the group_by and facet_by variables to factor to always be considered categorical. Fx Year is automatically loaded as numeric by R, but it should be considered categorical. 
  ## Grouping a heatmap by a continuous variable doesn't make sense 
  if(!is.null(group_by)) {
    metadata[group_by] <- lapply(metadata[group_by], factor)
  }
  
  if(!is.null(facet_by)) {
    if(is.null(group_by)) {
      group_by <- facet_by
    }
    metadata[facet_by] <- lapply(metadata[facet_by], factor)
  }
  
  ## Scale the data by a selected metadata sample variable
  if (!is.null(scale_by)){
    variable <- as.numeric(metadata[,scale_by])
    abund <- t(t(abund)*variable)
  }
  
  if (isTRUE(normalise)) {
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
    tidyr::gather(key = Sample, value = Abundance, -Display) %>% as.data.table()
  
  abund3 <- abund3[, "sum":=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    as.data.frame()
  
  ## Add group information
  
  if(!is.null(facet_by)){
    ogroup <- group_by
    group_by <- c(group_by, facet_by)
  }
  
  suppressWarnings(
    if (!is.null(group_by)){
      if (length(group_by) > 1){
        grp <- data.frame(Sample = metadata[,1], Group = apply(metadata[,group_by], 1, paste, collapse = " ")) 
        oldGroup <- unique(cbind.data.frame(metadata[,group_by], Group = grp$Group))
      } else{
        grp <- data.frame(Sample = metadata[,1], Group = metadata[,group_by]) 
      }
      abund3$Group <- grp$Group[match(abund3$Sample, grp$Sample)]
      abund5 <- abund3
    } else{ abund5 <- data.frame(abund3, Group = abund3$Sample)}
  )
  
  ## Take the average to group level
  
  if (measure == "mean"){
    abund6 <- data.table(abund5)[, Abundance:=mean(sum), by=list(Display, Group)] %>%
      setkey(Display, Group) %>%
      unique() %>% 
      as.data.frame()
    TotalCounts <- group_by(abund6, Display) %>%
      summarise(Abundance = sum(Abundance)) %>%
      arrange(desc(Abundance))
  }
  
  if (measure == "max"){
    abund6 <- data.table(abund5)[, Abundance:=max(sum), by=list(Display, Group)] %>%
      setkey(Display, Group) %>%
      unique() %>% 
      as.data.frame()
    TotalCounts <- group_by(abund6, Display) %>%
      summarise(Abundance = max(Abundance)) %>%
      arrange(desc(Abundance))
  }  
  
  if (measure == "median"){
    abund6 <- data.table(abund5)[, Abundance:=median(sum), by=list(Display, Group)] %>%
      setkey(Display, Group) %>%
      unique() %>% 
      as.data.frame()
    TotalCounts <- group_by(abund6, Display) %>%
      summarise(Abundance = median(Abundance)) %>%
      arrange(desc(Abundance))
  }
  
  if(!is.null(sort_by)) {
    if(is.null(group_by)) {
      if(!any(sort_by %in% abund6$Sample))
        stop("Can't find \"", sort_by, "\" among sample names", call. = FALSE)
      TotalCounts <- filter(abund6, Sample == sort_by) %>%
        arrange(desc(Abundance))
    } else if(!is.null(group_by)) {
      if(any(sort_by %in% abund6$Group) & any(sort_by %in% abund6$Sample)) {
        stop(paste0(sort_by, " is both found among samples and in the group_by variable (", group_by, "). Cannot sort by both a sample and a group."), call. = FALSE)
      } else if(!any(sort_by %in% abund6$Group) & any(sort_by %in% abund6$Sample)) {
        TotalCounts <- filter(abund6, Sample == sort_by) %>%
          arrange(desc(Abundance))
      } else if(any(sort_by %in% abund6$Group) & !any(sort_by %in% abund6$Sample)) {
        TotalCounts <- filter(abund6, Group == sort_by) %>%
          arrange(desc(Abundance))
      } else if(!any(sort_by %in% abund6$Group) & !any(sort_by %in% abund6$Sample))
        stop("Can't find \"", sort_by, "\" among sample or group names", call. = FALSE)
    }
  }
  
  ## Subset to X most abundant levels
  if (is.numeric(tax_show)){
    if (tax_show > nrow(TotalCounts)){
      tax_show <- nrow(TotalCounts)
    }
    abund7 <- filter(abund6, Display %in% unique(TotalCounts$Display)[1:tax_show])
  }
  
  ## Subset to a list of level names
  if (!is.numeric(tax_show)){
    if (tax_show != "all"){
      abund7 <- filter(abund6, Display %in% tax_show)
    }
    ### Or just show all  
    if (tax_show == "all"){
      tax_show <- nrow(TotalCounts)  
      abund7 <- filter(abund6, Display %in% TotalCounts$Display[1:tax_show]) 
    }
  }
  abund7 <- as.data.frame(abund7)
  
  ## Normalise to a specific group or sample (The Abundance of the group is set as 1)  
  if(!is.null(normalise_by)){
    temp <- data.table::dcast(abund7, Display~Group, value.var = "Abundance")
    temp1 <- cbind.data.frame(Display = temp$Display, temp[,-1]/temp[,normalise_by])   
    abund7 <- data.table::melt(temp1, id.var = "Display", value.name="Abundance", variable.name="Group")
  } 
  
  ## Order.y
  if (is.null(order_y_by)){
    abund7$Display <- factor(abund7$Display, levels = rev(unique(TotalCounts$Display)))
  }
  if (!is.null(order_y_by)){
    if ((length(order_y_by) == 1) && (order_y_by != "cluster")){
      temp1 <- filter(abund7, Group == order_y_by) %>%
        group_by(Display) %>%
        summarise(Mean = mean(Abundance)) %>%
        arrange(desc(Mean))
      
      abund7$Display <- factor(abund7$Display, levels = rev(temp1$Display))
    }
    if (length(order_y_by) > 1){
      abund7$Display <- factor(abund7$Display, levels = order_y_by)
    }
    if ((length(order_y_by) == 1) && (order_y_by == "cluster")){
      if (is.null(max_abundance)){max_abundance <- max(abund7$Abundance)}
      tdata <- mutate(abund7, 
                      Abundance = ifelse(Abundance < min_abundance, min_abundance, Abundance),
                      Abundance = ifelse(Abundance > max_abundance, max_abundance, Abundance))
      tdata <- data.table::dcast(tdata, Display~Group, value.var = "Abundance")
      rownames(tdata) <- tdata$Display
      tdata2 <- tdata[,-1]
      tclust <- hclust(dist(tdata2))
      tnames <- levels(droplevels(tdata$Display))[tclust$order]
      abund7$Display <- factor(abund7$Display, levels = tnames)
    }
  }
  
  ## Order.x
  if (!is.null(order_x_by)){
    if ((length(order_x_by) == 1) && (order_x_by != "cluster")){
      temp1 <- filter(abund7, Display == order_x_by) %>%
        group_by(Group) %>%
        summarise(Mean = mean(Abundance)) %>%
        arrange(desc(Mean))
      abund7$Group <- factor(abund7$Group, levels = as.character(temp1$Group))
    }    
    if (length(order_x_by) > 1){
      abund7$Group <- factor(abund7$Group, levels = order_x_by)
    }
    if ((length(order_x_by) == 1) && (order_x_by == "cluster")){
      if (is.null(max_abundance)){max_abundance <- max(abund7$Abundance)}
      tdata <- mutate(abund7, 
                      Abundance = ifelse(Abundance < min_abundance, min_abundance, Abundance),
                      Abundance = ifelse(Abundance > max_abundance, max_abundance, Abundance))
      tdata <- data.table::dcast(tdata, Display~Group, value.var = "Abundance")
      rownames(tdata) <- tdata$Display
      tdata2 <- tdata[,-1]
      tclust <- hclust(dist(t(tdata2)))
      tnames <- tclust$labels[tclust$order]
      abund7$Group <- factor(abund7$Group, levels = tnames) 
    }
  }
  
  ## Handle NA values
  if(plot_na == FALSE){ plot_na <- "grey50" }else{ if(!is.null(color_vector)) {plot_na <-color_vector[1]} else {plot_na <-"#67A9CF"}}  
  
  ## Scale to percentages if not normalised and scaled
  
  if (length(group_by) > 1 ){ abund7 <- merge(abund7, oldGroup)}
  
  if (is.null(min_abundance)){
    min_abundance <- ifelse(min(abund7$Abundance) > 0.001, min(abund7$Abundance), 0.001)
  }
  if (is.null(max_abundance)){
    max_abundance <- max(abund7$Abundance)
  }
  
  ## Define the output 
  if (!isTRUE(textmap)) {
    ## Make a heatmap style plot
    heatmap <- ggplot(abund7, aes_string(x = "Group", y = "Display", label = formatC("Abundance", format = "f", digits = 1))) +     
      geom_tile(aes(fill = Abundance), colour = "white", size = 0.5) +
      theme(axis.text.y = element_text(size = 12, color = "black", vjust = 0.4),
            axis.text.x = element_text(size = 10, color = "black", vjust = 0.5, angle = 90, hjust = 1),
            axis.title = element_blank(),
            text = element_text(size = 8, color = "black"),
            axis.line = element_blank(),
            #axis.ticks.length = unit(1, "mm"),
            plot.margin = unit(c(1,1,1,1), "mm"),
            title = element_text(size = 8),
            panel.background = element_blank())
    
    ## Get colorpalette for colorscale or set default
    if (!is.null(color_vector)){
      color.pal = color_vector
    } else {
      color.pal = rev(RColorBrewer::brewer.pal(3, "RdBu"))
    }
    
    if (plot_values == TRUE){
      abund8 <- abund7
      abund8$Abundance <- round(abund8$Abundance, round)
      heatmap <- heatmap + geom_text(data = abund8, size = plot_values_size, colour = "grey10", check_overlap = TRUE) +
        theme(legend.position = "none")
    }
    if (is.null(plot_legendbreaks)){
      heatmap <- heatmap +scale_fill_gradientn(colours = color.pal, trans = plot_colorscale, na.value=plot_na, oob = scales::squish, limits = c(min_abundance, max_abundance))
    }
    if (!is.null(plot_legendbreaks)){
      heatmap <- heatmap +scale_fill_gradientn(colours = color.pal, trans = plot_colorscale, breaks=plot_legendbreaks, na.value=plot_na , oob = scales::squish, limits = c(min_abundance, max_abundance))
    }
    
    
    if (is.null(normalise_by)){
      heatmap <- heatmap + labs(x = "", y = "", fill = "% Read\nAbundance")  
    }
    if (!is.null(normalise_by)){
      heatmap <- heatmap + labs(x = "", y = "", fill = "Relative")  
    }
    
    if(!is.null(facet_by)){
      if(length(ogroup) > 1){
        heatmap$data$Group <- apply(heatmap$data[,ogroup], 1, paste, collapse = " ")  
      } else{
        heatmap$data$Group <- heatmap$data[,ogroup]
      }
      
      if(plot_values == TRUE){
        if(length(ogroup) > 1){
          heatmap$layers[[2]]$data$Group <- apply(heatmap$layers[[2]]$data[,ogroup], 1, paste, collapse = " ")  
        } else{
          heatmap$layers[[2]]$data$Group <- heatmap$layers[[2]]$data[,ogroup]
        }
      }
      heatmap <- heatmap + facet_grid(reformulate(facet_by), scales = "free_x", space = "free")
      heatmap <- heatmap + theme(strip.text = element_text(size = 10))
    }
    
    #Return a function grid next to the heatmap with known functions about the Genera
    if(isTRUE(plot_functions)) {
      # Retrieve the genus names from the plot
      names <- data.frame(do.call('rbind', strsplit(levels(droplevels(heatmap$data$Display)),'; ',fixed=TRUE)))
      names <- data.frame(Genus = names[,which(c(tax_add, tax_aggregate) == "Genus")])
      names$Genus <- as.character(names$Genus)
      
      # Merge the genus and function information
      nameFunc <- merge(x = names, y = function_data[,c("Genus",functions)], all.x = TRUE, all.y = FALSE) 
      nameFunc[is.na(nameFunc)] <- "NT"
      nameFuncM <- data.table::melt(nameFunc, id.vars = "Genus", value.name = "Value", variable.name = "Function")
      nameFuncM$Value <- factor(nameFuncM$Value, levels = c("POS", "VAR", "NEG", "NT"))
      nameFuncM$Genus <- factor(nameFuncM$Genus, levels = names$Genus)
      
      functions_plot <- ggplot(nameFuncM, aes(x = Function, y = Genus, color = Value)) +
        geom_point(size = functions_point_size) +
        scale_color_manual(values = c("#31a354", "orange", "#f03b20", "grey90"), drop = FALSE) +
        theme(axis.text.x = element_text(size = 10, color = "black", angle = 90, hjust = 1, vjust = 0.4),
              axis.text.y = element_blank(),
              axis.title = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size = 10),
              axis.ticks.length = unit(1, "mm"),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_line(color = "grey95"),
              legend.key = element_blank()
        )
      return(cowplot::plot_grid(heatmap, functions_plot, ncol = 2, rel_widths = rel_widths, align = "h", axis = "tb"))
    } else if(!isTRUE(plot_functions)) {
      return(heatmap)
    }
  } else if (isTRUE(textmap)) {
    #raw text heatmap data frame
    textmap <- abund7[,c("Display", "Abundance", "Group"), drop = FALSE] %>% 
      unique() %>%
      spread(key = Group, value = Abundance)
    if(isTRUE(plot_functions)) {
      textmap <- merge(cbind(textmap, 
                             Genus = as.character(data.frame(do.call('rbind', strsplit(levels(droplevels(textmap$Display)), '; ', fixed=TRUE)))[,which(c(tax_add, tax_aggregate) == "Genus")])), 
                       function_data[,c("Genus", switch(isTRUE(plot_functions), functions)), drop = FALSE],
                       by = "Genus",
                       all.x = TRUE,
                       all.y = FALSE,
                       fill = NA)
    }
    textmap <- textmap %>%
      arrange(desc(droplevels(Display)))
    textmap <- data.frame(textmap[,-c(which(colnames(textmap) %in% c("Display", "Genus"))), drop = FALSE], 
                          row.names = textmap$Display, 
                          check.names = FALSE)
    return(textmap)
  }
}
