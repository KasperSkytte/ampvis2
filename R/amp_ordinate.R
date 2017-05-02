#' Generates a ggplot2 style ordinate plot of amplicon data
#'
#' A wrapper around the vegan package to make beautiful ggplot2 ordination plots to analyse and compare microbial community compositions. The input data must be a list loaded with amp_load() from the ampvis package, but any OTU table-like matrix (OTU's in rows, sampleID's in columns and abundances in the corresponding cells, aka a contingency table) alongside metadata for the same samples can be used. Simply choose an ordination type and a plot is returned.
#'
#' @usage amp_ordinate(data)
#'
#' @param data (required) Data list loaded with amp_load() containing the elements: abund, metadata and tax.
#' @param filter_species Remove low abundant OTU's across all samples below this threshold in percent. Recommended minimum: 0.1 pct (default: 0.1). 
#' @param type Ordination type; Principal Components Analysis(PCA), Redundancy Analysis(RDA), non-metric Multidimensional Scaling(NMDS), metric Multidimensional Scaling(MMDS, aka PCoA), Correspondence Analysis(CA), Canonical Correspondence Analysis(CCA) or Detrended Correspondence Analysis(DCA) (Default: PCA)
#' @param metric Distance metric used for the distance-based ordination methods (nMDS/PCoA), any of the following: "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "jsd" (Jensen-Shannon Divergence), "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis" or simply "none" or "sqrt". See details in ?vegdist, and for JSD http://enterotype.embl.de/enterotypes.html.
#' @param transform Transform the abundance table with the decostand() function, fx "normalize", "chi.square", "hellinger" or "sqrt", see details in ?decostand. Using the hellinger transformation is always a good choice and is recommended for PCA/RDA/nMDS/PCoA to obtain a more ecologically meaningful result (learn about the double-zero problem).
#' @param constrain Variable(s) in the metadata for constrained analyses (RDA and CCA). Multiple variables can be provided by a vector, fx c("Year", "Temperature"), but keep in mind the more variables the more the result will be similar to unconstrained analysis.
#' @param x_axis Which axis from the ordination results to plot as the first axis. Have a look at the $screeplot with output = "detailed" to validate axes (default: 1).
#' @param y_axis Which axis from the ordination results to plot as the second axis. Have a look at the $screeplot with output = "detailed" to validate axes (default: 2).
#' 
#' @param sample_color Color sample points by a variable in the metadata.   
#' @param sample_color_order Order the colors in color by a vector. 
#' @param sample_label Label sample points by a variable in the metadata.
#' @param sample_label_size Sample labels text size (default: 4).
#' @param sample_label_segment_color Sample labels repel-segment color (default: "black").
#' @param sample_shape Shape sample points by a variable in the metadata.       
#' @param sample_colorframe Frame the points with a polygon colored by the color argument (default: F).
#' @param sample_colorframe_label Label by a variable in the metadata.
#' @param sample_trajectory Make a trajectory between sample points by a variable in the metadata.
#' @param sample_trajectory_group Make a trajectory between sample points by the trajectory argument, but within individual groups.
#' @param sample_plotly Enable sample point howevering to display specific metadata or "all" (default: none).
#' 
#' @param species_plot Plot species points (default: F).
#' @param species_nlabels Number of most extreme species labels to plot.
#' @param species_label_taxonomy Taxonomic level by which to label the species points (default: "Genus").
#' @param species_label_size Size of the species text labels (default: 3).
#' @param species_label_color Color of the species text labels (default: "grey10).
#' @param species_rescale rescale species (default: F).
#' @param species_size Size of the species points (default: 2).
#' @param species_shape The shape of the points, fx 1 for hollow circles or 20 for dots (default: 20).
#' @param species_plotly Enable species specific point howevering to display taxonomy (default: F).
#' 
#' @param envfit_factor A vector of factor variables from the sample data used for envfit to the model
#' @param envfit_numeric A vector of numerical variables from the sample data used for envfit to the model
#' @param envfit_signif_level The significance treshold for displaying envfit parameters (default: 0.001).
#' @param envfit_textsize Size of the envfit text on the plot (default: 3).
#' @param envfit_color Color of the envfit text on the plot (default: "darkred").
#' @param envfit_numeric_arrows_scale Scale the size of the numeric arrows (default: 1).
#' @param envfit_show Show the results on the plot (default: T).
#' 
#' @param repel Repel all labels to prevent cluttering of plots (default: T).
#' @param opacity Opacity of all plotted points and sample_colorframe opacity, 0:invisible, 1:opaque (default: 0.8).
#' @param tax.empty Option to add "best" classification or just the "OTU" name to each "OTU" (default: best).
#' @param output "plot" or "detailed"; output as list with additional information(model, scores, inputmatrix etc) or just the plot (default: "plot).
#' @param ... Pass additional arguments to the vegan ordination functions, fx rda(...), cca(...), metaMDS(...), see vegan help.
#'       
#'         
#' @return A ggplot2 object or a list with the ggplot2 object and associated dataframes.
#' 
#' @export
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_ordinate<- function(data, filter_species = 0.1, type = "PCA", metric = NULL, transform = NULL, constrain = NULL, x_axis = 1, y_axis = 2, 
                        sample_color = NULL, sample_color_order = NULL, sample_shape = NULL, sample_colorframe = FALSE, sample_colorframe_label = NULL, 
                        sample_label = NULL, sample_label_size = 4, sample_label_segment_color = "black", sample_trajectory = NULL, sample_trajectory_group = sample_trajectory, sample_plotly = NULL,
                        species_plot = FALSE, species_nlabels = 0, species_label_taxonomy = "Genus", species_label_size = 3, species_label_color = "grey10", species_rescale = FALSE, species_size = 2, species_shape = 20, species_plotly = F,
                        envfit_factor = NULL, envfit_numeric = NULL, envfit_signif_level = 0.001, envfit_textsize = 3, envfit_color = "darkred", envfit_numeric_arrows_scale = 1, envfit_show = TRUE, 
                        repel = T, opacity = 0.8, tax.empty = "best", output = "plot", ...) {
  
  #Sanity check of options
  if(species_plotly == T & !is.null(sample_plotly)){
    stop("You can not use plotly for both species and samples in the same plot.")
  }
  if(species_plotly == T | !is.null(sample_plotly)){
    warning("Forcing repel = F in order to plotly to work.")
    repel <- F
  }
  
  #Check the data
  data <- amp_rename(data = data, tax.empty = tax.empty)
  
  #First transform to percentages
  abund_pct <- as.data.frame(sapply(data$abund, function(x) x/sum(x) * 100))
  rownames(abund_pct) <- rownames(data$abund) #keep rownames
  data$abund <- abund_pct
  
  #Then filter low abundant OTU's where ALL samples have below the threshold set with filter_species in percent
  data$abund <- data$abund[!apply(data$abund, 1, function(row) all(row <= filter_species)),] #remove low abundant OTU's 
  rownames(data$tax) <- data$tax$OTU
  data$tax <- data$tax[rownames(data$abund),] #same with taxonomy
  data$metadata <- data$metadata[colnames(data$abund),] #same with metadata
  
  #to fix user argument characters, so fx PCoA/PCOA/pcoa are all valid
  type <- tolower(type)
  output <- tolower(output)
  
  if(!is.null(metric)) {
    metric <- tolower(metric)
  } else if(is.null(metric)) {
    if(type == "nmds" | type == "mmds" | type == "pcoa" | type == "dca") {
      warning("No distance metric selected, using raw data. If this is not deliberate, please provide one with the argument: metric")
    }
    metric <- "none"
  }
  
  #data transformation with decostand()
  if(!is.null(transform)) {
    transform <- tolower(transform)
    if(transform == "sqrt") {
      data$abund <- t(sqrt(t(data$abund)))
    } else {
      data$abund <- t(decostand(t(data$abund), method = transform))
    }
  } 
  
  #Calculate distance matrix with vegdist()
  if (metric == "none") {
    inputmatrix <- t(data$abund)
  } else if(metric == "jsd") {
    #This is based on http://enterotype.embl.de/enterotypes.html
    #Abundances of 0 will be set to the pseudocount value to avoid 0-value denominators
    dist.JSD <- function(inMatrix, pseudocount=0.000001) {
      KLD <- function(x,y) sum(x *log(x/y))
      JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
      matrixColSize <- length(colnames(inMatrix))
      matrixRowSize <- length(rownames(inMatrix))
      colnames <- colnames(inMatrix)
      resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
      
      inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
      
      for(i in 1:matrixColSize) {
        for(j in 1:matrixColSize) { 
          resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                                 as.vector(inMatrix[,j]))
        }
      }
      colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
      as.dist(resultsMatrix)->resultsMatrix
      attr(resultsMatrix, "method") <- "dist"
      return(resultsMatrix) 
    }
    inputmatrix <- dist.JSD(data$abund)
  } else if(any(metric == c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis"))) {
    inputmatrix <- vegdist(t(data$abund), method = metric)
  }
  #################################### end of block ####################################
  
  #Generate data depending on the chosen ordination type
  if(type == "pca") {
    #make the model
    model <- rda(inputmatrix, ...)
    
    #axis (and data column) names
    x_axis_name <- paste0("PC", x_axis)
    y_axis_name <- paste0("PC", y_axis)
    
    #Calculate the amount of inertia explained by each axis
    totalvar <- round(model$CA$eig/model$tot.chi * 100, 1)
    
    #Calculate species- and site scores
    sitescores <- scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- scores(model, display = "species", choices = c(x_axis, y_axis))
    
  } else if(type == "rda") {
    if(is.null(constrain)) 
      stop("Argument constrain must be provided when performing constrained/canonical analysis.")
    #make the model
    codestring <- paste0("rda(inputmatrix~", paste(constrain, collapse = "+"), ", data$metadata, ...)") #function arguments written in the format "rda(x ~ y + z)" cannot be directly passed to rda(), now user just provides a vector
    model <-  eval(parse(text = codestring))
    
    #axes depend on the results
    x_axis_name <- paste0("RDA", x_axis)
    if (model$CCA$rank <= 1){
      y_axis_name <- "PC1"
      
      #Calculate the amount of inertia explained by each axis
      totalvar <- c(round(model$CCA$eig/model$tot.chi * 100, 1), #constrained of total space
                    round(model$CA$eig/model$tot.chi * 100, 1)  #UNconstrained of total space
      )
      constrainedvar <- c(round(model$CA$eig/model$CA$tot.chi * 100, 1)) #UNconstrained of total UNconstrained space
    } else if (model$CCA$rank > 1) {
      y_axis_name <- paste0("RDA", y_axis)
      
      #Calculate the amount of inertia explained by each axis
      totalvar <- c(round(model$CCA$eig/model$tot.chi * 100, 1), #constrained of total space
                    round(model$CA$eig/model$tot.chi * 100, 1)  #UNconstrained of total space
      )
      constrainedvar <- c(round(model$CCA$eig/model$CCA$tot.chi * 100, 1)) #constrained of total constrained space
    }
    
    #Calculate species- and site scores
    sitescores <- scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if(type == "nmds") {
    #make the model
    model <- metaMDS(inputmatrix, trace = FALSE, ...)
    
    #axis (and data column) names
    x_axis_name <- paste0("NMDS", x_axis)
    y_axis_name <- paste0("NMDS", y_axis)
    
    #Calculate species- and site scores
    #Speciesscores may not be available with MDS
    sitescores <- scores(model, display = "sites")
    if(!length(model$species) > 1) {
      speciesscores <- warning("Speciesscores are not available.")
    } else {
      speciesscores <- scores(model, display = "species", choices = c(x_axis, y_axis))
    }
  } else if(type == "mmds" | type == "pcoa") {
    #make the model
    model <- pcoa(inputmatrix, ...)
    
    #axis (and data column) names
    x_axis_name <- paste0("PCo", x_axis)
    y_axis_name <- paste0("PCo", y_axis)
    
    #Calculate the percentage of eigenvalues explained by the axes
    totalvar <- round(model$values$Relative_eig * 100, 1)
    names(totalvar) <- c(paste0("PCo", seq(1:length(totalvar))))
    
    #Calculate species- and site scores
    #Speciesscores are not available with pcoa
    sitescores <- as.data.frame(model$vectors)
    colnames(sitescores) <- c(paste0("PCo", seq(1:length(sitescores))))
    speciesscores <- warning("Speciesscores are not available.")
  } else if(type == "ca") {
    #make the model
    model <- cca(inputmatrix, ...)
    
    #axis (and data column) names
    x_axis_name <- paste0("CA", x_axis)
    y_axis_name <- paste0("CA", y_axis)
    
    #Calculate the percentage of eigenvalues explained by the axes
    totalvar <- round(model$CA$eig/model$tot.chi * 100, 1)
    
    #Calculate species- and site scores
    sitescores <- scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if(type == "cca") {
    if(is.null(constrain)) 
      stop("Argument constrain must be provided when performing constrained/canonical analysis.")
    #make the model
    codestring <- paste0("cca(inputmatrix~", paste(constrain, collapse = "+"), ", data$metadata, ...)") #function arguments written in the format "rda(x ~ y + z)" cannot be directly passed to rda(), now user just provides a vector
    model <-  eval(parse(text = codestring))
    
    #axes depend on the results
    x_axis_name <- paste0("CCA", x_axis)
    if (model$CCA$rank <= 1){
      y_axis_name <- "CA1"
      
      #Calculate the amount of inertia explained by each axis
      totalvar <- c(round(model$CCA$eig/model$tot.chi * 100, 1), #constrained of total space
                    round(model$CA$eig/model$tot.chi * 100, 1)  #UNconstrained of total space
      )
      constrainedvar <- c(round(model$CA$eig/model$CA$tot.chi * 100, 1)) #UNconstrained of total UNconstrained space
    } else if (model$CCA$rank > 1) {
      y_axis_name <- paste0("CCA", y_axis)
      
      #Calculate the amount of inertia explained by each axis
      totalvar <- c(round(model$CCA$eig/model$tot.chi * 100, 1), #constrained of total space
                    round(model$CA$eig/model$tot.chi * 100, 1)  #UNconstrained of total space
      )
      constrainedvar <- c(round(model$CCA$eig/model$CCA$tot.chi * 100, 1)) #constrained of total constrained space
    }
    
    #Calculate species- and site scores
    sitescores <- scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if(type == "dca") {
    #make the model
    model <- decorana(inputmatrix, ...)
    
    #axis (and data column) names
    x_axis_name <- paste0("DCA", x_axis)
    y_axis_name <- paste0("DCA", y_axis)
    
    #Calculate the percentage of eigenvalues explained by the axes
    #totalvar <- round(model$CA$eig/model$CA$tot.chi * 100, 1)
    
    #Calculate species- and site scores
    sitescores <- scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- scores(model, display = "species", choices = c(x_axis, y_axis))
  }
  #################################### end of block ####################################
  
  #Make data frames for ggplot
  dsites <- cbind.data.frame(data$metadata, sitescores)
  
  if (!is.null(sample_color_order)) {
    dsites[, sample_color] <- factor(dsites[, sample_color], levels = sample_color_order)
  }
  
  if(length(speciesscores) > 1) {
    dspecies <- merge(data.frame(speciesscores, OTU = rownames(speciesscores)), data$tax, by.x = "OTU")
    dspecies$dist <- dspecies[, x_axis_name]^2 + dspecies[, y_axis_name]^2
    dspecies <- arrange(dspecies, desc(dist))
    rownames(dspecies) <- dspecies$OTU
    if (species_rescale == TRUE) {
      maxx <- max(abs(dsites[, x_axis_name]))/max(abs(dspecies[,x_axis_name]))
      dspecies[, x_axis_name] <- dspecies[, x_axis_name] * maxx * 0.8
      maxy <- max(abs(dsites[, y_axis_name]))/max(abs(dspecies[,y_axis_name]))
      dspecies[, y_axis_name] <- dspecies[, y_axis_name] * maxy * 0.8
    }
  } else {
    dspecies = NULL
  }
  
  #Generate a nice ggplot with the coordinates from scores
  plot <- ggplot(dsites,
                 aes_string(x = x_axis_name,
                            y = y_axis_name,
                            color = sample_color,
                            shape = sample_shape))
  
  #Generate a color frame around the chosen color group
  if(sample_colorframe == TRUE) {
    if(is.null(sample_color)) stop("Please provide the argument sample_color")
    splitData <- split(dsites, dsites[, sample_color]) %>% 
      lapply(function(df) {
        df[chull(df[, x_axis_name], df[, y_axis_name]), ]
      })
    hulls <- do.call(rbind, splitData)
    plot <- plot + geom_polygon(data = hulls, aes_string(fill = sample_color, group = sample_color), alpha = 0.2*opacity)
  }
  
  # Add points and plotly functionality for samples
  if (!is.null(sample_plotly)){
    if(length(sample_plotly) > 1){
      data_plotly <- apply(data$metadata[,sample_plotly], 1, paste, collapse = "<br>")  
    } else if(sample_plotly == "all"){
      data_plotly <- apply(data$metadata[,], 1, paste, collapse = "<br>")  
    } else{
      data_plotly <- paste0(sample_plotly,": ",data$metadata[,sample_plotly])
    }
    plot <- plot +
      geom_point(size = 2, alpha = opacity,
                 aes(text = data_plotly)) + #HER
      theme_minimal() +
      theme(axis.line = element_line(colour = "black", size = 0.5))
  } else{
    plot <- plot +
      geom_point(size = 2, alpha = opacity) +
      theme_minimal() +
      theme(axis.line = element_line(colour = "black", size = 0.5))
  }
  
  #Only eigenvalue-based ordination methods can be displayed with % on axes
  if(type == "pca" | type == "ca" | type == "pcoa" | type == "mmds") {
    plot <- plot +
      xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "%]", sep = "")) + 
      ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "%]", sep = ""))
  } else if(type == "rda" | type == "cca") {
    if(model$CCA$rank > 1) {
      plot <- plot + 
        xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "% / ", constrainedvar[x_axis_name], "%]", sep = "")) +
        ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "% / ", constrainedvar[y_axis_name], "%]", sep = ""))
    } else if(model$CCA$rank <= 1) {
      plot <- plot + 
        xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "%]", sep = "")) +
        ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "% / ", constrainedvar[y_axis_name], "%]", sep = ""))
    }
  } else if (type == "nmds") {
    plot <- plot +
      annotate("text", size = 3, x = Inf, y = Inf, hjust = 1, vjust = 1, label = paste0("Stress value = ", round(model$stress, 3)))
  }
    
  #Plot species points
  if (species_plot == TRUE) {
    if(species_plotly == T){
      data_plotly <- paste("Kingdom: ", data$tax[,1],"<br>",
                           "Phylum: ", data$tax[,2],"<br>",
                           "Class: ", data$tax[,3],"<br>",
                           "Order: ", data$tax[,4],"<br>",
                           "Family: ", data$tax[,5],"<br>",
                           "Genus: ", data$tax[,6],"<br>",
                           "Species: ", data$tax[,7],"<br>",
                           "OTU: ", data$tax[,8],sep = "")
      plot <- plot + 
        geom_point(data = dspecies,
                   color = "darkgrey",
                   shape = species_shape,
                   size = species_size-1,
                   alpha = opacity,
                   aes(text = data_plotly))
    } else{
      plot <- plot + 
        geom_point(data = dspecies,
                   color = "darkgrey",
                   shape = species_shape,
                   size = species_size,
                   alpha = opacity) 
    }
    
  }
  
  #Plot text labels
  if (!is.null(sample_colorframe_label)) {
    temp <- data.frame(group = dsites[, sample_colorframe_label], 
                       x = dsites[, x_axis_name],
                       y = dsites[, y_axis_name]) %>% 
      group_by(group) %>%
      summarise(cx = mean(x), cy = mean(y)) %>% 
      as.data.frame()
      temp2 <- merge(dsites, temp,
                     by.x = sample_colorframe_label, 
                     by.y = "group")
      temp3 <- temp2[!duplicated(temp2[, sample_colorframe_label]), ]
      if (repel == T){plot <- plot +geom_text_repel(data = temp3, aes_string(x = "cx", y = "cy", label = sample_colorframe_label), size = 3,color = "black",fontface = 2)}
      else{plot <- plot +geom_text(data = temp3, aes_string(x = "cx", y = "cy", label = sample_colorframe_label), size = 3,color = "black",fontface = 2)}
  }
  
  #sample_trajectory
  if (!is.null(sample_trajectory)) {
    traj <- dsites[order(dsites[, sample_trajectory]), ]
    plot <- plot + geom_path(data = traj, aes_string(group = sample_trajectory_group))
  }
  
  #Sample point labels
  if(!is.null(sample_label)) {
    if (repel == T){plot <- plot + geom_text_repel(aes_string(label = sample_label),size = sample_label_size, color = "grey40", segment.color = sample_label_segment_color)}
    else{plot <- plot + geom_text(aes_string(label = sample_label),size = sample_label_size, color = "grey40", segment.color = sample_label_segment_color)}
  }
  
  #Plot species labels
  if (species_nlabels > 0) {
    if (repel == T){plot <- plot +geom_text_repel(data = dspecies[1:species_nlabels,], aes_string(x = x_axis_name, y = y_axis_name, label = species_label_taxonomy),colour = species_label_color, size = species_label_size,fontface = 4,inherit.aes = FALSE)}
    else{plot <- plot +geom_text(data = dspecies[1:species_nlabels,], aes_string(x = x_axis_name, y = y_axis_name, label = species_label_taxonomy),colour = species_label_color, size = species_label_size,fontface = 4,inherit.aes = FALSE)}
  }
  
  ######## Fit environmental variables ########
  # Categorial fitting
  if(!is.null(envfit_factor)) {
    evf_factor_model <- envfit(model,
                               data$metadata[,envfit_factor, drop = FALSE],
                               permutations = 999,
                               choices = c(x_axis_name, y_axis_name)
    )
    evf_factor_data <- data.frame(Name = rownames(evf_factor_model$factors$centroids),
                                  Variable = evf_factor_model$factors$var.id,
                                  evf_factor_model$factors$centroids,
                                  pval = evf_factor_model$factors$pvals
    ) %>% subset(pval <= envfit_signif_level)
    if (nrow(evf_factor_data) > 0 & envfit_show == TRUE) {
      if (repel == T){plot <- plot + geom_text_repel(data = evf_factor_data,aes_string(x = x_axis_name, y = y_axis_name, label = "Name"), colour = envfit_color, inherit.aes = FALSE, size = envfit_textsize, fontface = "bold")}
      else{plot <- plot + geom_text(data = evf_factor_data,aes_string(x = x_axis_name, y = y_axis_name, label = "Name"), colour = envfit_color, inherit.aes = FALSE, size = envfit_textsize, fontface = "bold")}
    }
    if (nrow(evf_factor_data) == 0) {
      warning("No environmental variables fit below the chosen significant level.")
    }
  } else {
    evf_factor_model <- NULL
  }
  
  # Numerical fitting
  if (!is.null(envfit_numeric)) {
    evf_numeric_model <- envfit(model,
                                data$metadata[,envfit_numeric, drop = FALSE],
                                permutations = 999,
                                choices = c(x_axis_name, y_axis_name)
    )
    evf_numeric_data <- data.frame(Name = rownames(evf_numeric_model$vectors$arrows),
                                   evf_numeric_model$vectors$arrows * sqrt(evf_numeric_model$vectors$r) * envfit_numeric_arrows_scale,
                                   pval = evf_numeric_model$vectors$pvals
    ) %>% subset(pval <= envfit_signif_level)
    if (nrow(evf_numeric_data) > 0 & envfit_show == TRUE) {
      plot <- plot + geom_segment(data = evf_numeric_data,
                                  aes_string(x = 0,
                                             xend = x_axis_name,
                                             y = 0,
                                             yend = y_axis_name
                                  ),
                                  arrow = arrow(length = unit(3, "mm")),
                                  colour = "darkred",
                                  size = 1,
                                  inherit.aes = FALSE) + 
        geom_text(data = evf_numeric_data,
                  aes_string(x = x_axis_name,
                             y = y_axis_name,
                             label = "Name"),
                  colour = envfit_color,
                  inherit.aes = FALSE,
                  size = envfit_textsize,
                  hjust = 1.2,
                  vjust = 1.2,
                  fontface = "bold"
        )
    } 
    if (nrow(evf_numeric_data) == 0) {
      warning("No environmental variables fit below the chosen significant level.")
    }
  } else {
    evf_numeric_model <- NULL
  }
  
  #################################### end of block ####################################
  
  #return plot or additional details
  if(!is.null(sample_plotly)){
    ggplotly(plot, tooltip = "text") %>% 
      layout(showlegend = FALSE)
  } 
  else if(species_plotly == T){
    ggplotly(plot, tooltip = "text") %>% 
      layout(showlegend = FALSE)
    }
  else if(output == "plot"){
    return(plot)
  }
  else if(output == "detailed"){
    if (type == "nmds") {
      screeplot <- NULL
    } else {
      ### screeplot ###
      #the data for it
      if (type == "mmds" | type == "pcoa") {
        if (length(model$values$Relative_eig) > 10) {
          unconstrained_eig <- model$values$Relative_eig[1:10]*100
        } else {
          unconstrained_eig <- model$values$Relative_eig*100
        }
        #the scree plot
        screeplot <- ggplot(data.frame(axis = factor(as.character(c(1:length(unconstrained_eig))), levels = c(1:length(unconstrained_eig))), eigenvalues = unconstrained_eig), aes(x = axis, y = eigenvalues)) +
          geom_col() +
          #geom_text(label = round(eigenvalues, 2), vjust = -1, size = 3)  + #Can't get it to work
          theme_minimal() +
          xlab("Axis (max. 10 axes will be shown)") +
          ylab("Eigenvalue in percent of total inertia")
      } else {
        unconstrained_eig <- model$CA$eig/model$tot.chi*100
        constrained_eig <- model$CCA$eig/model$tot.chi*100
        if (length(constrained_eig) > 10) {
          constrained_eig <- constrained_eig[1:10]
        }
        if (length(unconstrained_eig) > 10) {
          unconstrained_eig <- unconstrained_eig[1:10]
        }
        eigenvalues <- c(constrained_eig, unconstrained_eig) #constrained combined with unconstrained
        #the scree plot
        screeplot <- ggplot(data.frame(axis = factor(names(eigenvalues), levels = names(eigenvalues)), eigenvalues = eigenvalues), aes(x = axis, y = eigenvalues)) +
          geom_col() +
          geom_text(label = round(eigenvalues, 2), vjust = -1, size = 3)  +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          xlab("Axis (max. 10 axes will be shown)") +
          ylab("Eigenvalue in percent of total inertia")
      }
    }
  return(list(plot = plot,
              screeplot = screeplot,
              model = model,
              dsites = dsites,
              dspecies = dspecies,
              evf_factor_model = evf_factor_model,
              evf_numeric_model = evf_numeric_model)
         ) 
  }
}
