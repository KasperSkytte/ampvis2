#' Venn diagram of core OTUs
#'
#' Calculates the number of "core" OTUs shared by groups given thresholds for how frequent the OTUs should be above a certain abundance. Also returns the average abundance of the OTUs in a particular group.
#'
#' @usage amp_venn(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param group_by Group the data based on a sample variable.
#' @param cut_a Abundance cutoff in percent. (\emph{default:} \code{0.1})
#' @param cut_f Frequency cutoff in percent. (\emph{default:} \code{80})
#' @param text_size Size of the plotted text. (\emph{default:} \code{5})
#' @param raw (\emph{logical}) Display raw input instead of converting to percentages. (\emph{default:} \code{FALSE}) 
#' @param detailed_output (\emph{logical}) Return additional details or not. If \code{TRUE}, it is recommended to save to an object and then access the additional data by \code{View(object$data)}. (\emph{default:} \code{FALSE})
#' 
#' @return A ggplot2 object.
#' 
#' @export
#' @import tidyr
#' @import ggplot2
#' @import dplyr
#' @examples 
#' #Load example data
#' data("AalborgWWTPs")
#' 
#' #Venn diagram grouped by WWTP
#' amp_venn(AalborgWWTPs, group_by = "Plant")
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_venn <- function(data, 
                     group_by = NULL,
                     cut_a = 0.1, 
                     cut_f = 80, 
                     text_size = 5, 
                     raw = FALSE,
                     detailed_output = FALSE){
  
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  OTU <- data[["tax"]]["OTU"]
  metadata <- data[["metadata"]]
  
  if (raw == F){
    abund <- as.data.frame(sapply(abund, function(x) x/sum(x)*100))
  }
  
  ## Test for number of groups
  if (length(levels(metadata[,group_by])) > 3){
    stop(paste("Only up to 3 different groups in the group_by variable is supported. The chosen group_by variable has", as.character(length(unique(metadata[,group_by]))), "different groups:\n",
               paste(unique(as.character(metadata[,group_by])), collapse = ", ")))
  }
  
  ## Select grouping variable
  
  colnames(metadata)[1] <- "SeqID"
  if (!is.null(group_by)){
    metadata <- metadata[,c("SeqID", group_by)]
    colnames(metadata)[2] <- "GRP"  
  } else {
    metadata <- data.frame(SeqID = metadata[,1], GRP = "Core")
  }
  
  ## Add OTU names to the abundance information
  abund1 <- cbind.data.frame(abund, OTU)
  
  ## Melt the dataframe for subsequent processing
  abund2 <- tidyr::gather(data = abund1, key = SeqID, value = Abundance, -OTU)
  
  ## Merge metadata information with the abundance data
  abund3 <- merge(abund2, metadata, by = "SeqID")
  
  ## Add frequent abundant column
  abund4 <- mutate(abund3, 
                   Freq = ifelse(Abundance > cut_a, 1, 0))
  
  ## Evaluate if OTUs are part of the core
  abund5 <- group_by(abund4, OTU, GRP) %>%
    summarise(Abundance = mean(Abundance),
              cFreq = sum(Freq)/n()*100,
              Core = ifelse(sum(Freq)/n()*100 >= cut_f, 1, 0))
  
  ## Convert back into matrix format
  a <- tidyr::spread(abund5[,c("OTU","GRP", "Core")], key = GRP, value = Core)
  
  ################### PLOT #############################
  
  ## 1 group
  if(ncol(a) == 2){
    c_A <- subset(a, a[,2] == 1)
    c_n <- subset(a, a[,2] == 0)
    a_A <- subset(abund5, OTU %in% c_A$OTU) %>% group_by(OTU) %>% 
           summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_n <- subset(abund5, OTU %in% c_n$OTU) %>% group_by(OTU) %>% 
           summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    
    AD <- data.frame(counts = c(A = nrow(c_A), N = nrow(c_n)), 
                     abund = round(c(A = a_A$Sum, N = a_n$Sum), 1))
    
    ## Plot  
    p <- ggplot(data.frame(), aes(x=0, y=0)) +
           annotate("text", x=c(0, 0), y = c(0,-0.5), 
                    label = c(paste(AD[1,1], "\n(", AD[1,2],if_else(raw == TRUE, "", "%"),")", sep = ""), 
                            paste("Non-core: ", AD[2,1]," (", AD[2,2],if_else(raw == TRUE, "", "%"),")", sep = "")), 
                    size = text_size) +
           annotate("text", x=c(0), y = 0.45, label = colnames(a)[2], size = text_size) +
           xlim(-0.65,0.65) +
           ylim(-0.65,0.65) +
           annotate("path", 
                    x=00.4*cos(seq(0,2*pi,length.out=100)),
                    y=0+0.4*sin(seq(0,2*pi,length.out=100))) +
           theme(panel.grid.major=element_blank(), 
                 panel.grid.minor=element_blank(), 
                 axis.text=element_blank(),
                 axis.title=element_blank(),
                 axis.ticks=element_blank(),
                 panel.border=element_blank(),
                 panel.background = element_blank(),
                 legend.key = element_blank(),
                 plot.margin = unit(c(0,0,0,0), "mm"))
    
    ## Generate lists of species in each group
    
    ot <- cbind.data.frame(tax, abund) %>%
           mutate(Shared = "Non-core") %>%
           mutate(Shared = ifelse(OTU %in% as.character(unique(c_A$OTU)), colnames(a)[2], Shared))
    
  res <- list(plot = p, 
                A = as.character(unique(c_A$OTU)),
                Noncore = as.character(unique(c_n$OTU)),
                Otutable = ot)
    
    names(res)[2] <- colnames(a)[2]
    
  }
  
  ## 2 groups
  if(ncol(a) == 3){
    c_AB <- subset(a, a[,2] == 1 & a[,3] == 1)
    c_A <- subset(a, a[,2] == 1 & a[,3] == 0)
    c_B <- subset(a, a[,2] == 0 & a[,3] == 1)
    c_n <- subset(a, a[,2] == 0 & a[,3] == 0)
    
    a_AB <- subset(abund5, OTU %in% c_AB$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_A <- subset(abund5, OTU %in% c_A$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_B <- subset(abund5, OTU %in% c_B$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_n <- subset(abund5, OTU %in% c_n$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    
    AD <- data.frame(counts = c(A = nrow(c_A), AB = nrow(c_AB), B = nrow(c_B), N = nrow(c_n)), 
                     abund = round(c(A = a_A$Sum, AB = a_AB$Sum, B = a_B$Sum, N = a_n$Sum), 1))
    
    p <- ggplot(data.frame(), aes(x=c(-0.2,0.2), y=0)) +
           annotate("text", x=c(-0.4, 0, 0.4, 0), y = c(0,0,0,-0.5), 
                    label = c(paste(AD[1:3,1], "\n(", AD[1:3,2],if_else(raw == TRUE, "", "%"),")", sep = ""), 
                              paste("Non-core: ", AD[4,1]," (", AD[4,2],if_else(raw == TRUE, "", "%"),")", sep = "")), 
                    size = text_size) +
           annotate("text", x=c(-0.2, 0.2), y = 0.45, label = colnames(a)[2:3], size = text_size) +
           xlim(-0.65,0.65) +
           ylim(-0.65,0.65) +
           annotate("path",
                     x=0.2+0.4*cos(seq(0,2*pi,length.out=100)),
                     y=0+0.4*sin(seq(0,2*pi,length.out=100))) +
           annotate("path",
                     x=-0.2+0.4*cos(seq(0,2*pi,length.out=100)),
                     y=0+0.4*sin(seq(0,2*pi,length.out=100))) +
           theme(panel.grid.major=element_blank(), 
                 panel.grid.minor=element_blank(), 
                 axis.text=element_blank(),
                 axis.title=element_blank(),
                 axis.ticks=element_blank(),
                 panel.border=element_blank(),
                 panel.background = element_blank(),
                 legend.key = element_blank(),
                 plot.margin = unit(c(0,0,0,0), "mm"))
    
    ot <- cbind.data.frame(tax, abund) %>%
      mutate(Shared = "Non-core") %>%
      mutate(Shared = ifelse(OTU %in% as.character(unique(c_A$OTU)), colnames(a)[2], Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_B$OTU)), colnames(a)[3], Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_AB$OTU)), "Core", Shared))
    
    res <- list(plot = p, 
                A = as.character(unique(c_A$OTU)),
                B = as.character(unique(c_B$OTU)),
                Core = as.character(unique(c_AB$OTU)),
                Noncore = as.character(unique(c_n$OTU)),
                Otutable = ot)
    
    names(res)[2:4] <- c(colnames(a)[2],
                         colnames(a)[3], 
                         "Core")
    
  }
  
  ## 3 groups
  if(ncol(a) == 4){
    c_ABC <- subset(a, a[,2] == 1 & a[,3] == 1 & a[,4] == 1)
    c_AB  <- subset(a, a[,2] == 1 & a[,3] == 1 & a[,4] == 0)
    c_AC  <- subset(a, a[,2] == 1 & a[,3] == 0 & a[,4] == 1)
    c_BC  <- subset(a, a[,2] == 0 & a[,3] == 1 & a[,4] == 1)
    c_A   <- subset(a, a[,2] == 1 & a[,3] == 0 & a[,4] == 0)
    c_B   <- subset(a, a[,2] == 0 & a[,3] == 1 & a[,4] == 0)
    c_C   <- subset(a, a[,2] == 0 & a[,3] == 0 & a[,4] == 1)
    c_n   <- subset(a, a[,2] == 0 & a[,3] == 0 & a[,4] == 0)
    
    
    a_ABC <- subset(abund5, OTU %in% c_ABC$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_AB  <- subset(abund5, OTU %in% c_AB$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_AC  <- subset(abund5, OTU %in% c_AC$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_BC  <- subset(abund5, OTU %in% c_BC$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_A   <- subset(abund5, OTU %in% c_A$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_B   <- subset(abund5, OTU %in% c_B$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_C   <- subset(abund5, OTU %in% c_C$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_n   <- subset(abund5, OTU %in% c_n$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    
    AD <- data.frame(counts = c(ABC = nrow(c_ABC), AB = nrow(c_AB),AC = nrow(c_AC), BC = nrow(c_BC), 
                                A = nrow(c_A), B = nrow(c_B), C = nrow(c_C), N = nrow(c_n)), 
                     abund = round(c(ABC = a_ABC$Sum, AB = a_AB$Sum, AC = a_AC$Sum, BC = a_BC$Sum, 
                                     A = a_A$Sum, B = a_B$Sum, C = a_C$Sum, N = a_n$Sum), 1))
    
    p <- ggplot(data.frame(), aes(x=c(-0.2,0.2), y=0)) +
      annotate("text", x=c(0, 0, -0.25, 0.25, -0.4, 0.4, 0, 0.5), y = c(0.05, 0.4, -0.05, -0.05, 0.3, 0.3, -0.3, -0.6), 
               label = c(paste(AD[1:7,1], "\n(", AD[1:7,2],if_else(raw == TRUE, "", "%"),")", sep = ""), 
                         paste("Non-core:\n", AD[8,1]," (", AD[8,2],if_else(raw == TRUE, "", "%"),")", sep = "")), 
               size = text_size) +
      annotate("text", x=c(-0.2, 0.2, 0), y = c(0.65, 0.65, -0.65), label = colnames(a)[2:4], size = text_size) +
      xlim(-0.65,0.65) +
      ylim(-0.65,0.65) +
      annotate("path",
               x=0.2+0.4*cos(seq(0,2*pi,length.out=100)),
               y=0.2+0.4*sin(seq(0,2*pi,length.out=100))) +
      annotate("path",
               x=-0.2+0.4*cos(seq(0,2*pi,length.out=100)),
               y=0.2+0.4*sin(seq(0,2*pi,length.out=100))) +
      annotate("path",
               x=0+0.4*cos(seq(0,2*pi,length.out=100)),
               y=-0.2+0.4*sin(seq(0,2*pi,length.out=100))) +  
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            axis.text=element_blank(),
            axis.title=element_blank(),
            axis.ticks=element_blank(),
            panel.border=element_blank(),
            panel.background = element_blank(),
            legend.key = element_blank(),
            plot.margin = unit(c(0,0,0,0), "mm"))
    
    ot <- cbind.data.frame(tax, abund) %>%
      mutate(Shared = "Non-core") %>%
      mutate(Shared = ifelse(OTU %in% as.character(unique(c_ABC$OTU)), "Core", Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_AB$OTU)), paste(colnames(a)[2], colnames(a)[3], sep = "_"), Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_AC$OTU)), paste(colnames(a)[2], colnames(a)[4], sep = "_"), Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_BC$OTU)), paste(colnames(a)[3], colnames(a)[4], sep = "_"), Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_A$OTU)), colnames(a)[2], Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_B$OTU)), colnames(a)[3], Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_C$OTU)), colnames(a)[4], Shared))
    
    res <- list(plot = p, 
                Core = as.character(unique(c_ABC$OTU)),
                AB = as.character(unique(c_AB$OTU)),
                AC = as.character(unique(c_AC$OTU)),
                BC = as.character(unique(c_BC$OTU)),
                A = as.character(unique(c_A$OTU)),
                B = as.character(unique(c_B$OTU)),
                C = as.character(unique(c_C$OTU)),
                Noncore = as.character(unique(c_n$OTU)),
                Otutable = ot)
    
    names(res)[2:8] <- c("Core",
                         paste(colnames(a)[2], colnames(a)[3], sep = "_"),
                         paste(colnames(a)[2], colnames(a)[4], sep = "_"),
                         paste(colnames(a)[3], colnames(a)[4], sep = "_"),
                         colnames(a)[2],
                         colnames(a)[3],
                         colnames(a)[4]
    )
    
  }  
  
  ## Export data
  if (detailed_output)
    return(res)
  if (!detailed_output)
    return(p)
}
