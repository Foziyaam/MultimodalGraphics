# Load Libraries
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(reshape2)
library(lattice)
library(plyr)
library(scales)
library(directlabels)

# Define the S4 class as ScatteredPlot
setClass("ScatteredPlot",
         representation(
           metadata_filepath = "character",
           metadata = "data.frame",
           plot = "ANY"
         )
)

# Constructor method for the class
setMethod("initialize",
          signature(.Object = "ScatteredPlot"),
          function(.Object, metadata_filepath, timePointLevels = c("TP-1", "TP-2", "TP-3", "TP-4"), pValueColumn = "p", qValueColumn = "q", expressionColumnName = "log2fc",
                   highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue =  1.301) {
            #pValueColumn = pValueColumn
            .Object@metadata_filepath <- metadata_filepath
            .Object@metadata <- process_metadata(metadata_filepath, pValueColumn = pValueColumn, qValueColumn = qValueColumn, expressionColumnName = expressionColumnName,
                                                 highLog2fc = highLog2fc, lowLog2fc = lowLog2fc, negLog10pValue =  negLog10pValue)
            .Object@plot <- NULL
            return(.Object)
          }
)

# Define a method to create the plot
setGeneric("createPlot", function(object, ...) standardGeneric("createPlot"))

setMethod("createPlot",
          signature(object = "ScatteredPlot"),
          function(object, color1 = "cornflowerblue", color2 = "grey", color3="indianred", 
                   highLog2fc = 0.585, lowLog2fc = -0.585, expressionDirection = "regulation",
                   timeVariable="reg_time_org") {
            if (is.null(object@plot)) {
              object@plot <- create_plot(object@metadata, color1 = color1, color2 = color2, color3=color3,
                                         highLog2fc = highLog2fc, lowLog2fc = lowLog2fc,
                                         expressionDirection = expressionDirection,
                                         timeVariable = timeVariable)
            }
            return(object)
          }
)

# Define a method to print the plot
setMethod("show",
          signature(object = "ScatteredPlot"),
          function(object) {
            object <- createPlot(object)
            print(object@plot)
          }
)

# Metadata Processing Function
process_metadata <- function(filepath, timePointLevels, pValueColumn = "p", qValueColumn = "q", expressionColumnName = "log2fc",
                             highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue =  1.301) {
  metadata <- read.csv(filepath)
  
  metadata$neglog10p = -log10(metadata[[pValueColumn]])
  metadata$neglog10q = -log10(metadata[[qValueColumn]])
  
  #this for p < 0.05 and FC > 1.5 (linear fold change)
  metadata$color_flag <-
    ifelse(
      metadata$log2fc > highLog2fc &
        metadata$neglog10p > negLog10pValue, 1, ifelse(metadata$log2fc < lowLog2fc &
                                                metadata$neglog10p > negLog10pValue,-1, 0)
    )
  
  metadata$timePoint <- factor(
    metadata$timePoint, levels = timePointLevels
  )
  # View(metadata[,-c(14,15,17)])
  # write.csv(metadata[,-c(14,15,17)],"scattered_data.csv",row.names = F)
  return(metadata)
}


# Plotting Function
create_plot <- function(metadata, color1 = "cornflowerblue", color2 = "grey", color3="indianred",
                        highLog2fc = 0.585, lowLog2fc = -0.585, expressionDirection = "regulation",
                        timeVariable="reg_time_org") {
  colorFlag = "color_flag"
  gp_obj <- ggplot(data = metadata, aes(
    x = log2fc, y = neglog10p, color = as.factor(metadata[[colorFlag]])
  )) +
    geom_point(alpha = 0.5, size = 1.75) +
    theme(legend.position = "none") + geom_jitter() +
    scale_color_manual(values = c(color1, color2, color3)) +
    labs(
      x = expression("log2 (fold change)"), y = expression("-log10 (p-value)")
    ) + facet_grid(organ ~ timePoint, space = "free") + theme_bw() +
    theme(
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12)
    )
  
  gp_obj2 = gp_obj + theme(
    axis.text = element_text(size = 12),axis.title.x = element_text(size = 12, face =
                                                                      "bold"),title = element_text(size = 12, face = "bold"), strip.text.x = element_text(
                                                                        size = 12, face = "bold",colour = "black", angle = 0
                                                                      )
  )
  gp_obj3 = gp_obj2 + theme(
    axis.text = element_text(size = 12),axis.title.y = element_text(size = 12, face =
                                                                      "bold"),title = element_text(size = 12, face = "bold"), strip.text.y = element_text(
                                                                        size = 12, face = "bold", colour = "black", angle = 90
                                                                      )
  )
  p <- gp_obj3 + theme(legend.position = "right") 
  
  p1 <-
    p + theme(strip.background = element_rect(
      fill = "white", color = color2, linewidth = 1
    ))
  
  p2 <- p1 + theme(axis.title.x = element_text(size =
                                                 12, face = "bold")) + theme(axis.title.x = element_text(size = 12, face =
                                                                                                           "bold")) + scale_color_manual(
                                                                                                             name = "log2 (fold change)", 
                                                                                                             values = c(color1, color2, color3), 
                                                                                                             labels = c(paste("<", lowLog2fc), paste(lowLog2fc,"to", highLog2fc), paste(">", highLog2fc))
                                                                                                           ) + scale_fill_hue()
  
  metadata2 <- metadata[metadata[[colorFlag]] == 1, ]
  metadata3 <- metadata[metadata[[colorFlag]] == -1, ]
  metadata4 <- rbind(metadata2, metadata3)
  
  metadata5 <- metadata4[metadata4[[expressionDirection]] == "up", ]
  metadata6 <- metadata4[metadata4[[expressionDirection]] == "down", ]
  
  meta.cor1 <- ddply(.data = metadata5,
                     .(organ, timePoint),
                     .fun = function(x) {
                       summarize(x, n1 = paste(length(x[[timeVariable]])))
                     }
  )
  
  p3 <- p2 + geom_text(data = meta.cor1, aes(x = 4, y = 9, label = n1),
                       colour = color3, inherit.aes = FALSE, parse = TRUE
  )
  
  meta.cor2 <- ddply(.data = metadata6,
                     .(organ, timePoint),
                     .fun = function(x) {
                       summarize(x, n2 = paste(length(x[[timeVariable]])))
                     }
                     
  )
  
  p4 <- p3 + geom_text(data = meta.cor2, aes(x = -4, y = 9, label = n2),
                       colour = color1, inherit.aes = FALSE, parse = FALSE
  )
  return(p4 + theme(legend.position = "bottom"))
}
