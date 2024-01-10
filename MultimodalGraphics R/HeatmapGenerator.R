# Load required packages
library(ComplexHeatmap)
library(seriation)
library(circlize)
library(gridtext)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(GetoptLong)
library(methods)
library(reshape)
library(dplyr)

# Define the S4 Class
setClass("HeatmapGenerator",
         slots = list(
           data_matrix = "matrix",
           small_mat = "ANY",
           small_mat_pv = "ANY",
           selected_columns = "numeric",
           selected_columns_pv = "numeric",
           heatmap_object = "ANY"
         )
)

# Constructor Method
setGeneric("newHeatmapGenerator", function(data_matrix, selected_columns, selected_columns_pv) {
  standardGeneric("newHeatmapGenerator")
})

setMethod("newHeatmapGenerator", signature(data_matrix = "matrix", selected_columns = "numeric", selected_columns_pv = "numeric"),
          function(data_matrix, selected_columns, selected_columns_pv) {
            new("HeatmapGenerator",
                data_matrix = data_matrix,
                selected_columns = selected_columns,
                selected_columns_pv = selected_columns_pv,
                small_mat = NULL,
                small_mat_pv = NULL,
                heatmap_object = NULL)
          }
)


# Method for Creating Row Annotations
setGeneric("createRowAnnotation", function(object, ...) {
  standardGeneric("createRowAnnotation")
})

setMethod("createRowAnnotation", signature(object = "HeatmapGenerator"),
          function(object, data, textColumn, valueColumn, pValueColumn, pValueThreshold1, 
                   textJustification = "right", lUnit, tValue, fValue, order = "textFirst", checkNull = FALSE,
                   fold_col_fun = colorRamp2(c(-2, 0, 2), c("cyan", "white", "orange"))) {
            
            
            # Check and use the external data instead of object@data_matrix
            if (is.null(data) || ncol(data) == 0) {
              stop("Data not provided or is empty")
            }
            
            if (!is.null(pValueColumn)) {
              if (checkNull) {
                pch_values <- ifelse(!is.na(data[[pValueColumn]]) & data[[pValueColumn]] < pValueThreshold1, tValue, fValue)
              } else {
                pch_values <- ifelse(data[[pValueColumn]] < pValueThreshold1, tValue, fValue)
              }
            }
            
            ha_annotation <- if (order == "textFirst") {
              rowAnnotation(
                text = anno_text(data[[textColumn]], gp = gpar(fontsize = 10, fontface = "bold"), just = textJustification, location = unit(lUnit, "npc")),
                .= anno_simple(data[[valueColumn]], col = fold_col_fun, na_col = "white", pch = pch_values)
              )
            } else if (order == "valueFirst") {
              rowAnnotation(
                .= anno_simple(data[[valueColumn]], col = fold_col_fun, na_col = "white", pch = pch_values),
                text = anno_text(data[[textColumn]], gp = gpar(fontsize = 10, fontface = "bold"), just = textJustification, location = unit(lUnit, "npc"))
              )
            } else if (order == "textOnly") {
              rowAnnotation(
                text = anno_text(data[[textColumn]], gp = gpar(fontsize = 8, fontface = "bold"), just = textJustification, location = unit(lUnit, "npc"))
              )
            } else if (order == "valueOnly") {
              rowAnnotation(
                .= anno_simple(data[[valueColumn]], col = fold_col_fun, na_col = "white", pch = pch_values)
              )
            }
            
            return(ha_annotation)
          }
)


# Method for Data Processing
setGeneric("calculate_small_mat_pv", function(object) {
  standardGeneric("calculate_small_mat_pv")
})

setMethod("calculate_small_mat_pv", "HeatmapGenerator", function(object) {
  if (length(object@selected_columns) < 2) {
    stop("At least two columns must be selected")
  }
  
  # Apply the seriate function
  o1 <- seriate(dist(object@data_matrix[, object@selected_columns]), method = "GW")
  o2 <- seriate(dist(t(object@data_matrix[, object@selected_columns])), method = "GW")
  
  # Create small_mat and small_mat_pv based on the seriated columns
  object@small_mat <- object@data_matrix[, object@selected_columns]
  object@small_mat_pv <- object@data_matrix[, object@selected_columns_pv]
  
  return(object)
})

# Method for Main Heatmap Generation
setGeneric("generate_heatmap", function(object, ...) {
  standardGeneric("generate_heatmap")
})

setMethod("generate_heatmap", "HeatmapGenerator", function(object, lowThreshold = 0.05, lowColor = "black", 
                                                           highThreshold = 0.1, highColor = "yellow", borderColor = "white",
                                                           col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red")), 
                                                           annotationSide = "none") {
  # Setting up heatmap arguments except for data
  heatmapArgs <- list(name = "log2FC", 
                      col = col_fun,
                      show_heatmap_legend = F,
                      column_km = 2, 
                      column_gap = unit(2, "mm"),
                      show_column_dend = F,
                      show_row_dend = F,
                      cluster_columns = F,
                      row_title = NULL,
                      column_title = NULL,
                      column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                      row_names_side = "left",
                      row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                      column_names_rot = 45,
                      row_km = 2,
                      rect_gp = gpar(col = borderColor, lwd = 2),
                      layer_fun = function(j, i, x, y, w, h, fill) {
                        ind_mat = restore_matrix(j, i, x, y)
                        for(ir in seq_len(nrow(ind_mat))) {
                          # start from the second column
                          for(ic in seq_len(ncol(ind_mat))) {
                            ind = ind_mat[ir, ic] # previous column
                            v = object@small_mat_pv[i[ind], j[ind]]
                            
                            grid.points(x[ind], y[ind], 
                                        pch = 16, gp = gpar(col = ifelse(v < lowThreshold, lowColor, ifelse(v>=lowThreshold && v<highThreshold, highColor, NA))), size = unit(1, "mm"))
                            
                          }
                        }
                      })
  
  # Adding the annotation parameter based on user choice
  if (annotationSide == "right") {
    heatmapArgs$right_annotation <- annotation
  } else if (annotationSide == "left") {
    heatmapArgs$left_annotation <- annotation
  }
  
  # Create the heatmap using do.call to pass the additional arguments
  object@heatmap_object <- do.call("Heatmap", c(list(object@small_mat), heatmapArgs))
  
  return(object)
})


# Method for Creating Row Annotations
setGeneric("addHeatmap", function(object, ...) {
  standardGeneric("addHeatmap")
})

setMethod("addHeatmap", signature(object = "HeatmapGenerator"),
          function(object, data, annotation, annotationSide = "right", borderColor = "white", wUnit = 10,
                   col_fun_dmr = colorRamp2(c(-0.015, 0, 0.01), c("green", "white", "magenta"))) {
            
            # Setting up heatmap arguments except for data
            heatmapArgs <- list(show_row_names = F, 
                                column_title = "",
                                column_title_side = "top",
                                column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                column_title_rot = 0,
                                show_heatmap_legend = F,
                                rect_gp = gpar(col = borderColor, lwd = 2),
                                cluster_rows = F,
                                cluster_columns = F, 
                                show_column_names = F,
                                col = col_fun_dmr, 
                                column_names_rot = 45,
                                width = unit(wUnit, "mm"))
            
            # Adding the annotation parameter based on user choice
            if (annotationSide == "right") {
              heatmapArgs$right_annotation <- annotation
            } else if (annotationSide == "left") {
              heatmapArgs$left_annotation <- annotation
            } else {
              stop("Invalid annotation side specified. Choose 'left' or 'right'.")
            }
            
            # Create the heatmap using do.call to pass the additional arguments
            otherHeatmap <- do.call("Heatmap", c(list(data), heatmapArgs))
            
            return(otherHeatmap)
          }
)

