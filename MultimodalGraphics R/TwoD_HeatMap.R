# Load required libraries
library(ggplot2)
library(RColorBrewer)

# Define the TwoD_HeatMap S4 class
setClass(
  "TwoD_HeatMap",
  slots = list(
    data = "data.frame",
    title = "character",
    x_label = "character",
    y_label = "character",
    color_palette = "character",
    breaks = "numeric"
  )
)

# Create a constructor function for the TwoD_HeatMap class
TwoD_HeatMap <- function(data, title = "Heatmap", x_label = "X Label", y_label = "Y Label",
                         color_palette = "RdYlBu", breaks = seq(-1, 1, 0.5)) {
  new("TwoD_HeatMap", data = data, title = title, x_label = x_label, y_label = y_label,
      color_palette = color_palette, breaks = breaks)
}

# Define a method for plotting the heatmap
setGeneric("plot_heatmap", function(object, ...) standardGeneric("plot_heatmap"))

setMethod("plot_heatmap", signature(object = "TwoD_HeatMap"), function(object, lowColor = "yellow", highColor = "red", borderColor="grey60", mySize = object@data$number_of_genes) {
  data <- object@data
  title <- object@title
  x_label <- object@x_label
  y_label <- object@y_label
  color_palette <- object@color_palette
  breaks <- object@breaks
  
  # Calculate neglog10p
  data$neglog10p <- -log10(data$p)
  
  # Define color palette
  color <- colorRampPalette(rev(brewer.pal(n = 7, name = color_palette)))(100)
  
  # Create the base ggplot object
  base_plot <- ggplot(data, aes(y = signaling, x = tissue)) +
    geom_tile(aes(fill = Activation_z_score), colour = borderColor) +
    scale_fill_gradientn(colours = color, breaks = breaks, labels = scales::comma_format()) +
    geom_point(aes(colour = neglog10p, size = mySize)) +
    scale_color_gradient(low = lowColor, high = highColor) +
    scale_size(range = c(1, 10)) +
    labs(
      x = x_label,
      y = y_label,
      title = title,
      fill = "Activation z-score",
      colour = "-log10(p-value)",
      size = "Number of Genes"
    ) +
    facet_grid(facets = . ~ timePoint, scales = "free_x", space = "free") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.border = element_rect(fill = NA, colour = "grey80", size = 0.6),
      axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 18, face = "bold"),
      title = element_text(size = 18),
      strip.text.x = element_text(size = 14, face = "bold", colour = "black", angle = 0)
    )
  
  # Display the plot
  print(base_plot)
})
