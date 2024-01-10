# Example usage:
# Create a TwoD_HeatMap object

setwd("/cloud/project/rstudio")
# Load or source the TwoD_HeatMap.R file
source("TwoD_HeatMap.R")

#beh <-read.csv("signaling3.csv",sep = ",", header = T)
beh <-read.csv("milk2.csv",sep = ",", header = T)
#View(beh)
heatmap_obj <- TwoD_HeatMap(data = beh, title = "Trauma Exposure and Post-Trauma Days",
                            x_label = "Brain Region", y_label = "Neuronal Signaling, Synaptic Plasticity and Neurogenesis",
                            color_palette = "RdYlBu", breaks = seq(-1, 1, 0.5))

# Plot the heatmap
plot_heatmap(heatmap_obj, lowColor = "green", highColor = "blue", borderColor="grey60", mySize = heatmap_obj@data$number_of_genes)

