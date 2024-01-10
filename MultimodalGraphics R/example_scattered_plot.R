setwd("/cloud/project/rstudio")

# Load or source the ScatteredPlotClass.R file
source("ScatteredPlot.R")

# Create an instance of the ScatteredPlot class
metadata_filepath <- "scattered_data.csv"
scattered_plot <- new("ScatteredPlot", metadata_filepath = metadata_filepath,  pValueColumn = "p", 
                      qValueColumn = "q", expressionColumnName = "log2fc",
                      highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue =  1.301)

# Call methods on the object
scattered_plot <- createPlot(scattered_plot, color1 = "cornflowerblue", color2 = "grey", color3="indianred",
                             highLog2fc = 0.585, lowLog2fc = -0.585,
                             expressionDirection = "regulation",
                             timeVariable="reg_time_org")  # Create the plot

# Print the plot
scattered_plot  # This will call the 'show' method and display the plot

