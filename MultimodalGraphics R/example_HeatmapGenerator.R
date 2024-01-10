setwd("/cloud/project/rstudio")

# Load or source the HeatmapGenerator.R file
source("HeatmapGenerator.R")

load_and_preprocess_data <- function(filepath) {
  # ... load and preprocess data ...
  combinedpathways <- read.csv(filepath, row.names = 1)
  oxidoinflammation <- combinedpathways[combinedpathways$pathway %in% c("inflammatory response", "oxidative stress"),]
  oxidoinflammation_matrix <- as.matrix(oxidoinflammation[,-30])
  return(oxidoinflammation_matrix)
}

# Read Main Heatmap data and Generate Heatmap
data_matrix <- load_and_preprocess_data("heatmap combined pathways proteins postmortem sbc marin_good.csv")
heatmap_gen <- newHeatmapGenerator(data_matrix, c(1,5,3,4,2,22,16,18,20), c(6,10,8,9,7,23,17,19,21))

heatmap_gen <- calculate_small_mat_pv(heatmap_gen)
heatmap_gen <- generate_heatmap(heatmap_gen, lowThreshold = 0.05, lowColor = "black", 
                                highThreshold = 0.1, highColor = "yellow", borderColor = "white",
                                col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red")), 
                                annotationSide = "none")

# Read Annotation data
data_path <- "GWAS-pQTL-miR pathways_v4_good_MVP.csv"
gwaspqtlmirpathways_mvp <- read.csv(data_path, row.names = 1)
gwaspqtlmirpathways = gwaspqtlmirpathways_mvp
gwaspqtlmiroxoinflammation = gwaspqtlmirpathways[gwaspqtlmirpathways$pathway=='inflammatory response'|gwaspqtlmirpathways$pathway=='oxidative stress',]

# Apply function to create annotations
ha_miR_oxoinflm <- createRowAnnotation(heatmap_gen, 
                                       gwaspqtlmiroxoinflammation,
                                       textColumn = "miR_ID", 
                                       valueColumn = "logFC_mir", 
                                       pValueColumn = "P.Value_mir", 
                                       pValueThreshold1 = 0.05, 
                                       textJustification = "left",
                                       lUnit = 0.05, 
                                       tValue = "*", 
                                       fValue = "", 
                                       order = "valueFirst", 
                                       checkNull = TRUE)

ha_gwas_mvp_ptsd_oxinfl <- createRowAnnotation(heatmap_gen, gwaspqtlmiroxoinflammation,
                                               textColumn = "SNP_ptsd_mvp", 
                                               valueColumn = "GWAS_log2FC", 
                                               pValueColumn = "P_GWAS_MVP_ptsd", 
                                               pValueThreshold1 = 0.0007, 
                                               textJustification = "right",
                                               lUnit = 1, 
                                               tValue = "***", 
                                               fValue = "**", 
                                               order = "textFirst", 
                                               checkNull = FALSE)

ha_gwas_pgc_ptsd_oxinfl <- createRowAnnotation(heatmap_gen, gwaspqtlmiroxoinflammation,
                                               textColumn = "SNP_ptsd_pgc", 
                                               valueColumn = "GWAS_log2FC", 
                                               pValueColumn = "P_GWAS_pgc_ptsd", 
                                               pValueThreshold1 = 0.0005, 
                                               textJustification = "right",
                                               lUnit = 1, 
                                               tValue = "***", 
                                               fValue = "**", 
                                               order = "textFirst", 
                                               checkNull = FALSE)

ha_pqtl_sbc_oxinfl <- createRowAnnotation(heatmap_gen, gwaspqtlmiroxoinflammation,
                                          textColumn = "SNP_pQTL_SBC", 
                                          valueColumn = "pQTL_log2FC", 
                                          pValueColumn = "P_pQQTL_SBC", 
                                          pValueThreshold1 = 0.0001, 
                                          textJustification = "right",
                                          lUnit = 1, tValue = "***", 
                                          fValue = "*", order = "textFirst", 
                                          checkNull = FALSE)

ha_pqtl_fcc_oxinfl <- createRowAnnotation(heatmap_gen, gwaspqtlmiroxoinflammation,
                                          textColumn = "SNP_pQTL_FCC", 
                                          valueColumn = "pQTL_log2FC", 
                                          pValueColumn = "P_pQQTL_FCC", 
                                          pValueThreshold1 = 0.0001, 
                                          textJustification = "right",
                                          lUnit = 1, tValue = "***", 
                                          fValue = "*", order = "textFirst", 
                                          checkNull = FALSE)

#
# --------------- heatmap for dmr-proteins -------------------------
#

dmrprotein = read.csv("DMR inflammation angiogenesis DMR-protein_v4_good.csv", row.names = 1)

# subset the DMR file for the inflammatory and oxidative stress pathways
dmroxoinflammation = dmrprotein[dmrprotein$pathway=="inflammatory response"|dmrprotein$pathway=="oxidative stress",]
dim(dmroxoinflammation)

dmrprotein_matrix = as.matrix(dmroxoinflammation[,2:4])

# change NA to zeros
dmrprotein_matrix[is.na(dmrprotein_matrix)] <- 0

# text annotation for the protein heatmap
ha_dmrprotein <- createRowAnnotation(heatmap_gen, dmroxoinflammation,
                                     textColumn = "chr_TSS1500_TSS200", 
                                     valueColumn = "", 
                                     pValueColumn = "", 
                                     pValueThreshold1 = 0, 
                                     textJustification = "left",
                                     lUnit = 0, tValue = "", 
                                     fValue = "", order = "textOnly", 
                                     checkNull = FALSE)

h2 <- addHeatmap(heatmap_gen, dmrprotein_matrix,ha_dmrprotein,"right")

#
#------------ protein and metabolites heatmap ------------------
#

# read metabolites data for the combined pathways
proteinmetabolite = read.csv("metabolites inflammation angiogenesis proteins-metabolites_v4_good.csv", row.names = 1)

# subset the metabolite data for oxidative stress and inflammatory response
metaboliteinflammationoxidative = proteinmetabolite[proteinmetabolite$pathway=="inflammatory response"|proteinmetabolite$pathway=="oxidative stress",]

# data frame to matrix
protmetabol_matrix = as.matrix(metaboliteinflammationoxidative[,2:6])

# change NA to zeros
protmetabol_matrix[is.na(protmetabol_matrix)] <- 0

# text annotation for the metabolite heatmap
ha_protmetabol <- createRowAnnotation(heatmap_gen, metaboliteinflammationoxidative,
                                      textColumn = "metabolites", 
                                      valueColumn = "", 
                                      pValueColumn = "", 
                                      pValueThreshold1 = 0, 
                                      textJustification = "right",
                                      lUnit = 1, tValue = "", 
                                      fValue = "", order = "textOnly", 
                                      checkNull = FALSE)


col_fun_metabol = colorRamp2(c(-0.5, 0, 0.4), c("dodgerblue1", "white", "deeppink1"))
h3 <- addHeatmap(heatmap_gen, protmetabol_matrix,ha_protmetabol, "left", "white", 15, col_fun_metabol)

# draw(heatmap_gen@heatmap_object)
### Main drawing
draw(h3 + heatmap_gen@heatmap_object + ha_miR_oxoinflm + h2 + ha_gwas_mvp_ptsd_oxinfl + ha_gwas_pgc_ptsd_oxinfl + 
       ha_pqtl_sbc_oxinfl + ha_pqtl_fcc_oxinfl, ht_gap = unit(c(3,3,3,3,3,3,3), "mm"), 
     main_heatmap = "log2FC", auto_adjust = FALSE)

