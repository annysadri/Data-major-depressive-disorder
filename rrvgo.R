# Script for GO term similarity and reduction using rrvgo
# By: Anny S. Adri - PhD student University of Sao Paulo (USP)
# References: https://bioconductor.org/packages/release/bioc/html/rrvgo.html

# Load necessary libraries
if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
BiocManager::install("rrvgo")

library(rrvgo)  # Load rrvgo package
library(stringr)  # Load stringr package
library(dplyr)  # Load dplyr for data manipulation
library(readxl)  # Load readxl to read Excel files
library(openxlsx)  # Load openxlsx to save Excel files

# Load enrichment data
enrichment_file <- read_excel("enrichment_file.xlsx")

# Calculating the similarity matrix and reducing GO terms
simMatrix <- calculateSimMatrix(enrichment_file$ID,                                 
                                orgdb="org.Hs.eg.db", 
                                ont="BP", 
                                method="Rel")

# Setting up the scores
scores <- setNames(-log10(enrichment_file$qvalue), enrichment_file$ID)

# Reducing GO terms based on similarity
reducedTerms <- reduceSimMatrix(simMatrix,                                 
                                scores, 
                                threshold = 0.7, 
                                orgdb="org.Hs.eg.db")

# Save reduced GO terms and similarity matrix
openxlsx::write.xlsx(reducedTerms, "reducedTerms.xlsx")
matriz <- as.data.frame(simMatrix)
openxlsx::write.xlsx(matriz, "simMatrix.xlsx")

# Set up heatmap colors
colors_scale <- c("#0077b6", "white", "brown2")
heatmap_colors <- colorRampPalette(colors_scale)(100)

# Save heatmap as TIFF
tiff("rrvgo_full.tiff", width = 15, height = 10, units = "in", res = 200)

heatmapPlot(simMatrix, reducedTerms, 
            annotateParent=TRUE, 
            annotationLabel="parentTerm", 
            fontsize=10, 
            col=heatmap_colors)
dev.off()


# ------------------Selected Terms------------------#######
# Merge tables for selected terms heatmap

keys_selected <- c("alpha-beta T cell activation", "alpha-beta T cell differentiation", 
                   "B cell activation involved in immune response", 
                   "cellular response to catecholamine stimulus", 
                   "cellular response to dopamine", 
                   "cellular response to interleukin-2", 
                   "chemical synaptic transmission, postsynaptic", 
                   "embryonic limb morphogenesis", 
                   "excitatory postsynaptic potential", 
                   "immunoglobulin production involved in immunoglobulin-mediated immune response", 
                   "interleukin-2-mediated signaling pathway", 
                   "isotype switching", 
                   "limb morphogenesis", 
                   "natural killer cell activation", 
                   "negative regulation of interleukin-10 production", 
                   "neural tube development", 
                   "neurotransmitter secretion", 
                   "positive regulation of adaptive immune response", 
                   "positive regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains", 
                   "positive regulation of leukocyte cell-cell adhesion", 
                   "positive regulation of leukocyte differentiation", 
                   "positive regulation of leukocyte proliferation", 
                   "positive regulation of lymphocyte activation", 
                   "positive regulation of lymphocyte proliferation", 
                   "positive regulation of mononuclear cell proliferation", 
                   "positive regulation of T cell activation", 
                   "positive regulation of T cell proliferation", 
                   "regulation of alpha-beta T cell differentiation", 
                   "regulation of leukocyte proliferation", 
                   "regulation of lymphocyte proliferation", 
                   "regulation of mononuclear cell proliferation", 
                   "regulation of neurotransmitter secretion", 
                   "regulation of postsynaptic membrane potential", 
                   "regulation of T cell proliferation", 
                   "regulation of T-helper 1 type immune response", 
                   "response to catecholamine", "response to dopamine", 
                   "response to interleukin-2", "response to monoamine", 
                   "somatic diversification of immune receptors via germline recombination within a single locus", 
                   "somatic diversification of immunoglobulins involved in immune response", 
                   "somatic recombination of immunoglobulin gene segments", 
                   "somatic recombination of immunoglobulin genes involved in immune response", 
                   "synaptic vesicle cycle", "synaptic vesicle exocytosis", 
                   "T cell proliferation", "vesicle-mediated transport in synapse")

# Subset the data based on selected keys
df <- enrichment_file[enrichment_file$Description %in% keys_selected, ]

# Calculating the similarity matrix for selected terms
simMatrix_selected <- calculateSimMatrix(df$ID,                                 
                                         orgdb="org.Hs.eg.db", 
                                         ont="BP", 
                                         method="Rel")

# Set up scores for selected terms
scores_selected <- setNames(-log10(df$qvalue), df$ID)

# Reduce selected GO terms
reducedTerms_selected <- reduceSimMatrix(simMatrix_selected,                                 
                                         scores_selected, 
                                         threshold=0.7, 
                                         orgdb="org.Hs.eg.db")

# Save reduced selected GO terms and similarity matrix
tiff("rrvgo_selected.tiff", width = 10, height = 6, units = "in", res = 200)

heatmapPlot(simMatrix_selected, reducedTerms_selected, 
            annotateParent=TRUE, 
            annotationLabel="parentTerm", 
            fontsize=10, 
            col=heatmap_colors)
dev.off()

#save tables
openxlsx::write.xlsx(reducedTerms_selected, "reducedTerms_selected.xlsx")
openxlsx::write.xlsx(as.data.frame(simMatrix_selected), "simMatrix_selected.xlsx")
