#Script Meta-analysis 
#By: Anny S. Adri - PhD Student - University of São Paulo (USP) 
#References: 10.18129/B9.bioc.MetaVolcanoR
#References: https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

# Load required libraries, installing them if necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github('kevinblighe/EnhancedVolcano')
BiocManager::install("MetaVolcanoR", dependencies = TRUE)
BiocManager::install('EnhancedVolcano')

# Load libraries
library(MetaVolcanoR)
library(ggplot2)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(openxlsx)
library(EnhancedVolcano)

# Create a list of dataframes for the meta-analysis
diffexplist <- list(
  "counts_file_1" = counts_file_1,
  "counts_file_2" = counts_file_2,
  "counts_file_3" = counts_file_3
)

# Run MetaVolcanoPlot and generate results
meta_degs_comb <- combining_mv(
  diffexp = diffexplist,
  pcriteria = 'pvalue',
  foldchangecol = 'log2FoldChange',
  genenamecol = 'Genes',
  geneidcol = NULL,
  metafc = 'Mean',
  metathr = 0.01,
  collaps = TRUE,
  jobname = "MetaVolcano",
  outputfolder = ".",
  draw = 'HTML'
)

# Convert results to dataframe and save
metaDegs <- meta_degs_comb@metaresult
openxlsx::write.xlsx(metaDegs, "metaDegs.xlsx")

# Remove duplicates and filter significant results
datawf <- metaDegs[!duplicated(metaDegs$Genes),]
sig <- datawf[datawf$metap < 0.05,]
up <- subset(sig, metafc > 0)
down <- subset(sig, metafc < 0)

# Save filtered tables
write.xlsx(sig, "sig_notconverted.xlsx")
write.xlsx(up, "up_notconverted.xlsx")
write.xlsx(down, "down_notconverted.xlsx")

#---------------Create Metavolcano plot-----------------------#

# Define colors based on significance
keyvals <- ifelse(
  metaDegs$metafc > 0 & metaDegs$metap < 0.05, 'brown2',
  ifelse(metaDegs$metafc < 0 & metaDegs$metap < 0.05, 'deepskyblue3', 'gray')
)
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'brown2'] <- 'up-regulated'
names(keyvals)[keyvals == 'gray'] <- 'nonDEGs'
names(keyvals)[keyvals == 'deepskyblue3'] <- 'down-regulated'

# Summary of key values
table(keyvals)

# Define plot limits
xmin <- sort(na.exclude(metaDegs$metafc))[1] - 0.5
h <- 1 + (-log10(min(na.omit(metaDegs$metap))))
xmax <- sort(na.omit(metaDegs$metafc), decreasing = TRUE)[2] + 0.5

# Identify top differentially expressed genes
degs <- which(names(keyvals) %in% c('up-regulated', 'down-regulated'))
topdegs <- c(
  head(order(metaDegs$metafc[degs]), n = 10),
  head(order(metaDegs$metafc[degs], decreasing = TRUE), n = 10),
  head(order(metaDegs$metap[degs]), n = 10)
)

# Plot using EnhancedVolcano
p <- EnhancedVolcano(metaDegs,
                     lab = metaDegs$Genes,
                     x = "metafc",
                     y = "metap",
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     xlim = c(xmin, xmax),
                     ylim = c(-0.05, h),
                     xlab = bquote(~metafc ~ "fold change"),
                     ylab = bquote(-log10() ~ italic(p)*"adj"),
                     pointSize = 1,
                     selectLab = metaDegs$Genes[degs][topdegs],
                     labSize = 4,
                     cutoffLineType = 'twodash',
                     legendPosition = "right",
                     legendLabSize = 13,
                     colAlpha = 0.9,
                     labCol = 'grey21',
                     border = 'full',
                     borderColour = 'grey41',
                     borderWidth = 0.6,
                     drawConnectors = FALSE,
                     widthConnectors = 0.5,
                     typeConnectors = "open",
                     lengthConnectors = unit(0.01, "npc"),
                     arrowheads = FALSE,
                     cutoffLineCol = 'grey31',
                     cutoffLineWidth = 0.4,
                     colCustom = keyvals,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE
) + 
  ggplot2::annotate("label", 
                    x = c(xmin + 2, 0, xmax - 2), 
                    y = h - 0.5, 
                    label = table(keyvals)[c(2, 3, 1)], 
                    col = c("deepskyblue3", "grey21", "brown2"))

# Save the plot
ggsave("metavolcano.tiff", plot = p, device = "tiff", width = 7, height = 5)

