#Differential Expression Analysis
#By: Anny S. Adri 12/06/2023 - PhD Student - University of São Paulo (USP)
#Reference: Michael I. Love, Simon Anders, and Wolfgang Huber (04/30/2024)

# Load required libraries
library(readr)
library(dplyr)
library(GEOquery)
library(tidyverse)
library(DESeq2)
library(openxlsx)
library(pheatmap)

# Generalized function to load and process any dataset
process_dataset <- function(count_file, metadata_file, geo_accession = NULL, output_prefix, design_formula = ~ 1) {
  
  # Import count data and metadata
  counts_data <- read_delim(count_file, delim = ";", escape_double = FALSE, trim_ws = TRUE)
  metadata <- openxlsx::read.xlsx(metadata_file, colNames = TRUE, rowNames = TRUE)
  
  # Optionally import GEO data if a GEO accession is provided
  if (!is.null(geo_accession)) {
    gse <- getGEO(geo_accession)[[1]]
    sample_info <- pData(gse)
    metadata <- cbind(metadata, sample_info)
  }
  
  # Ensure column names in count data match row names in metadata
  if (!all(colnames(counts_data) %in% rownames(metadata))) {
    stop("Column names in count data do not match row names in metadata.")
  }
  
  # Construct DESeqDataSet using a flexible design formula
  dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                                colData = metadata, 
                                design = design_formula)
  
  # Pre-filter: remove genes with low counts across samples
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Perform variance stabilizing transformation and plot PCA
  vsd <- vst(dds)
  DESeq2::plotPCA(vsd, intgroup = colnames(metadata)[1])  # Plot PCA with the first metadata column (or adjust as needed)
  
  # Relevel the first level of the primary group (if specified)
  if ("Group" %in% colnames(metadata)) {
    dds$Group <- relevel(dds$Group, ref = "control")
  }
  
  # Run DESeq
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Differential expression with flexible contrast
  resultsNames <- resultsNames(dds)
  Contrast.DEG.results <- results(dds, 
                                  independentFiltering = TRUE, 
                                  alpha = 0.05)
  
  # Export results
  res_sig <- subset(Contrast.DEG.results, padj < 0.05)
  resl2fc_up <- subset(res_sig, log2FoldChange > 0)
  resl2fc_down <- subset(res_sig, log2FoldChange < 0)
  
  write.xlsx(res_sig, paste0(output_prefix, "_significant_transcripts.xlsx"), rowNames = TRUE)
  write.xlsx(resl2fc_up, paste0(output_prefix, "_transcripts_up.xlsx"), rowNames = TRUE)
  write.xlsx(resl2fc_down, paste0(output_prefix, "_transcripts_down.xlsx"), rowNames = TRUE)
  
  
# Example call of the function for a dataset
process_dataset("count_file.csv", "metadata_file.xlsx", geo_accession = "GSEXXXXX", output_prefix = "output", design_formula = ~ Group + Gender)
process_dataset("count_file.csv", "metadata_file.xlsx", geo_accession = "GSEXXXXX", output_prefix = "output", design_formula = ~ Group + Sex)



