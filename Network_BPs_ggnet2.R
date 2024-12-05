# Biological Process Network Analysis Script
# Author: Adriel Nobile and Anny Adri - PhD Student, University of Sao Paulo (USP)
# Description: This script performs GO term network analysis on biological processes with p < 0.05.
# References: https://briatte.github.io/ggnet/

# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

# Load packages
library(readxl)
library(dplyr)
library(tidyr)
library(GGally)
library(ggnet)
library(network)
library(sna)
library(ggplot2)
library(RColorBrewer)
library(intergraph)
library(svglite)
library(openxlsx)

# Set working directory (adjust path as needed)
setwd("/path/to/your/working/directory")

# Load input files (adjust file names as needed)
down <- read_excel("down_enrich.xlsx", sheet = 2)
up <- read_excel("up_enrich.xlsx", sheet = 2)

# Add columns indicating 'up-regulated' or 'down-regulated'
down$DEGs <- "down-regulated"
up$DEGs <- "up-regulated"

# Rename columns for consistency
colnames(down_enrich_1)[2] <- "Systems"
both <- rbind(down_enrich_1, up_enrich_1)

# Remove text within parentheses in the 'Term' column
both <- both %>%
  mutate(Term = gsub("\\s*\\([^\\)]+\\)", "", Term))

# Split 'Genes' column entries by ';' into separate rows
both_separated <- both %>%
  separate_rows(Genes, sep = ";")

# Filter out rows where Systems is "Others"
both_separated <- both_separated %>% filter(Systems != "Others")

# Prepare network data
dfnet <- both_separated[c(1, 2, 10, 11)] # Selecting relevant columns
dfnet_total <- dfnet %>% distinct(Term, Genes, Systems, DEGs, .keep_all = TRUE)

# Create network object
net <- network(dfnet_total[, c("Term", "Genes")])

# Prepare unique data frame for network attributes
u <- data.frame(Term = network.vertex.names(net))
bp_df <- dfnet[c(1, 2, 4)]
gene_df <- dfnet[c(3, 2, 4)]
colnames(gene_df)[1] <- colnames(bp_df)[1]
u <- unique(left_join(u, rbind(bp_df, gene_df)))

# Create interaction matrix for system type classification
u %>%
  mutate(Immune = ifelse(Systems == "Immune", 1, 0),
         Nervous = ifelse(Systems == "Nervous", 1, 0)) %>%
  group_by(Term, DEGs) %>%
  summarise(Immune = sum(Immune), Nervous = sum(Nervous), .groups = 'drop') %>%
  mutate(condition = case_when(
    Immune == 1 & Nervous == 1 ~ "Neuroimmune",
    Immune == 1 & Nervous == 0 ~ "Immune",
    Immune == 0 & Nervous == 1 ~ "Nervous",
    .default = "Others"
  )) -> u

# Add condition and DEGs as network vertex attributes
net %v% "Groups" <- as.character(u$condition)
net %v% "DEGs" <- as.character(u$DEGs)

# Define processes of interest for network labeling
processos_selecionados <- c(
  "Cytokine-Mediated Signaling Pathway", "Vesicle-Mediated Transport In Synapse",
  "Synaptic Vesicle Recycling", "regulation of leukocyte chemotaxis",
  "T Cell Apoptotic Process", "positive regulation of leukocyte chemotaxis",
  "Response To Type I Interferon", "Regulation Of Lymphocyte Proliferation",
  "Regulation Of Cytokine Production", "Positive Regulation Of Leukocyte Chemotaxis",
  "Neural Tube Development", "Negative Regulation Of Chemokine Production",
  "Interleukin-1-Mediated Signaling Pathway", "Interleukin-27-Mediated Signaling Pathway",
  "Antiviral Innate Immune Response", "Vesicle-Mediated Transport In Synapse"
)

# Create labels for selected processes
net %v% "label" <- ifelse(network.vertex.names(net) %in% processos_selecionados |
                            !network.vertex.names(net) %in% unique(both$Term),
                          network.vertex.names(net), "")

# Plot the network with selected labels
p <- ggnet2(net,
            size.cut = 50,
            size = "degree",
            label.size = 3.3,
            shape = "Group",
            color = "Type",
            palette = c("Immune" = "#9BBEC8",
                        "Nervous" = "#435585",
                        "Neuroimmune" = "#C5705D",
                        "BPs" = "grey40"),
            edge.color = "gray80",
            legend.size = 5,
            alpha = 0.85,
            mode = "fruchtermanreingold",
            label = net %v% "label")
print(p)

# Save plot (adjust file name as needed)
ggsave(filename = "network_biological_processes.svg", plot = p, width = 11, height = 8)

# Save output tables
write.xlsx(u, "network_table_input_figure.xlsx")
write.xlsx(dfnet, "network_table_input_analysis.xlsx")
