# Synapse Network Analysis Script
# This script processes synaptic gene ontology (GO) terms and generates a network plot for the synapse enrichment analysis.
# Author: Adriel Nobile and Anny Adri - PhD Student, University of Sao Paulo (USP)
# References: https://briatte.github.io/ggnet/

# Load required libraries
library(dplyr)
library(tidyr)
library(network)
library(GGally)
library(ggplot2)
library(scales)
library(sna)
library(ggnet)
library(RColorBrewer)
library(readxl)
library(openxlsx)
library(stringr)

# Load data
data <- read.xlsx("data.xlsx")

# Create a copy of the data for processing
dfnet <- data

# Separate Genes in individual rows
dfnet <- dfnet %>% 
  separate_rows(Genes, sep = ";")

# Extract GO ID and create a new 'ID' column, removing it from 'Term'
dfnet <- dfnet %>% 
  mutate(ID = str_extract(Term, "\\(GO:\\d+\\)"),
         Term = str_replace(Term, " \\(GO:\\d+\\)", ""))

# Remove parentheses from GO ID column
dfnet$ID <- str_replace_all(dfnet$ID, "[()]", "")

# Filter out rows where System is 'Others'
dfnet <- dfnet %>% filter(System != "Others")

# Define "Neuroimmune" as a type based on Genes present in both "Immune" and "Nervous" systems
dfnet <- dfnet %>% 
  mutate(Type = case_when(
    Genes %in% dfnet$Genes[System == "Immune"] & Genes %in% dfnet$Genes[System == "Nervous"] ~ "Neuroimmune",
    TRUE ~ as.character(System)
  ))

# Create unique identifiers by concatenating 'Genes' and 'Type'
dfnet <- dfnet %>% 
  mutate(Genes = as.character(Genes), 
         Type = as.character(Type), 
         cgene = paste(Genes, Type, sep = " "))

# Create unique identifiers by concatenating 'Term' and 'System'
dfnet <- dfnet %>% 
  mutate(Term = as.character(Term), 
         System = as.character(System), 
         cterm = paste(Term, System, sep = " "))

# Remove duplicates based on concatenated 'Term' and 'Gene'
net_data <- dfnet %>% distinct(cterm, cgene, .keep_all = TRUE)

# Select relevant columns for the network data
net_data <- net_data[, c(2, 11, 13, 14, 15)]

# Rename columns for clarity
names(net_data)[1] <- "System"
names(net_data)[2] <- "Groups"
names(net_data)[4] <- "Genes"
names(net_data)[5] <- "Term"

# Create a network object using 'Term' and 'Genes' columns
net <- network(net_data[, c("Term", "Genes")], directed = FALSE)

# Create a list of vertex attributes (Type)
vertex_color <- unique(net_data[, c("Genes", "Type")])
vertex_color <- setNames(vertex_color$Type, vertex_color$Genes)

# Assign 'Type' attribute to network nodes
set.vertex.attribute(net, "Type", vertex_color[network.vertex.names(net)])

# Define colors for each 'Type'
df_color_dicio <- data.frame("Type" = c("Neuroimmune", "Immune", "Nervous", "BPs"), 
                             "color" = c("#C5705D", "#9BBEC8", "#435585", "gray50"))
color_palette <- setNames(df_color_dicio$color, df_color_dicio$Type)

# Assign colors based on 'Type'
df_color <- data.frame("Type" = get.vertex.attribute(net, "Type"))
df_color %>% 
  mutate(Type = ifelse(is.na(Type), "BPs", Type)) -> df_color
set.vertex.attribute(net, "Type", df_color$Type)

# Create a list of vertex shapes based on 'Groups'
vertex_shape <- unique(net_data[, c("Genes", "Groups")])
vertex_shape <- setNames(vertex_shape$Groups, vertex_shape$Genes)

# Assign 'Groups' attribute to network nodes
set.vertex.attribute(net, "Groups", vertex_shape[network.vertex.names(net)])

# Define shapes for each 'Group'
df_shape_dicio <- data.frame("Groups" = c("down", "up", "BPs"), 
                             "shape" = c(19, 17, 18))
shape_palette <- setNames(df_shape_dicio$shape, df_shape_dicio$Groups)

# Assign shapes based on 'Group'
df_shape <- data.frame("Groups" = get.vertex.attribute(net, "Groups"))
df_shape %>% 
  mutate(Groups = ifelse(is.na(Groups), "BPs", Groups)) -> df_shape
set.vertex.attribute(net, "Groups", df_shape$Groups)

# Create edge attributes based on 'Groups'
vertex_edge <- unique(net_data[, c("Genes", "Groups")])
vertex_edge <- setNames(vertex_edge$Groups, vertex_edge$Genes)
set.edge.attribute(net, "Groups", vertex_edge[network.vertex.names(net)])

# Define edge colors for each 'Group'
df_edge_dicio <- data.frame("Groups" = c("down", "up", "BPs"), 
                            "edges" = c("#5696c7", "#d74a49", "gray90"))
edge_palette <- setNames(df_edge_dicio$edges, df_edge_dicio$Groups)

# Assign edge colors
df_edge <- data.frame("Groups" = get.edge.attribute(net, "Groups"))
df_edge %>% 
  mutate(Groups = ifelse(is.na(Groups), "BPs", Groups)) -> df_edge
set.edge.attribute(net, "Groups", df_edge$Groups)

# Update edge colors based on 'Groups'
direction_colors <- ifelse(get.edge.attribute(net, "Groups") == "up", "#d74a49",
                           ifelse(get.edge.attribute(net, "Groups") == "down", "#5696c7",
                                  "gray90"))
set.edge.attribute(net, "Groups", direction_colors)

# Define biological processes to highlight
processos_selecionados <- c(
  "Cytokine-Mediated Signaling Pathway", "Vesicle-Mediated Transport In Synapse", 
  "Synaptic Vesicle Recycling", "Regulation of Leukocyte Chemotaxis", 
  "T Cell Apoptotic Process", "Positive Regulation of Leukocyte Chemotaxis", 
  "Response To Type I Interferon", "Regulation Of Lymphocyte Proliferation", 
  "Regulation Of Cytokine Production", "Positive Regulation Of Leukocyte Chemotaxis", 
  "Neural Tube Development", "Negative Regulation Of Chemokine Production", 
  "Interleukin-1-Mediated Signaling Pathway", "Interleukin-27-Mediated Signaling Pathway", 
  "Antiviral Innate Immune Response", "Vesicle-Mediated Transport In Synapse"
)

# Prepare node labels based on selected processes
teste <- as.data.frame(network.vertex.names(net))
teste$`network.vertex.names(net)` <- gsub(" Immune| Nervous| Neuroimmune", 
                                          "", teste$network.vertex.names)
names(teste)[1] <- "names"
network.vertex.names(net) <- teste$names

# Customize label sizes based on node type
nodes_names <- network.vertex.names(net)
label_size <- sapply(nodes_names, function(label) {
  if (label %in% processos_selecionados) { 4 }
  else if (label %in% genes) { 3 }
  else { 0 }
})

# Plot network
p <- ggnet2(net, 
            size = "degree", 
            shape = "Groups", 
            shape.palette = shape_palette, 
            edge.color = "gray90", 
            color = "Type", 
            palette = color_palette) + 
  theme(legend.position = "right") + 
  geom_text(aes(label = nodes_names), size = label_size)

# Save plot as SVG
ggsave(filename = "network_bps_pbmc.svg", plot = p, width = 11, height = 8)
