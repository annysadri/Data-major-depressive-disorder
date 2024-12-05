# GO Term Analysis Script
# Author: Débora Albuquerque, MSc student, University of Sao Paulo (USP)
# Author: Anny S. Adri, PhD Student, University of Sao Paulo (USP)
# Description: This script performs data preparation, GO term similarity analysis, and visualization.

# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Load packages
library(readr)
library(openxlsx)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggraph)
library(ggrepel)
library(forcats)
library(readxl)

# Set working directory (adjust the path as needed)
setwd("/Users/annysilvaadri/Documents/Transcriptoma certo/Meus peipers/Paper 1/Novas análises Synapse/Metanalise PBMC/Appyter")

# Load input files (adjust file names as needed)
down <- read_excel("down_enrich_tudo_input.xlsx")
up <- read_excel("up_enrich_tudo_input.xlsx")

# Combine data
combined_data <- rbind(down, up)

# Load additional base file (adjust file name as needed)
base_data <- read.delim("appyter_base.txt")

# Extract GO term IDs
base_data$ID <- str_extract(base_data$term, "\\(GO:[0-9]+\\)")
combined_data$ID <- str_extract(combined_data$Term, "\\(GO:[0-9]+\\)")

# Clean up GO term IDs by removing parentheses
base_data$ID <- gsub("[()]", "", base_data$ID)
combined_data$ID <- gsub("[()]", "", combined_data$ID)

# Merge datasets on GO term IDs
merged_data <- merge(base_data, combined_data, by = "ID")

# Data transformation
merged_data$cluster <- as.factor(gsub("Cluster ", "", merged_data$cluster))
merged_data$P.value <- ifelse(merged_data$P.value < 0.05, 1, 0)
merged_data <- merged_data %>%
  group_by(Systems, cluster) %>%
  mutate(filtered = sum(P.value))

merged_data$alpha <- ifelse(merged_data$filtered >= 150, 1, 0.1)

# Define color mapping for Systems (adjust as needed)
color_mapping <- c("Nervous" = "#324a8a", "Immune" = "#68a4b5", "Others" = "#EEEDEB")
merged_data$color <- color_mapping[merged_data$Term]

# Filter data for plotting and calculate means
filtered_data <- merged_data %>% 
  filter(cluster %in% c(0, 13, 1, 20, 19, 2, 10, 8, 16, 17)) %>%
  group_by(cluster) %>%
  summarise(x_mean = mean(x), y_mean = mean(y))

# Adjust transparency based on filtering
merged_data$alpha <- scales::rescale(sqrt(merged_data$filtered), to = c(0.1, 1))

# Generate plot
plot <- ggplot(merged_data, aes(x = x, y = y, color = Systems, alpha = alpha, shape = Systems)) +
  geom_point(size = 7) + 
  geom_label_repel(data = filtered_data, 
                   mapping = aes(x = x_mean, y = y_mean, label = paste("Cluster", cluster)),
                   inherit.aes = FALSE,
                   size = 12,
                   segment.color = "gray90",
                   label.padding = 0.15,
                   box.padding = 0.5,
                   point.padding = 0.5,
                   min.segment.length = 0.1,
                   max.overlaps = 10,
                   label.size = 0.1) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = c("Nervous" = 19, "Immune" = 19, "Others" = 1)) +
  scale_alpha_continuous(range = c(0.1, 1)) +  
  labs(x = "UMAP1", y = "UMAP2", color = "Groups") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.key.size = unit(0.3, "cm"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40))

# Save plot (adjust file name as needed)
ggsave(filename = "clusteres_leiden_3.tiff",
       plot = plot,
       device = "tiff",
       width = 13,
       height = 8,
       units = "in",
       dpi = 400)

# Save output data to Excel file (adjust file name as needed)
write.xlsx(merged_data, "tabela_input_figure_appyter.xlsx")
