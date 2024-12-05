# Binomial Regression Analysis
# Author: Anny Adri - PhD Student, University of Sao Paulo (USP)

# Genes of interest for analysis
genes <- c("ADORA3", "FBXO2", "FXYD6", "RPS28", "RPS9")

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(broom)
library(gridExtra)
library(readxl)
library(writexl)

# Load data
cathomas_regressao <- read_excel("cathomas_regressao.xlsx")
trang_regressao <- read_excel("trang_regressao.xlsx")

# Filter data to include only relevant columns
filtered_data_cathomas <- cathomas_regressao %>% select(Groups, all_of(genes))
filtered_data_trang <- trang_regressao %>% select(Groups, all_of(genes))

# Save filtered data as Excel files
write.xlsx(filtered_data_cathomas, "input_cathomas.xlsx")
write.xlsx(filtered_data_trang, "input_trang.xlsx")

# Plotting and Analysis function
plot_logistic <- function(gene, data, results_df) {
  # Perform logistic regression for the gene
  formula <- as.formula(paste("as.numeric(Groups == 'MDD') ~", gene))
  model <- glm(formula, data = data, family = binomial)
  
  # Extract p-value for the gene
  p_value <- summary(model)$coefficients[2, "Pr(>|z|)"]
  
  # Add results to global dataframe
  results_df <<- rbind(results_df, data.frame(Gene = gene, P_Value = p_value))
  
  # Create the plot
  ggplot(data, aes_string(x = gene, y = "as.numeric(Groups == 'MDD')", color = "Groups")) +
    geom_point() +
    stat_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, color = "grey40") +
    labs(x = "", y = "Probability", title = paste(gene, "\np-value =", round(p_value, 6))) +
    scale_color_manual(values = c("MDD" = '#b492d1', "control" = '#8bb1c2')) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 16),  # Axis x title font size
      axis.title.y = element_text(size = 16),  # Axis y title font size
      axis.text.x = element_text(size = 16, angle = 45, hjust = 1),  # X-axis labels font size
      axis.text.y = element_text(size = 16),
      legend.position = "none"  # Remove legend
    )
}

########### CATHOMAS DATA ANALYSIS ############################

# Initialize an empty dataframe for storing results
results_cathomas <- data.frame(Gene = character(), P_Value = numeric(), stringsAsFactors = FALSE)

# Apply logistic regression for each gene and plot
lapply(genes, function(gene) plot_logistic(gene, filtered_data_cathomas, results_cathomas))

# Filter significant results (p < 0.05)
significant_results_cathomas <- results_cathomas[results_cathomas$P_Value < 0.05, ]

# Apply function to generate plots for significant results only
plots_cathomas <- lapply(significant_results_cathomas$Gene, function(gene) plot_logistic(gene, filtered_data_cathomas, results_cathomas))

# Arrange the plots in a grid and save as TIFF and SVG
p_cathomas <- do.call("grid.arrange", c(plots_cathomas, ncol = 5))
ggsave("plots_cathomas_5_genes.tiff", plot = p_cathomas, width = 16, height = 4, dpi = 300, compression = "lzw")
ggsave("plots_cathomas_5_genes.svg", plot = p_cathomas, width = 16, height = 4, dpi = 300)

# Save significant results as an Excel file
write.xlsx(significant_results_cathomas, "results_cathomas_significant.xlsx")

# Display significant results
print(significant_results_cathomas)

########### TRANG DATA ANALYSIS ############################

# Gene for analysis in Trang dataset
genes_trang <- "LRFN4"

# Initialize an empty dataframe for storing results
results_trang <- data.frame(Gene = character(), P_Value = numeric(), stringsAsFactors = FALSE)

# Apply logistic regression and plot for the specified gene
lapply(genes_trang, function(gene) plot_logistic(gene, filtered_data_trang, results_trang))

# Filter significant results (p < 0.05)
significant_results_trang <- results_trang[results_trang$P_Value < 0.05, ]

# Apply function to generate plots for significant results only
plots_trang <- lapply(significant_results_trang$Gene, function(gene) plot_logistic(gene, filtered_data_trang, results_trang))

# Arrange the plots in a grid and save as TIFF
p_trang <- do.call("grid.arrange", c(plots_trang, ncol = 5))
ggsave("plots_trang.tiff", plot = p_trang, width = 15.1, height = 3.3, dpi = 300, compression = "lzw")

# Save significant results as an Excel file
write.xlsx(significant_results_trang, "results_trang_significant.xlsx")

# Display significant results
print(significant_results_trang)

########### LEGEND PLOT ############################

# Create a plot with just the legend
legend_plot <- ggplot(filtered_data_trang, aes_string(x = genes_trang, y = "as.numeric(Groups == 'MDD')", color = "Groups")) +
  geom_point() +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, color = "grey40") +
  scale_color_manual(values = c("MDD" = '#b492d1', "Control" = '#8bb1c2')) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  ) + 
  guides(color = guide_legend(title = "Groups"))

# Extract legend only
legend_only <- legend_plot + geom_blank()
print(legend_only)

# Save the legend plot as a TIFF file
ggsave("legenda.tiff", plot = legend_only, width = 5, height = 5, dpi = 300, compression = "lzw")
