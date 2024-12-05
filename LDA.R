# LDA for Meta-analysis MDD Blood DEGs Synapse
# Author: Anny Adri - PhD Student, University of Sao Paulo (USP)
# References: https://www.r-bloggers.com/2021/05/linear-discriminant-analysis-in-r/


# Set working directory
setwd("/path/to/your/directory")

# Load necessary libraries
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(klaR)
library(psych)
library(MASS)
library(ggord)
library(devtools)
library(svglite)
library(ggplot2)
library(caret)
library(e1071)
library(openxlsx)

# Load data
df <- read_excel("/path/to/tabela_input_tudo_dotplot_ontologias_sinapse.xlsx", sheet = 2)
counts_GSE85855_Cathomas <- read_excel("/path/to/data_GSE85855_Cathomas.xlsx")
counts_PMID_Trang <- read.delim("/path/to/Counts_PMID_Trang_rnaseq.tab")
metadata_GSE85855 <- read_excel("/path/to/metadata_GSE85855.xlsx")
metadata_PMID_Trang <- read_excel("/path/to/metadata_PMID_Trang.xlsx")

# Process the genes from the synapse table
df_genes <- df %>%
  mutate(Genes = str_split(Genes, ";")) %>%
  unnest(Genes) %>%
  distinct(Genes)

# Rename 'Genes' column for merging
counts_GSE85855_Cathomas <- counts_GSE85855_Cathomas %>% rename(Genes = gene)
counts_PMID_Trang <- counts_PMID_Trang %>% rename(Genes = gene)

# Merge gene data with count tables
cathomas <- merge(df_genes, counts_GSE85855_Cathomas, by = "Genes")
trang <- merge(df_genes, counts_PMID_Trang, by = "Genes")

###------------------CATHOMAS Analysis------------------###

# Set working directory for Cathomas data
setwd("/path/to/Cathomas")

# Data pre-processing for MDD immune system analysis
df_teste_cathomas <- cathomas

# Filter out low-count values
col_zero <- colSums(df_teste_cathomas < 10)
row_zero <- rowSums(df_teste_cathomas < 10)
excl_col <- which(col_zero >= 12)
excl_row <- which(row_zero >= 5)
df_teste1_cathomas <- df_teste_cathomas[-excl_row, -excl_col]

# Transpose the data
df_teste1_cathomas <- as.data.frame(t(df_teste1_cathomas))
colnames(df_teste1_cathomas) <- df_teste1_cathomas[1,]
df_teste1_cathomas <- df_teste1_cathomas[-1,]
df_teste1_cathomas[1:62] <- as.data.frame(sapply(df_teste1_cathomas[1:62], as.numeric))

# Add group labels
colnames(df_teste1_cathomas) <- make.unique(colnames(df_teste1_cathomas))
df_teste1_cathomas <- tibble::rownames_to_column(df_teste1_cathomas, "Groups")

# Map MDD and control groups
group_map <- setNames(as.character(metadata_GSE85855$Group), metadata_GSE85855$id)
df_teste1_cathomas$Groups <- group_map[as.character(df_teste1_cathomas$Groups)]

# Perform Linear Discriminant Analysis (LDA)
set.seed(123)
ind <- sample(2, nrow(df_teste1_cathomas), replace = TRUE, prob = c(0.7, 0.3))
training <- df_teste1_cathomas[ind == 1,]
testing <- df_teste1_cathomas[ind == 2,]

# Linear Discriminant Analysis
linear <- lda(Groups ~ ., data = training, prior = c(7,3)/10)

# Predictions and plotting
p <- predict(linear, training)

# Plot histogram of LDA results
tiff("BarplotLD1_cathomas.tiff", width = 700, height = 1100, units = "px", res = 200)
ldahist(data = p$x[,1], g = training$Groups, col = "#E7E6E6")
dev.off()

# Create dataframe with LDA results
resultados_lda_cathomas <- data.frame(Group = training$Groups, LD1 = p$x[,1])

# Plot density plot for groups
cores <- c("MDD" = "#C8A4E6", "control" = "#83BCD4")
density_plot <- ggplot(resultados_lda_cathomas, aes(x = LD1, fill = factor(Group))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = cores) +
  labs(title = "Cathomas, F. et al., 2022", x = "LD1", y = "Density", fill = "Group") +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_line(color = "gray86", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(size = 15),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14)
  ) +
  scale_x_continuous(limits = c(-4, 4))

# Save and display the plot
ggsave("histogram_cathomas.tiff", density_plot, device = "tiff", scale = 1, width = 4, height = 3, units = "in")
print(density_plot)

# Confusion matrix and accuracy for training data
p1 <- predict(linear, training)$class
tab <- table(Predicted = p1, Actual = training$Groups)
accuracy_training <- sum(diag(tab)) / sum(tab)

# Confusion matrix and accuracy for testing data
p2 <- predict(linear, testing)$class
tab1 <- table(Predicted = p2, Actual = testing$Groups)
accuracy_testing <- sum(diag(tab1)) / sum(tab1)

# Create a dataframe for confusion matrix results
results <- data.frame(
  Metric = c("Training Accuracy", "Testing Accuracy"),
  Accuracy = c(accuracy_training, accuracy_testing),
  True_Positive_Training = c(tab[2, 2], NA),
  True_Negative_Training = c(tab[1, 1], NA),
  False_Positive_Training = c(tab[2, 1], NA),
  False_Negative_Training = c(tab[1, 2], NA),
  True_Positive_Testing = c(NA, tab1[2, 2]),
  True_Negative_Testing = c(NA, tab1[1, 1]),
  False_Positive_Testing = c(NA, tab1[2, 1]),
  False_Negative_Testing = c(NA, tab1[1, 2])
)

# Save results to an Excel file
write_xlsx(results, "lda_results_cathomas.xlsx")

# 1. Get LDA coefficients
coefficients_cathomas <- linear$scaling

# 2. Create a dataframe with genes and coefficients
resultados_genes_cathomas <- data.frame(Gene = rownames(coefficients_cathomas), Coefficient = as.vector(coefficients_cathomas))

# 3. Sort genes by magnitude of coefficients
resultados_genes_cathomas <- resultados_genes_cathomas %>%
  arrange(desc(abs(Coefficient)))

# 4. Save results to Excel
write.xlsx(resultados_genes_cathomas, "resultados_genes_lda.xlsx", rowNames = TRUE)

# Plot barplot for genes
comuns_pos <- comuns %>% filter(`LD1 Cathomas` > 0) %>% arrange(desc(`LD1 Cathomas`))
comuns_neg <- comuns %>% filter(`LD1 Cathomas` <= 0) %>% arrange(desc(`LD1 Cathomas`))
comuns_ord <- bind_rows(comuns_pos, comuns_neg)

p1 <- ggplot(comuns_ord, aes(x = reorder(Gene, `LD1 Cathomas`), y = `LD1 Cathomas`)) +
  geom_col(fill = "grey85", color = "black", size = 0.1) +
  labs(x = "Genes", y = "LD1 value", title = "Cathomas, F. et al., 2022") +
  theme_bw() +
  ylim(-5.583e-02, 0.04) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(size = 15),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14, angle = 90, hjust = 1)
  )

# Save and display barplot
ggsave("barplot_gene_cathomas.tiff", plot = p1, device = "tiff", width = 10.5, height = 4, units = "in")
print(p1)

###------------------Trang Analysis------------------###


# Set working directory for Trang data
setwd("/path/to/Trang")

# Data processing for Trang
df_teste_trang <- trang

# Filter low-count values for Trang
col_zero <- colSums(df_teste_trang < 10)
row_zero <- rowSums(df_teste_trang < 10)
excl_col <- which(col_zero >= 12)
excl_row <- which(row_zero >= 5)
df_teste1_trang <- df_teste_trang[-excl_row, -excl_col]

# Transpose and further processing
df_teste1_trang <- as.data.frame(t(df_teste1_trang))
colnames(df_teste1_trang) <- df_teste1_trang[1,]
df_teste1_trang <- df_teste1_trang[-1,]
df_teste1_trang[1:62] <- as.data.frame(sapply(df_teste1_trang[1:62], as.numeric))

# Add group labels
colnames(df_teste1_trang) <- make.unique(colnames(df_teste1_trang))
df_teste1_trang <- tibble::rownames_to_column(df_teste1_trang, "Groups")

# Map group names from metadata
group_map <- setNames(as.character(metadata_PMID_Trang$Group), metadata_PMID_Trang$id)
df_teste1_trang$Groups <- group_map[as.character(df_teste1_trang$Groups)]

# Perform LDA for Trang
set.seed(123)
ind <- sample(2, nrow(df_teste1_trang), replace = TRUE, prob = c(0.7, 0.3))
training <- df_teste1_trang[ind == 1,]
testing <- df_teste1_trang[ind == 2,]

# Linear Discriminant Analysis for Trang
linear <- lda(Groups ~ ., data = training, prior = c(7,3)/10)

# Predictions for Trang and plot results
p <- predict(linear, training)
tiff("BarplotLD1_trang.tiff", width = 700, height = 1100, units = "px", res = 200)
ldahist(data = p$x[,1], g = training$Groups, col = "#E7E6E6")
dev.off()

# Create dataframe for LDA results
resultados_lda_trang <- data.frame(Group = training$Groups, LD1 = p$x[,1])

# Plot density for Trang groups
density_plot <- ggplot(resultados_lda_trang, aes(x = LD1, fill = factor(Group))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = cores) +
  labs(title = "Trang, 2015", x = "LD1", y = "Density", fill = "Group") +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_line(color = "gray86", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(size = 15),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14)
  ) +
  scale_x_continuous(limits = c(-4, 4))

# Save plot and print
ggsave("histogram_trang.tiff", density_plot, device = "tiff", scale = 1, width = 4, height = 3, units = "in")
print(density_plot)

# Confusion matrix and accuracy for Trang
p1 <- predict(linear, training)$class
tab <- table(Predicted = p1, Actual = training$Groups)
accuracy_training <- sum(diag(tab)) / sum(tab)

p2 <- predict(linear, testing)$class
tab1 <- table(Predicted = p2, Actual = testing$Groups)
accuracy_testing <- sum(diag(tab1)) / sum(tab1)

# Save confusion matrix results for Trang
results <- data.frame(
  Metric = c("Training Accuracy", "Testing Accuracy"),
  Accuracy = c(accuracy_training, accuracy_testing),
  True_Positive_Training = c(tab[2, 2], NA),
  True_Negative_Training = c(tab[1, 1], NA),
  False_Positive_Training = c(tab[2, 1], NA),
  False_Negative_Training = c(tab[1, 2], NA),
  True_Positive_Testing = c(NA, tab1[2, 2]),
  True_Negative_Testing = c(NA, tab1[1, 1]),
  False_Positive_Testing = c(NA, tab1[2, 1]),
  False_Negative_Testing = c(NA, tab1[1, 2])
)

# Save results for Trang to Excel
write_xlsx(results, "lda_results_trang.xlsx")