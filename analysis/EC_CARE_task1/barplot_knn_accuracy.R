library(ggplot2)
library(readr)
library(gridExtra)

setwd("//wsl$/Ubuntu/home/rashedul/project/ProSEC/analysis/EC_CARE_task1/")

# Read both summary CSVs
df1 <- read.csv("knn_accuracy_summary_sample_30-50_5000.csv")
df2 <- read.csv("knn_accuracy_summary_sample_30-50_2000.csv")

# Standardize model names for both
df1$Model <- gsub("^esm2large$", "esm2_large", df1$Model)
df1$Model <- gsub("^esm2small$", "esm2_small", df1$Model)
df1$Model <- gsub("^protbert$", "prot_bert", df1$Model)
df1$Model <- gsub("^prott5$", "prot_t5", df1$Model)

df2$Model <- gsub("^esm2large$", "esm2_large", df2$Model)
df2$Model <- gsub("^esm2small$", "esm2_small", df2$Model)
df2$Model <- gsub("^protbert$", "prot_bert", df2$Model)
df2$Model <- gsub("^prott5$", "prot_t5", df2$Model)

method_order <- c("ProtSEC", "esm2_large", "esm2_small", "prot_t5", "prot_bert")
df1$Model <- factor(df1$Model, levels = method_order)
df2$Model <- factor(df2$Model, levels = method_order)

# Make sure TestSet factor levels are in the same order for both dataframes
testset_order <- c("30-50_protein_test.csv", "30_protein_test.csv")
df1$TestSet <- factor(df1$TestSet, levels = testset_order)
df2$TestSet <- factor(df2$TestSet, levels = testset_order)

p1 <- ggplot(df1, aes(x = Model, y = Accuracy, fill = TestSet)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  labs(title = "Sample 5000", x = "", y = "Accuracy (%)") +
  scale_fill_manual(values = paired_hex) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "gray", size = 1),
        strip.background = element_rect(fill="white"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black", vjust = 0.5, hjust=1))

p2 <- ggplot(df2, aes(x = Model, y = Accuracy, fill = TestSet)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  labs(title = "Sample 2000", x = "", y = "Accuracy (%)") +
  scale_fill_manual(values = paired_hex) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "gray", size = 1),
        strip.background = element_rect(fill="white"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black", vjust = 0.5, hjust=1))

# Arrange plots side by side
combined_plot <- grid.arrange(p1, p2, ncol = 2)

# Save combined plot
ggsave("../plots/knn_accuracy_barplot_combined.pdf", plot = combined_plot, width = 12, height = 6)


