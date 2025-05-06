library(tidyverse)

# Check current working directory
getwd()

# evoprot
df1 <- read.table("../benchmark/PA-SigPro/acc_percent_sigprot.txt")[, c(1, 3)]
colnames(df1) <- c("filename", "value")

# plm 
df2 <- read.table("../benchmark/AI-Embedding/acc_percent_AI_embd.txt")[, c(1, 3)]
colnames(df2) <- c("filename", "value")

# 
df = rbind(df1, df2)

# Extract 'parameter' and split into 'type' and 'group'
df$parameter <- sub("_uniprot.*", "", df$filename)
df <- df %>%
  separate(parameter, into = c("type", "group"), sep = "_")

# Extract and convert size to "K" format
df$sizeK <- sub(".*_(\\d+)_.*", "\\1", df$filename)
df$sizeK_num <- as.numeric(df$sizeK) / 1000
df$sizeK <- paste0(df$sizeK_num, "K")

# Order sizeK factor by numeric value
df$sizeK <- factor(df$sizeK, levels = paste0(sort(unique(df$sizeK_num)), "K"))

# Create plot
p <- df %>%
  ggplot(aes(x = sizeK, y = value, fill = group)) +
  geom_col(position = "dodge") +
  facet_grid(~type) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Database",
    y = "Accuracy (%)",
    fill = "Group")

# Display the plot (optional)
print(p)

# Save plot as PDF
ggsave("../../plots/dimreduct_distfunc.pdf", plot = p, width = 8, height = 5)
ggsave("../../plots/dimreduct_distfunc.png", plot = p, width = 8, height = 5, dpi = 300)

df %>%
  group_by(type, group) %>%
  summarise(avg_value = mean(value), .groups = "drop") %>%
  arrange(desc(avg_value))

