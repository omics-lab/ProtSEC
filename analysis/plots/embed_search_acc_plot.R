library(tidyverse)
library(ggsci)

# Check current working directory
getwd()

# evoprot
df1 <- read.table("../benchmark/PA-SigPro/acc_percent_sigprot.txt")[, c(1, 3)]
colnames(df1) <- c("filename", "value")

# plm 
df2 <- read.table("../benchmark/AI-Embedding/acc_percent_AI_embd.txt")[, c(1, 3)]
colnames(df2) <- c("filename", "value")

# blastp 
df3 <- read.table("../benchmark/blastp/acc_percent_blastp.txt")[, c(1, 3)]
colnames(df3) <- c("filename", "value")

# 
df = rbind(df1, df2, df3)

df$new_column <- str_extract(df$filename, "(?<=sprot_).*(?=_result)")
df$new_column <- str_replace(df$new_column, fixed(".fasta_db.pkl_PA-SigPro"), "_ComProt")
df$new_column <- str_remove(df$new_column, "^5000_vs_uniprot_sprot_")
df$new_column <- str_remove(df$new_column, "_db_blastp$")

df <- separate(df, new_column, into = c("dataset_size", "method"), sep = "_", extra = "merge") %>%
  filter(dataset_size != 80000)

df$dataset_size <- df$dataset_size %>%
  as.numeric() %>%
  sapply(function(x) ifelse(x >= 1000, paste0(x / 1000, "K"), as.character(x)))

df$dataset_size = factor(df$dataset_size, levels = c("5K",  "10K", "20K", "40K"))

# Create plot
p = df %>%
  ggplot(aes(x = method, y = value, fill = method)) +
  geom_col(position = "dodge") +
  facet_grid(~dataset_size) +
  scale_fill_jama() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "",
    y = "Protein similary search accuracy (%)")

# Display the plot (optional)
print(p)

# Save plot as PDF
ggsave("embed_search_acc.pdf", plot = p, width = 8, height = 5)
ggsave("embed_search_acc.png", plot = p, width = 8, height = 5, dpi = 300)

