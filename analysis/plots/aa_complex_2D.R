library(tidyverse)
library(ggsci)



# Check current working directory
getwd()

AA = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'Z')
Real = c(0.47301727329604, 0.8538094956668423, 0.11828499793341991, 0.9747755230028201, 0.8199014520800844, 0.1606248385548573, 0.6757813367907463, 0.7655475910683724, 0.22845927297235719, 0.6434182427720703, 0.2227242052954308, 0.3110071300912875, 0.885539733325962, 0.882696763865174, 0.7096238831354776, 0.6417589349110779, 0.6306863593294083, 0.5269649336844214, 0.2944342578286395, 0.06376658936235352, 0.3176178122104376, 0.7492495377328181)
Imaginary = c(0.7232536733082316, 0.5697059160384558, 1, 0.5499076221718956, 0.43284167492105957, 0.3406523494456186, 0.9784576691706325, 0.1006221953193782, 0.651919464341051, 0.45469533326941014, 0.5739188104178984, 0.5104547567892527, 0.6498448427150636, 0.9408035208616283, 0.39104919843056674, 0.29915271002578964, 0.6716984414123757, 0.7836565106841363, 0.6813887581916396, 0, 0.1755879467197631, 0.4700787706762518)
hydrophobic = c("A", "V", "I", "L", "M", "F", "Y", "W")
  
df <- data.frame(AA, Real, Imaginary)
df$HydroGroup <- ifelse(df$AA %in% hydrophobic, "Hydrophobic", "non-Hydrophobic")

df$HydroGroup = factor(df$HydroGroup, levels = c("non-Hydrophobic", "Hydrophobic"))

# plot 
p = df %>%
    ggplot(aes(x = Real, y = Imaginary, color = HydroGroup)) +
    geom_point(size = 3) +
    geom_text(aes(label = AA), vjust = -1, size = 3) +
    scale_color_jama() +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
  
    )
p

# Save plot as PDF
ggsave("aa_complex_2d.pdf", plot = p, width = 7, height = 5)
ggsave("aa_complex_2d.png", plot = p, width = 7, height = 5, dpi = 300)

# Ref
# https://en.wikipedia.org/wiki/Proteinogenic_amino_acid#/media/File:Proteinogenic_Amino_Acid_Table.png
# https://en.wikipedia.org/wiki/Amino_acid#/media/File:ProteinogenicAminoAcids.svg 

