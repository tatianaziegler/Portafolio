# =========================================================
# Manhattan Plot - Selection Scan (FST)
# POP1 vs POP2
# =========================================================

library(ggplot2)
library(dplyr)

# =========================================================
# 1. LOAD DATA
# =========================================================
fst_data <- read.table("fst_data.txt", header = TRUE)

# Expected columns:
# SNP | CHR | BP | FST

fst_data$CHR <- as.factor(fst_data$CHR)

# =========================================================
# 2. ORDER CHROMOSOMES
# =========================================================
fst_data <- fst_data %>%
  arrange(CHR, BP) %>%
  mutate(pos_index = row_number())

# =========================================================
# 3. CHROMOSOME CENTERS
# =========================================================
chromosome_center <- fst_data %>%
  group_by(CHR) %>%
  summarise(center = (max(pos_index) + min(pos_index)) / 2)

# =========================================================
# 4. MANHATTAN PLOT
# =========================================================
ggplot(fst_data, aes(x = pos_index, y = FST, color = CHR)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_x_continuous(
    label = chromosome_center$CHR,
    breaks = chromosome_center$center
  ) +
  theme_minimal() +
  labs(
    title = "Manhattan Plot - Genetic Differentiation (FST)",
    x = "Chromosome",
    y = "FST"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
