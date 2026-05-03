# ============================================================
# FROH by chromosome comparison (pop1 vs pop2)
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. LOAD DATA
# ------------------------------------------------------------
FROH_chrF <- read.csv("FROH_CHROMOSOMEF.txt")
FROH_chrCRA <- read.csv("FROH_CHROMOSOMECRA.txt")

# ------------------------------------------------------------
# 2. FUNCTION TO SUMMARIZE FROH
# ------------------------------------------------------------
summarise_froh <- function(df) {
  df %>%
    pivot_longer(
      cols = starts_with("Chr_"),
      names_to = "Chr",
      values_to = "value"
    ) %>%
    group_by(Chr, group) %>%
    summarise(
      MEAN = mean(value, na.rm = TRUE),
      SE = sd(value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}

summary_F <- summarise_froh(FROH_chrF)
summary_CRA <- summarise_froh(FROH_chrCRA)

# ------------------------------------------------------------
# 3. COMBINE DATASETS
# ------------------------------------------------------------
df_combinado <- bind_rows(summary_CRA, summary_F)

df_combinado$Chr <- as.numeric(gsub("Chr_", "", df_combinado$Chr))

# ------------------------------------------------------------
# 4. PLOT
# ------------------------------------------------------------
p <- ggplot(df_combinado, aes(x = Chr, y = MEAN, fill = group)) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(
    aes(ymin = MEAN - SE, ymax = MEAN + SE),
    position = position_dodge(0.9),
    width = 0.3,
    colour = "black"
  ) +
  scale_fill_manual(values = c("purple", "orange")) +
  scale_x_continuous(breaks = 1:29) +
  labs(
    x = "Chromosome",
    y = expression("Mean F"["ROH"])
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.grid = element_blank()
  )

print(p)

# ------------------------------------------------------------
# 5. SAVE FIGURE 
# ------------------------------------------------------------
ggsave("FROH_by_chromosome.png",
       plot = p,
       width = 10,
       height = 6,
       dpi = 300)
