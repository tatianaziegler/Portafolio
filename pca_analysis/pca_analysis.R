# =========================================================
# Population Structure Analysis (PCA)
# POP1 vs POP2 Populations
# Author: Tatiana Ziegler
# =========================================================

# ==============================
# Libraries
# ==============================
library(ggplot2)
library(tidyverse)
library(plotly)

# =========================================================
# 1. SNP DATA PREPARATION (PLINK assumed pre-run)
# =========================================================
# NOTE: Genotypes obtained using PLINK (--recode A)

pop1 <- read.table("POP1.raw", header = TRUE)
pop2 <- read.table("POP2.raw", header = TRUE)

# Remove metadata columns
pop1_geno <- pop1[, -(1:6)]
pop2_geno <- pop2[, -(1:6)]

# =========================================================
# 2. SNP FILTERING
# =========================================================
pop1_geno <- as.matrix(pop1_geno)
pop2_geno <- as.matrix(pop2_geno)

pop1_geno[is.na(pop1_geno)] <- 0
pop2_geno[is.na(pop2_geno)] <- 0

# Keep only SNPs with variance
combined_geno <- rbind(pop1_geno, pop2_geno)

var_snp <- apply(combined_geno, 2, var)
combined_geno <- combined_geno[, var_snp > 0]

# =========================================================
# 3. PCA ANALYSIS
# =========================================================
pca_result <- prcomp(combined_geno, center = TRUE, scale. = TRUE)

var_exp <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

# =========================================================
# 4. SAMPLE ANNOTATION
# =========================================================
pop_info <- c(
  rep("POP1", nrow(pop1)),
  rep("POP2", nrow(pop2))
)

pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  PC3 = pca_result$x[,3],
  Population = pop_info
)

# =========================================================
# 5. 2D PCA PLOT
# =========================================================
ggplot(pca_df, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme_minimal() +
  labs(
    title = "Population Structure Analysis (PCA)",
    x = paste0("PC1 (", round(var_exp[1], 2), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 2), "%)")
  ) +
  theme(legend.position = "top")

# =========================================================
# 6. 3D PCA PLOT
# =========================================================
fig <- plot_ly(
  pca_df,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Population,
  type = "scatter3d",
  mode = "markers"
)

fig <- fig %>% layout(
  title = "3D PCA – Population Structure",
  scene = list(
    xaxis = list(title = paste0("PC1 (", round(var_exp[1], 2), "%)")),
    yaxis = list(title = paste0("PC2 (", round(var_exp[2], 2), "%)")),
    zaxis = list(title = paste0("PC3 (", round(var_exp[3], 2), "%)"))
  )
)

fig
