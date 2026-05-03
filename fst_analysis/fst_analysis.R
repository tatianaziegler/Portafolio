# =========================================================
# FST Analysis – Selection Signatures in Goat Populations
# Author: Tatiana Ziegler
# =========================================================

# ==============================
# Libraries
# ==============================
library(vcfR)
library(adegenet)
library(hierfstat)
library(ggplot2)
library(dplyr)

# ==============================
# Load VCF data
# ==============================
pop1_vcf <- read.vcfR("POP1.vcf")
pop2_vcf <- read.vcfR("POP2.vcf")

# ==============================
# Convert to genind format
# ==============================
pop1_genind <- vcfR2genind(pop1_vcf)
pop2_genind <- vcfR2genind(pop2_vcf)

# ==============================
# Identify shared SNPs
# ==============================
common_snps <- intersect(colnames(pop1_genind@tab),
                         colnames(pop2_genind@tab))

pop1_genind <- pop1_genind[, common_snps]
pop2_genind <- pop2_genind[, common_snps]

# ==============================
# Merge populations
# ==============================
combined_tab <- rbind(pop1_genind@tab, pop2_genind@tab)

combined_genind <- df2genind(
  combined_tab,
  pop = rep(c("POP1", "POP2"),
            c(nInd(pop1_genind), nInd(pop2_genind))),
  ploidy = 2,
  ncode = 1
)

# ==============================
# Convert to hierfstat format
# ==============================
genpop_data <- genind2hierfstat(combined_genind)

# ==============================
# FST estimation
# ==============================
fst_results <- varcomp.glob(
  genpop_data,
  levels = pop(combined_genind)
)

# ==============================
# Extract per-locus FST
# ==============================
fst_loc_df <- data.frame(
  SNP = rownames(fst_results$loc)[-1],
  Fst = fst_results$loc[-1, 1]
)

# ==============================
# Save results
# ==============================
write.table(fst_loc_df,
            "fst_per_locus.txt",
            row.names = FALSE,
            sep = "\t")

# ==============================
# FST distribution plot
# ==============================
p1 <- ggplot(fst_loc_df, aes(x = Fst)) +
  geom_histogram(binwidth = 0.001,
                 fill = "steelblue",
                 alpha = 0.7) +
  theme_minimal() +
  labs(title = "FST Distribution (POP1 vs POP2)",
       x = "FST",
       y = "Frequency")

ggsave("fst_distribution.png", p1, width = 6, height = 4)

# ==============================
# Outlier detection (top 1%)
# ==============================
threshold <- quantile(fst_loc_df$Fst, 0.99)

fst_outliers <- subset(fst_loc_df, Fst >= threshold)

write.table(fst_outliers,
            "fst_outliers.txt",
            row.names = FALSE,
            sep = "\t")
