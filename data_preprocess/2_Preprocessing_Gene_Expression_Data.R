# -*- UTF-8 -*-
# !/usr/bin/env R4.2
# by zhiyuan

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
install.packages("snp.plotter")

library("snp.plotter")
library(dplyr)

# Read expression data
ex_df <- read.table("/work/data/gene.expressed.txt", row.names = 1, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# Filter expression data (TPM > 1 in 10% lines)
ex_df_filter <- ex_df[(rowSums(ex_df > 1) >= ncol(ex_df) * 0.1), ]

# Normal quantile transformation
ex_df_normal <- t(apply(ex_df_filter, 1, function(i) (qqnorm(i, plot = F))$x))
colnames(ex_df_normal) <- colnames(ex_df_filter)

# Save transformed data
write.table(ex_df_normal, file = "gene.expressed.normal.txt", sep = "\t", quote = FALSE)

# Expression confounder
## 1 Expression PCA factor
spearman_cor <- cor(ex_df_normal, method = "spearman")
pca <- princomp(spearman_cor, cor = T)
pca_sum <- summary(pca, loadings = T)
s <- pca_sum$sdev
proption <- s^2 / sum(s^2)
num_pc <- seq(1:length(proption))
pca_pve <- cbind(num_pc, proption)
head(pca_pve)
plot(pca_pve[1:15, 1], pca_pve[1:15, 2] * 100, type = "o", main = "contribution_variance_pheno368_spearman", xlab = "PCAs", ylab = "% variance")
scores <- pca_sum$scores
pca.scores <- scores[, 1:10] # select top 10 PCs
colnames(pca.scores.df) <- c("PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15")

## 2 Population structure
system("plink --bfile /work/data/gene --pca 5")
pca_pop <- read.table("plink.eigenvec", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
pca_pop <- pca_pop[, 2:6]
colnames(pca_pop) <- c("sample", "PC1", "PC2", "PC3", "PC4", "PC5")

pca.scores.df$sample <- rownames(pca.scores)
fin_cov = t(merge(pca.scores.df, pca_pop, by.x = "sample", by.y = "sample"))

write.table(fin_cov, "fin_cov_10HF_5PS.txt", sep = "\t", quote = FALSE, row.names = FALSE)
