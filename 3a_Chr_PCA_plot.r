#!/usr/bin/env Rscript
# 3a_Chr_PCA_plot.r
# Plots plink PCA output colored by species
# usage: Rscript 3a_Chr_PCA_plot.r <plink_prefix> <sample_info.txt> <outdir> <label>

# Load required packages, installing if necessary
required_packages <- c("ggplot2", "dplyr", "readr")
installed_packages <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
sampleinfo_file <- args[2]
outdir <- args[3]
label <- args[4]

# read plink output
eigenvec <- read.table(paste0(prefix, ".eigenvec"), header = FALSE)
colnames(eigenvec) <- c("FID", "Sample", paste0("PC", 1:(ncol(eigenvec) - 2)))
eigenval <- scan(paste0(prefix, ".eigenval"))

# percent variation explained (of the PCs computed)
pv <- round((eigenval / sum(eigenval)) * 100, digits = 2)

# read species assignments
labs <- read.table(sampleinfo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

pca <- eigenvec %>%
  left_join(labs, by = "Sample") %>%
  mutate(Species = ifelse(is.na(Species), "Unknown", Species))

if (any(pca$Species == "Unknown")) {
  cat("Warning: samples without a species assignment:\n")
  cat(pca$Sample[pca$Species == "Unknown"], sep = "\n")
}

pca$Species <- factor(pca$Species)

# colors for each species
species_colors <- c(BCCH = "#1b9e77", CACH = "#d95f02", MOCH = "#7570b3", BOCH = "#e7298a", Unknown = "grey60")

write.table(pca, file.path(outdir, paste0("pca_", label, "_species.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_pcs <- function(x, y) {
  ggplot(pca, aes(x = .data[[paste0("PC", x)]], y = .data[[paste0("PC", y)]], color = Species)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_color_manual(values = species_colors) +
    labs(x = paste0("PC", x, " (Percent Variation = ", pv[x], "%)"),
         y = paste0("PC", y, " (Percent Variation = ", pv[y], "%)"),
         title = paste0("PCA, ", label)) +
    theme_classic()
}

pdf(file = file.path(outdir, paste0("pca_", label, "_species.PC1_PC2.pdf")), width = 7, height = 6, useDingbats = FALSE)
print(plot_pcs(1, 2))
dev.off()

pdf(file = file.path(outdir, paste0("pca_", label, "_species.PC1_PC3.pdf")), width = 7, height = 6, useDingbats = FALSE)
print(plot_pcs(1, 3))
dev.off()

pdf(file = file.path(outdir, paste0("pca_", label, "_species.PC2_PC3.pdf")), width = 7, height = 6, useDingbats = FALSE)
print(plot_pcs(2, 3))
dev.off()

# scree plot
pdf(file = file.path(outdir, paste0("pca_", label, "_scree.pdf")), width = 6, height = 4, useDingbats = FALSE)
barplot(pv, names.arg = paste0("PC", seq_along(pv)), ylab = "Percent Variation", las = 2)
dev.off()
