#!/usr/bin/env Rscript
# B1_PCA_plot3d.r
# 3D PCA (PC1, PC2, PC3) colored by species, from plink output
# usage: Rscript B1_PCA_plot3d.r <plink_prefix> <sample_info.txt> <outdir> <label>
# outputs:
#   pca_<label>_species.3D.html  interactive, rotatable (open in a browser)
#   pca_<label>_species.3D.pdf   static, four viewing angles
# install with: micromamba install -c conda-forge r-plotly r-htmlwidgets r-scatterplot3d r-dplyr pandoc

required_packages <- c("plotly", "htmlwidgets", "scatterplot3d", "dplyr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Missing R package: ", pkg, ". Install with: micromamba install -c conda-forge r-", pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
sampleinfo_file <- args[2]
outdir <- args[3]
label <- args[4]

# read plink output, works with plink 1.9 (no header) and plink2 (#FID IID PC1... or #IID PC1...)
eigenvec_file <- paste0(prefix, ".eigenvec")
first_line <- readLines(eigenvec_file, n = 1)
if (startsWith(first_line, "#")) {
  eigenvec <- read.table(eigenvec_file, header = TRUE, comment.char = "")
  colnames(eigenvec) <- sub("^X\\.", "", colnames(eigenvec))
  colnames(eigenvec)[colnames(eigenvec) == "IID"] <- "Sample"
} else {
  eigenvec <- read.table(eigenvec_file, header = FALSE)
  colnames(eigenvec) <- c("FID", "Sample", paste0("PC", 1:(ncol(eigenvec) - 2)))
}
eigenval <- scan(paste0(prefix, ".eigenval"), quiet = TRUE)
pv <- round((eigenval / sum(eigenval)) * 100, digits = 2)

labs <- read.table(sampleinfo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

pca <- eigenvec %>%
  left_join(labs, by = "Sample") %>%
  mutate(Species = ifelse(is.na(Species), "Unknown", Species))
pca$Species <- factor(pca$Species)

# same colors as the 2D plots
species_colors <- c(BCCH = "#1b9e77", CACH = "#d95f02", MOCH = "#7570b3", BOCH = "#e7298a", Unknown = "grey60")
species_colors <- species_colors[levels(pca$Species)]

xlab <- paste0("PC1 (", pv[1], "%)")
ylab <- paste0("PC2 (", pv[2], "%)")
zlab <- paste0("PC3 (", pv[3], "%)")

# interactive 3D, hover shows sample name
p <- plot_ly(pca, x = ~PC1, y = ~PC2, z = ~PC3,
             color = ~Species, colors = species_colors,
             text = ~Sample, hoverinfo = "text",
             type = "scatter3d", mode = "markers",
             marker = list(size = 4, opacity = 0.85)) %>%
  layout(title = paste0("PCA, ", label),
         scene = list(xaxis = list(title = xlab),
                      yaxis = list(title = ylab),
                      zaxis = list(title = zlab)))

html_file <- file.path(outdir, paste0("pca_", label, "_species.3D.html"))
# selfcontained needs pandoc; fall back to a html + _files folder without it
selfcontained <- nzchar(Sys.which("pandoc"))
saveWidget(p, html_file, selfcontained = selfcontained)

# static 3D, four viewing angles
pdf(file = file.path(outdir, paste0("pca_", label, "_species.3D.pdf")), width = 10, height = 10, useDingbats = FALSE)
par(mfrow = c(2, 2))
for (ang in c(30, 75, 120, 165)) {
  scatterplot3d(pca$PC1, pca$PC2, pca$PC3,
                color = species_colors[as.character(pca$Species)],
                pch = 16, angle = ang, type = "h", lty.hplot = 3,
                xlab = xlab, ylab = ylab, zlab = zlab,
                main = paste0("PCA, ", label, " (angle ", ang, ")"))
  legend("topright", legend = names(species_colors), col = species_colors, pch = 16, bty = "n", cex = 0.8)
}
dev.off()
