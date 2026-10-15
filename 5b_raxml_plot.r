#!/usr/bin/env Rscript
# B3_raxml_plot.r
# Plots a raxml-ng tree with tips colored by species and bootstrap support on nodes
# usage: Rscript B3_raxml_plot.r <tree_file> <sample_info.txt> <outdir> <label> [outgroup]
# outgroup: comma separated sample names; if empty, the tree is midpoint rooted

for (pkg in c("ape", "phangorn")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Missing R package: ", pkg, ". Install with: micromamba install -c conda-forge r-", pkg))
  }
}
library(ape)

args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
sampleinfo_file <- args[2]
outdir <- args[3]
label <- args[4]
outgroup <- ifelse(length(args) >= 5, args[5], "")

tree <- read.tree(tree_file)

# root
if (nchar(outgroup) > 0) {
  og <- intersect(strsplit(outgroup, ",")[[1]], tree$tip.label)
  tree <- root(tree, outgroup = og, resolve.root = TRUE, edgelabel = TRUE)   # edgelabel keeps support on the right branches
  root_note <- paste0("rooted on ", paste(og, collapse = ", "))
} else {
  tree <- phangorn::midpoint(tree, node.labels = "support")   # keeps support on the right branches
  root_note <- "midpoint rooted"
}
tree <- ladderize(tree)

labs <- read.table(sampleinfo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
species <- labs$Species[match(tree$tip.label, labs$Sample)]
species[is.na(species)] <- "Unknown"

species_colors <- c(BCCH = "#1b9e77", CACH = "#d95f02", MOCH = "#7570b3", BOCH = "#e7298a", Unknown = "grey60")
tip_cols <- species_colors[species]

ntip <- length(tree$tip.label)
pdf(file = file.path(outdir, paste0("tree_", label, ".pdf")),
    width = 8, height = max(6, ntip * 0.15), useDingbats = FALSE)
plot(tree, tip.color = tip_cols, cex = 0.5, label.offset = 0.001, no.margin = FALSE,
     main = paste0("raxml-ng, ", label, "\n", root_note))
# bootstrap support: show values >= 50
bs <- suppressWarnings(as.numeric(tree$node.label))
show <- !is.na(bs) & bs >= 50
if (any(show)) {
  nodelabels(text = bs[show], node = ((ntip + 1):(ntip + tree$Nnode))[show],
             frame = "none", adj = c(1.2, -0.4), cex = 0.4)
}
present <- unique(species)
legend("topleft", legend = present, col = species_colors[present], pch = 15, bty = "n", cex = 0.8)
add.scale.bar(cex = 0.6)
dev.off()

write.tree(tree, file.path(outdir, paste0("tree_", label, ".rooted.nwk")))