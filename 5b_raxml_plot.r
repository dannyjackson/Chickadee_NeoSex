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

# same tree with branch length values printed on each branch and a substitution-distance axis
pdf(file = file.path(outdir, paste0("tree_", label, ".brlen.pdf")),
    width = 10, height = max(6, ntip * 0.18), useDingbats = FALSE)
plot(tree, tip.color = tip_cols, cex = 0.5, label.offset = 0.001, no.margin = FALSE,
     main = paste0("raxml-ng, ", label, " (branch lengths)\n", root_note))
edgelabels(text = signif(tree$edge.length, 2), frame = "none", adj = c(0.5, -0.4), cex = 0.35, col = "grey30")
if (any(show)) {
  nodelabels(text = bs[show], node = ((ntip + 1):(ntip + tree$Nnode))[show],
             frame = "none", adj = c(1.2, 1.2), cex = 0.35, col = "firebrick")
}
legend("topleft", legend = c(present, "branch length", "bootstrap (>= 50)"),
       col = c(species_colors[present], "grey30", "firebrick"),
       pch = c(rep(15, length(present)), NA, NA), lty = c(rep(NA, length(present)), NA, NA),
       text.col = c(rep("black", length(present)), "grey30", "firebrick"), bty = "n", cex = 0.8)
axisPhylo(side = 1, cex.axis = 0.6)
mtext("Substitutions per site (SNP alignment)", side = 1, line = 2.2, cex = 0.7)
dev.off()

# branch lengths as a table: parent node, child node or tip, length, support of the child node
edges <- data.frame(parent = tree$edge[, 1],
                    child = tree$edge[, 2],
                    child_label = ifelse(tree$edge[, 2] <= ntip, tree$tip.label[pmin(tree$edge[, 2], ntip)], ""),
                    branch_length = tree$edge.length,
                    child_support = ifelse(tree$edge[, 2] > ntip, bs[pmax(tree$edge[, 2] - ntip, 1)], NA))
write.table(edges, file.path(outdir, paste0("tree_", label, ".branch_lengths.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.tree(tree, file.path(outdir, paste0("tree_", label, ".rooted.nwk")))