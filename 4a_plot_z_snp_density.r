#!/usr/bin/env Rscript
# SNP distribution along one chromosome.
# Args: <sites.tsv> <chrom> <window_bp>
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[1] else "z_sites.tsv"
chrom  <- if (length(args) >= 2) args[2] else "CM059970.1"
win    <- if (length(args) >= 3) as.numeric(args[3]) else 1e5

d <- fread(infile, header = FALSE, col.names = c("chrom", "pos", "AC", "cM"),
           colClasses = c("character", "numeric", "character", "character"))
d <- d[chrom == ..chrom]
if (nrow(d) == 0) stop("No sites for ", chrom, " in ", infile)

# Variable in this sample = any non-zero allele count (handles multiallelic "1,0")
d[, polymorphic := !grepl("^0(,0)*$", AC) & AC != "."]
d[, cM := suppressWarnings(as.numeric(cM))]
d[, bin := floor(pos / win) * win]

# Count per window, keeping empty windows as 0 so gaps show as gaps
all_bins <- data.table(bin = seq(0, max(d$bin), by = win))
dens <- d[, .(all = .N, poly = sum(polymorphic)), by = bin]
dens <- merge(all_bins, dens, by = "bin", all.x = TRUE)
dens[is.na(all), `:=`(all = 0L, poly = 0L)]
dens_long <- melt(dens, id.vars = "bin", variable.name = "set", value.name = "n")
dens_long[, set := factor(set, levels = c("all", "poly"),
                          labels = c("All sites in VCF", "Variable in these samples (AC > 0)"))]

# ---- summary printed to the terminal ----
gaps <- dens[all == 0]
cat(sprintf("%s: %s sites, %s variable (%.1f%%), %.1f Mb spanned\n", chrom,
            format(nrow(d), big.mark = ","), format(sum(d$polymorphic), big.mark = ","),
            100 * mean(d$polymorphic), (max(d$pos) - min(d$pos)) / 1e6))
cat(sprintf("Windows (%g kb) with zero sites: %d of %d\n", win / 1e3, nrow(gaps), nrow(dens)))
if (nrow(gaps)) {
  g <- gaps[, .(start = bin)][, run := cumsum(c(1, diff(start) != win))]
  g <- g[, .(from_Mb = min(start) / 1e6, to_Mb = (max(start) + win) / 1e6), by = run]
  g[, length_Mb := to_Mb - from_Mb]
  cat("Longest empty stretches:\n"); print(head(g[order(-length_Mb)], 5))
}

# ---- plots ----
ink <- "#3d3d3a"; muted <- "#8a8984"
cols <- c("All sites in VCF" = "#2a78d6",
          "Variable in these samples (AC > 0)" = "#eb6834")

theme_set(theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e5e0", linewidth = 0.3),
        axis.text = element_text(colour = muted), axis.title = element_text(colour = ink),
        plot.title = element_text(colour = ink, face = "bold"),
        legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(colour = ink)))

p1 <- ggplot(dens_long, aes(bin / 1e6, n, colour = set)) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = cols) +
  labs(x = NULL, y = sprintf("SNPs per %g kb", win / 1e3),
       title = sprintf("SNP density along %s", chrom))

p2 <- ggplot(dens, aes(bin / 1e6, ifelse(all > 0, poly / all, NA))) +
  geom_line(linewidth = 0.5, colour = "#eb6834") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "Share variable")

plots <- list(p1, p2)
if (any(!is.na(d$cM))) {
  map <- d[!is.na(cM)][seq(1, .N, length.out = min(.N, 20000))]   # thin for plotting
  p3 <- ggplot(map, aes(pos / 1e6, cM)) +
    geom_line(linewidth = 0.5, colour = ink) +
    labs(x = "Position (Mb)", y = "Genetic position (cM)")
  plots <- c(plots, list(p3))
} else {
  plots[[2]] <- plots[[2]] + labs(x = "Position (Mb)")
}

out <- sprintf("snp_density_%s_%gkb.png", chrom, win / 1e3)
ggsave(out, wrap_plots(plots, ncol = 1, heights = c(2, 1, 1.3)[seq_along(plots)]),
       width = 11, height = 2.6 * length(plots) + 1, dpi = 300, bg = "#fcfcfb")
cat("Saved:", out, "\n")