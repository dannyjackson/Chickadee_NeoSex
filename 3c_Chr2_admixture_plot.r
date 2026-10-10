#!/usr/bin/env Rscript
# B1_admixture_plot.r
# Plots ADMIXTURE ancestry proportions for each K, samples grouped by species
# usage: Rscript B1_admixture_plot.r <admix_dir> <admix_prefix> <sample_info.txt> <kmin> <kmax> <label>
# outputs (in admix_dir):
#   admix_<label>_allK.pdf       all K stacked, one row per K
#   admix_<label>_K<K>.pdf       one plot per K
#   admix_<label>_cv_error.pdf   cross-validation error by K
#   admix_<label>_Q_allK.txt     ancestry proportions with sample and species

required_packages <- c("ggplot2", "dplyr", "tidyr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Missing R package: ", pkg, ". Install with: micromamba install -c conda-forge r-", pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
admixdir <- args[1]
prefix <- args[2]
sampleinfo_file <- args[3]
kmin <- as.integer(args[4])
kmax <- as.integer(args[5])
label <- args[6]

# sample order in the Q files matches the .fam
fam <- read.table(file.path(admixdir, paste0(prefix, ".fam")), header = FALSE)
samples <- fam$V2

labs <- read.table(sampleinfo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
species_order <- c("BCCH", "CACH", "MOCH", "BOCH", "Unknown")

info <- data.frame(Sample = samples, stringsAsFactors = FALSE) %>%
  left_join(labs, by = "Sample") %>%
  mutate(Species = ifelse(is.na(Species), "Unknown", Species))

# order samples by species, then by ancestry at the highest K so similar individuals sit together
qmax <- read.table(file.path(admixdir, paste0(prefix, ".", kmax, ".Q")), header = FALSE)
info$main_cluster <- max.col(qmax)
info$main_prop <- apply(qmax, 1, max)
info <- info %>%
  mutate(Species = factor(Species, levels = intersect(species_order, unique(Species)))) %>%
  arrange(Species, main_cluster, desc(main_prop))
sample_order <- info$Sample

# read all K into long format
qlong <- bind_rows(lapply(kmin:kmax, function(k) {
  q <- read.table(file.path(admixdir, paste0(prefix, ".", k, ".Q")), header = FALSE)
  colnames(q) <- paste0("Cluster", seq_len(k))
  q$Sample <- samples
  q %>%
    pivot_longer(starts_with("Cluster"), names_to = "Cluster", values_to = "Ancestry") %>%
    mutate(K = paste0("K=", k))
})) %>%
  left_join(info %>% dplyr::select(Sample, Species), by = "Sample") %>%
  mutate(Sample = factor(Sample, levels = sample_order),
         K = factor(K, levels = paste0("K=", kmin:kmax)),
         Cluster = factor(Cluster, levels = paste0("Cluster", seq_len(kmax))))

write.table(qlong, file.path(admixdir, paste0("admix_", label, "_Q_allK.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cluster_colors <- setNames(hcl.colors(kmax, palette = "Set 2"), paste0("Cluster", seq_len(kmax)))
n_samples <- length(samples)
plot_width <- max(8, n_samples * 0.12)

admix_plot <- function(df) {
  ggplot(df, aes(x = Sample, y = Ancestry, fill = Cluster)) +
    geom_col(width = 1) +
    facet_grid(K ~ Species, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = cluster_colors) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
    labs(x = NULL, y = "Ancestry", title = paste0("ADMIXTURE, ", label)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
          panel.spacing.x = unit(0.15, "lines"),
          strip.background = element_blank(),
          legend.position = "none")
}

# all K together
pdf(file = file.path(admixdir, paste0("admix_", label, "_allK.pdf")),
    width = plot_width, height = 1.6 * (kmax - kmin + 1) + 1.5, useDingbats = FALSE)
print(admix_plot(qlong))
dev.off()

# one plot per K
for (k in kmin:kmax) {
  pdf(file = file.path(admixdir, paste0("admix_", label, "_K", k, ".pdf")),
      width = plot_width, height = 3, useDingbats = FALSE)
  print(admix_plot(filter(qlong, K == paste0("K=", k))))
  dev.off()
}

# cross-validation error
cv_file <- file.path(admixdir, paste0(prefix, ".cv_error.txt"))
if (file.exists(cv_file) && file.info(cv_file)$size > 0) {
  cv <- read.table(cv_file, header = FALSE, col.names = c("K", "CV_error"))
  pdf(file = file.path(admixdir, paste0("admix_", label, "_cv_error.pdf")), width = 5, height = 4, useDingbats = FALSE)
  print(ggplot(cv, aes(x = K, y = CV_error)) +
          geom_line() + geom_point(size = 2) +
          scale_x_continuous(breaks = cv$K) +
          labs(x = "K", y = "Cross-validation error", title = paste0("ADMIXTURE CV error, ", label)) +
          theme_classic())
  dev.off()
  cat("Best supported K (lowest CV error):", cv$K[which.min(cv$CV_error)], "\n")
}
