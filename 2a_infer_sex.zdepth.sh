#!/bin/bash
# -----------------------------------------------------------------------------
# Infer sex of each sample from mean read depth on the Z, normalized by the
# pseudoautosomal region (PAR, first PAR_END bp of the Z).
#   PAR        : POS <= PAR_END                      (diploid in both sexes)
#   gap        : PAR_END < POS <= PAR_END + BUFFER   (skipped for sex calls)
#   Z-specific : POS >  PAR_END + BUFFER
#   ratio = mean DP (Z-specific) / mean DP (PAR); males (ZZ) ~1, females (ZW) ~0.5
# Samples are assigned M/F by 2-cluster k-means on the ratio (all species pooled).
# Plots are split by species (from reference_lists/<sp>_samples.txt).
#
# Run (from anywhere):
#   nohup bash 2a_infer_sex.zdepth.sh > /nas/bcch_for_danny/analyses/sex_inference/infer_sex_Zdepth.log 2>&1 &
# Optional overrides: PAR_END=50000000 BUFFER=20000000 MINQ=30 WIN=500000
# -----------------------------------------------------------------------------
set -euo pipefail

PROJDIR=/nas/bcch_for_danny
VCF=${PROJDIR}/data/vcfs/Z_noW.all_ind.diploid.vcf.gz
OUTDIR=${PROJDIR}/analyses/sex_inference
REFLISTS=${PROJDIR}/reference_lists
SPECIES="bcch cach moch boch"     # expects ${REFLISTS}/<sp>_samples.txt
PAR_END=${PAR_END:-50000000}      # end of PAR control
BUFFER=${BUFFER:-20000000}        # gap after the PAR excluded from sex calls
MINQ=${MINQ:-30}
WIN=${WIN:-500000}                # window size for depth along the Z
mkdir -p ${OUTDIR}

echo "[$(date)] Summarizing per-sample depth from ${VCF}"
echo "  PAR: 1-${PAR_END}; gap: to $((PAR_END + BUFFER)); window: ${WIN}"

# Species key: sample <tab> species, from the per-species sample lists
: > ${OUTDIR}/species_key.tsv
for sp in ${SPECIES}; do
  f=${REFLISTS}/${sp}_samples.txt
  if [[ -s ${f} ]]; then
    awk -v sp=${sp} 'NF{print $1"\t"toupper(sp)}' ${f} >> ${OUTDIR}/species_key.tsv
  else
    echo "  WARNING: ${f} missing or empty"
  fi
done

# One pass over biallelic SNPs (missing GT/DP skipped):
#   stdout -> mean DP per sample in PAR and Z-specific regions
#   winout -> mean DP and site count per sample per window (whole Z)
bcftools view -m2 -M2 -v snps -i "QUAL>=${MINQ}" ${VCF} \
| bcftools query -f '%POS[\t%SAMPLE\t%GT\t%DP]\n' \
| awk -v par=${PAR_END} -v buf=${BUFFER} -v win=${WIN} \
      -v winout=${OUTDIR}/Z_window_depth.tsv 'BEGIN{OFS="\t"}
    {
      w = int(($1 - 1) / win)
      if ($1 <= par) r = "PAR"
      else if ($1 > par + buf) r = "Z"
      else r = ""
      for (i = 2; i <= NF; i += 3) {
        if ($(i+1) ~ /\./ || $(i+2) == ".") continue
        s[$i] = 1
        wn[$i, w]++; wd[$i, w] += $(i+2); if (w > maxw) maxw = w
        if (r != "") { n[$i, r]++; d[$i, r] += $(i+2) }
      }
    }
    END{
      print "sample", "n_PAR", "dp_PAR", "n_Z", "dp_Z", "dp_ratio"
      for (x in s) {
        p = n[x, "PAR"] ? d[x, "PAR"] / n[x, "PAR"] : "NA"
        z = n[x, "Z"]   ? d[x, "Z"]   / n[x, "Z"]   : "NA"
        print x, n[x, "PAR"] + 0, p, n[x, "Z"] + 0, z, (p != "NA" && z != "NA" && p > 0) ? z / p : "NA"
      }
      print "sample", "win_start", "win_end", "n_sites", "mean_dp" > winout
      for (x in s) for (k = 0; k <= maxw; k++) if (wn[x, k])
        print x, k * win + 1, (k + 1) * win, wn[x, k], wd[x, k] / wn[x, k] > winout
    }' \
| { IFS= read -r hdr; echo "${hdr}"; sort -k1,1; } > ${OUTDIR}/Z_PAR_depth.tsv

echo "[$(date)] Clustering and plotting"

Rscript - ${OUTDIR} ${PAR_END} ${BUFFER} <<'EOF'
args <- commandArgs(TRUE)
out <- args[1]; par_end <- as.numeric(args[2]); buffer <- as.numeric(args[3])
d <- read.delim(file.path(out, "Z_PAR_depth.tsv"))

## --- Species assignment ---
# Match VCF sample names to sample lists; names are compared after stripping any
# directory and extension (e.g. /path/S1.sorted.bam -> S1) so either form works.
norm <- function(x) sub("\\..*$", "", basename(x))
kf <- file.path(out, "species_key.tsv")
key <- if (file.size(kf) > 0) read.delim(kf, header = FALSE, col.names = c("id", "species"))
       else data.frame(id = character(), species = character())
key$id <- norm(key$id)
dups <- key$id[duplicated(key$id)]
if (length(dups)) cat("WARNING: samples in >1 species list:", unique(dups), "\n")
key <- key[!duplicated(key$id), ]
d$species <- key$species[match(norm(d$sample), key$id)]
d$species[is.na(d$species)] <- "unassigned"
if (any(d$species == "unassigned"))
  cat("WARNING: no species for:", d$sample[d$species == "unassigned"], "\n")
miss <- setdiff(key$id, norm(d$sample))
if (length(miss)) cat("Note: in sample lists but not in VCF:", miss, "\n")

## --- Sex calls: 2-cluster k-means on Z/PAR depth ratio (species pooled) ---
ok <- !is.na(d$dp_ratio)
set.seed(1)
km <- kmeans(d$dp_ratio[ok], centers = 2, nstart = 50)
fem <- which.min(km$centers)                      # lower-depth cluster = ZW
d$sex <- NA
d$sex[ok] <- ifelse(km$cluster == fem, "F", "M")
cut <- mean(km$centers)                           # boundary between clusters

d <- d[order(d$species, d$dp_ratio), ]
write.table(d, file.path(out, "sex_calls.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("Cluster centers (Z/PAR depth):", round(sort(km$centers), 3),
    " boundary:", round(cut, 3), "\n")
print(table(d$species, d$sex, useNA = "ifany"))

cols <- c(F = "#D55E00", M = "#0072B2")
sp_levels <- unique(d$species)
nsp <- length(sp_levels)

## --- Plot 1: Z/PAR ratio per sample, one panel per species ---
pdf(file.path(out, "sex_calls.pdf"), width = 3 + 2.2 * nsp, height = 4.5)
par(mfrow = c(1, nsp), mar = c(4, 4, 2.5, 0.5), oma = c(0, 0, 0, 0))
ylim <- c(0, max(1.3, d$dp_ratio, na.rm = TRUE))
for (sp in sp_levels) {
  ds <- d[d$species == sp & !is.na(d$dp_ratio), ]
  if (!nrow(ds)) { plot.new(); title(main = paste(sp, "(no ratio)")); next }
  plot(seq_len(nrow(ds)), ds$dp_ratio, pch = 21, bg = cols[ds$sex], cex = 1.2,
       ylim = ylim, xlim = c(0.5, max(nrow(ds), 1) + 0.5),
       main = sprintf("%s (F=%d, M=%d)", sp, sum(ds$sex == "F"), sum(ds$sex == "M")),
       xlab = "Samples (ranked)",
       ylab = if (sp == sp_levels[1]) "Z-specific / PAR mean depth" else "")
  abline(h = c(0.5, 1), lty = 2, col = "grey60")
  abline(h = cut, col = "red", lty = 3)
}
legend("bottomright", c("F (ZW)", "M (ZZ)", "k-means boundary"),
       pt.bg = c(cols, NA), pch = c(21, 21, NA), lty = c(NA, NA, 3),
       col = c("black", "black", "red"), bty = "n", cex = 0.8)
dev.off()

## --- Plot 2: raw windowed depth for every individual, one panel per species ---
w <- read.delim(file.path(out, "Z_window_depth.tsv"))
w <- merge(w, d[, c("sample", "species", "sex")], by = "sample")
w <- w[!is.na(w$sex), ]
w$mid_mb <- (w$win_start + w$win_end) / 2 / 1e6
w <- w[order(w$sample, w$mid_mb), ]
write.table(w, file.path(out, "Z_window_depth_sex.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

xr <- range(w$mid_mb)
pdf(file.path(out, "Z_window_depth.pdf"), width = 11, height = 2.6 * nsp + 0.6)
par(mfrow = c(nsp, 1), mar = c(2, 4.5, 2, 1), oma = c(2.5, 0, 0, 0))
for (sp in sp_levels) {
  ws <- w[w$species == sp, ]
  if (!nrow(ws)) { plot.new(); title(main = paste(sp, "(no sexed samples)")); next }
  ymax <- quantile(ws$mean_dp, 0.995, na.rm = TRUE) * 1.05   # trim spikes
  plot(NA, xlim = xr, ylim = c(0, ymax), xlab = "", ylab = "Mean DP per window",
       main = sprintf("%s (F=%d, M=%d)", sp,
                      length(unique(ws$sample[ws$sex == "F"])),
                      length(unique(ws$sample[ws$sex == "M"]))))
  rect(0, -1, par_end / 1e6, ymax * 2, col = "grey92", border = NA)
  rect(par_end / 1e6, -1, (par_end + buffer) / 1e6, ymax * 2,
       col = "grey97", border = NA)
  # draw females last so they are never hidden under the (usually more) males
  for (sx in c("M", "F")) for (s in unique(ws$sample[ws$sex == sx])) {
    wi <- ws[ws$sample == s, ]
    lines(wi$mid_mb, wi$mean_dp, col = adjustcolor(cols[sx], 0.6), lwd = 0.7)
  }
  box()
  if (sp == sp_levels[1])
    legend("topright", c("F (ZW)", "M (ZZ)", "PAR control", "excluded gap"),
           col = c(cols, "grey85", "grey95"), lwd = c(1.5, 1.5, 8, 8),
           bty = "n", cex = 0.8, ncol = 4)
}
mtext("Position on Z (Mb)", side = 1, outer = TRUE, line = 1)
dev.off()
EOF

echo "[$(date)] Done: ${OUTDIR}/sex_calls.tsv, sex_calls.pdf, Z_window_depth.pdf"