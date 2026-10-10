#!/bin/bash
# -----------------------------------------------------------------------------
# Infer sex of each sample from mean read depth on the Z, normalized by the
# pseudoautosomal region (PAR, first PAR_END bp of the Z).
#   ratio = mean depth (Z-specific, > PAR_END + BUFFER) / mean depth (PAR)
#   males (ZZ) ~ 1.0, females (ZW) ~ 0.5
# Samples are assigned M/F by 2-cluster k-means on the ratio.
# Also writes raw mean depth in windows across the whole Z and plots it by
# inferred sex.
#
# Run:
#   nohup bash infer_sex_Zdepth.sh > logs/infer_sex_Zdepth.log 2>&1 &
# Optional overrides: PAR_END=50000000 BUFFER=20000000 MINQ=30 WIN=500000
# -----------------------------------------------------------------------------
set -euo pipefail

PROJDIR=/nas/bcch_for_danny/
VCF=${PROJDIR}/data/vcfs/Z_noW.all_ind.diploid.vcf.gz
OUTDIR=${PROJDIR}/results/sex_inference
PAR_END=${PAR_END:-50000000}   # end of PAR control
BUFFER=${BUFFER:-20000000}      # sites near the PAR boundary are skipped (sex calls only)
MINQ=${MINQ:-30}
WIN=${WIN:-500000}             # window size for depth along the Z
mkdir -p ${OUTDIR} logs

echo "[$(date)] Summarizing per-sample depth from ${VCF}"

# One pass over biallelic SNPs (missing GT/DP skipped):
#   stdout -> mean DP per sample in PAR and Z-specific regions
#   WINOUT -> summed DP and site count per sample per window (whole Z)
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

Rscript - ${OUTDIR} ${PAR_END} <<'EOF'
args <- commandArgs(TRUE); out <- args[1]; par_end <- as.numeric(args[2])
d <- read.delim(file.path(out, "Z_PAR_depth.tsv"))
ok <- !is.na(d$dp_ratio)

## --- Sex calls: 2-cluster k-means on Z/PAR depth ratio ---
set.seed(1)
km <- kmeans(d$dp_ratio[ok], centers = 2, nstart = 50)
fem <- which.min(km$centers)                      # lower-depth cluster = ZW
d$sex <- NA
d$sex[ok] <- ifelse(km$cluster == fem, "F", "M")

write.table(d, file.path(out, "sex_calls.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("Cluster centers (Z/PAR depth):", round(sort(km$centers), 3), "\n")
print(table(d$sex, useNA = "ifany"))

cols <- c(F = "#D55E00", M = "#0072B2")

pdf(file.path(out, "sex_calls.pdf"), width = 6, height = 4.5)
o <- order(d$dp_ratio)
plot(seq_along(o), d$dp_ratio[o], pch = 21, bg = cols[d$sex[o]],
     xlab = "Samples (ranked)", ylab = "Z-specific / PAR mean depth",
     ylim = c(0, max(1.2, d$dp_ratio, na.rm = TRUE)))
abline(h = c(0.5, 1), lty = 2, col = "grey50")
legend("topleft", c("F (ZW)", "M (ZZ)"), pt.bg = cols, pch = 21, bty = "n")
dev.off()

## --- Windowed raw mean depth along the Z, by inferred sex ---
w <- read.delim(file.path(out, "Z_window_depth.tsv"))
w <- merge(w, d[, c("sample", "sex")], by = "sample")
w <- w[!is.na(w$sex), ]
w$mid_mb <- (w$win_start + w$win_end) / 2 / 1e6
w <- w[order(w$sample, w$mid_mb), ]

# Per-sex median raw depth per window
med <- aggregate(mean_dp ~ sex + mid_mb, data = w, FUN = median)
write.table(w, file.path(out, "Z_window_depth_sex.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

ymax <- quantile(w$mean_dp, 0.995, na.rm = TRUE) * 1.05   # trim extreme windows
pdf(file.path(out, "Z_window_depth.pdf"), width = 11, height = 4.5)
plot(NA, xlim = range(w$mid_mb), ylim = c(0, ymax),
     xlab = "Position on Z (Mb)", ylab = "Mean depth (DP) per window")
rect(0, -1, par_end / 1e6, ymax * 2, col = "grey95", border = NA)
for (s in unique(w$sample)) {
  ws <- w[w$sample == s, ]
  lines(ws$mid_mb, ws$mean_dp, col = adjustcolor(cols[ws$sex[1]], 0.25), lwd = 0.6)
}
for (sx in c("F", "M")) {
  m <- med[med$sex == sx, ]
  lines(m$mid_mb, m$mean_dp, col = cols[sx], lwd = 2.5)
}
abline(v = par_end / 1e6, lty = 3)
legend("topright", c("F (ZW) median", "M (ZZ) median", "PAR control"),
       col = c(cols, "grey85"), lwd = c(2.5, 2.5, 8), bty = "n", cex = 0.8)
box()
dev.off()
EOF

echo "[$(date)] Done: ${OUTDIR}/sex_calls.tsv, sex_calls.pdf, Z_window_depth.pdf"
