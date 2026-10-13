#!/bin/bash
# -----------------------------------------------------------------------------
# Infer sex of each sample in one species from read depth on the ancestral Z,
# normalized by the same sample's mean depth across the autosomes.
#   ratio = mean DP (ancestral Z: ANC_START-ANC_END) / mean DP (autosomes)
#   males (ZZ) ~ 1, females (ZW) ~ 0.5
# Samples are assigned M/F by 2-cluster k-means on the ratio (within species).
# Also plots raw mean DP in windows along the whole Z for every individual.
#
# Usage:
#   bash 2a_infer_sex.zdepth.sh <Z_VCF> <AUTOSOME_VCF> [SPECIES_LABEL]
#   SPECIES_LABEL defaults to the 2nd dot-field of the Z VCF name
#   (Z_noW.bcch.diploid.vcf.gz -> bcch).
#
# Run all four species (from 2_Infer_Sex.md):
#   PROJDIR=/nas/bcch_for_danny
#   for sp in bcch cach moch boch; do
#     mkdir -p ${PROJDIR}/analyses/sex_inference/${sp}
#     nohup bash 2a_infer_sex.zdepth.sh \
#       ${PROJDIR}/data/vcfs/Z_noW.${sp}.diploid.vcf.gz \
#       ${PROJDIR}/data/vcfs/autosomes.all_ind.diploid.vcf.gz \
#       > ${PROJDIR}/analyses/sex_inference/${sp}/infer_sex.log 2>&1 &
#   done
#   Autosomal depth is computed once (all samples) and cached; parallel runs
#   wait for each other on a lock instead of recomputing it.
#
# Optional env overrides:
#   ANC_START=65501052 ANC_END=146584261 MINQ=30 WIN=500000
#   AUTO_REGIONS=chr1,chr2   # restrict autosomal depth to some chromosomes
#                            # (faster; needs an indexed VCF)
# -----------------------------------------------------------------------------
set -euo pipefail

Z_VCF=${1:?Usage: $0 <Z_VCF> <AUTOSOME_VCF> [SPECIES_LABEL]}
AUTO_VCF=${2:?Usage: $0 <Z_VCF> <AUTOSOME_VCF> [SPECIES_LABEL]}
SP=${3:-$(basename ${Z_VCF} | cut -d. -f2)}

PROJDIR=/nas/bcch_for_danny
OUTBASE=${PROJDIR}/analyses/sex_inference
OUTDIR=${OUTBASE}/${SP}
ANC_START=${ANC_START:-65501052}   # ancestral Z
ANC_END=${ANC_END:-146584261}
MINQ=${MINQ:-30}
WIN=${WIN:-500000}
AUTO_REGIONS=${AUTO_REGIONS:-}
mkdir -p ${OUTDIR}

echo "[$(date)] Species: ${SP}"
echo "  Z VCF:        ${Z_VCF}"
echo "  Autosome VCF: ${AUTO_VCF}"
echo "  Ancestral Z:  ${ANC_START}-${ANC_END}; window: ${WIN}; QUAL>=${MINQ}"

# --- 1. Autosomal mean DP per sample (all samples, cached, shared by runs) ---
REGTAG=$( [[ -n ${AUTO_REGIONS} ]] && echo ".$(echo ${AUTO_REGIONS} | tr ',' '_')" || true )
AUTO_TSV=${OUTBASE}/autosome_depth.$(basename ${AUTO_VCF} .vcf.gz).Q${MINQ}${REGTAG}.tsv
(
  flock 9
  if [[ -s ${AUTO_TSV} ]]; then
    echo "[$(date)] Using cached autosomal depth: ${AUTO_TSV}"
  else
    echo "[$(date)] Computing autosomal depth (slow, done once): ${AUTO_TSV}"
    REG=""; [[ -n ${AUTO_REGIONS} ]] && REG="-r ${AUTO_REGIONS}"
    bcftools view ${REG} -m2 -M2 -v snps -i "QUAL>=${MINQ}" ${AUTO_VCF} \
    | bcftools query -f '[%SAMPLE\t%GT\t%DP\n]' \
    | awk 'BEGIN{OFS="\t"}
        $2 !~ /\./ && $3 != "." { n[$1]++; d[$1] += $3 }
        END{ print "sample", "n_auto", "dp_auto"
             for (s in n) print s, n[s], d[s] / n[s] }' \
    > ${AUTO_TSV}.tmp
    mv ${AUTO_TSV}.tmp ${AUTO_TSV}
  fi
) 9> ${AUTO_TSV}.lock

# --- 2. Z: ancestral-Z mean DP per sample + windowed DP along the whole Z ---
echo "[$(date)] Summarizing Z depth"
bcftools view -m2 -M2 -v snps -i "QUAL>=${MINQ}" ${Z_VCF} \
| bcftools query -f '%POS[\t%SAMPLE\t%GT\t%DP]\n' \
| awk -v a0=${ANC_START} -v a1=${ANC_END} -v win=${WIN} \
      -v winout=${OUTDIR}/Z_window_depth.tsv 'BEGIN{OFS="\t"}
    {
      w = int(($1 - 1) / win); anc = ($1 >= a0 && $1 <= a1)
      for (i = 2; i <= NF; i += 3) {
        if ($(i+1) ~ /\./ || $(i+2) == ".") continue
        s[$i] = 1
        wn[$i, w]++; wd[$i, w] += $(i+2); if (w > maxw) maxw = w
        if (anc) { n[$i]++; d[$i] += $(i+2) }
      }
    }
    END{
      print "sample", "n_Zanc", "dp_Zanc"
      for (x in s) print x, n[x] + 0, (n[x] ? d[x] / n[x] : "NA")
      print "sample", "win_start", "win_end", "n_sites", "mean_dp" > winout
      for (x in s) for (k = 0; k <= maxw; k++) if (wn[x, k])
        print x, k * win + 1, (k + 1) * win, wn[x, k], wd[x, k] / wn[x, k] > winout
    }' \
> ${OUTDIR}/Z_depth.tsv

# --- 3. Ratio, k-means sex calls, plots ---
echo "[$(date)] Clustering and plotting"
Rscript - ${OUTDIR} ${AUTO_TSV} ${SP} ${ANC_START} ${ANC_END} <<'EOF'
args <- commandArgs(TRUE)
out <- args[1]; auto_f <- args[2]; sp <- toupper(args[3])
anc <- as.numeric(args[4:5]) / 1e6

z    <- read.delim(file.path(out, "Z_depth.tsv"))
auto <- read.delim(auto_f)

# Z and autosome VCFs come from different BAM sets; match names after
# stripping directory and extension (/path/S1.sorted.bam -> S1).
norm <- function(x) sub("\\..*$", "", basename(x))
auto$id <- norm(auto$sample)
m <- match(norm(z$sample), auto$id)
z$n_auto  <- auto$n_auto[m]
z$dp_auto <- auto$dp_auto[m]
z$dp_ratio <- z$dp_Zanc / z$dp_auto
if (any(is.na(m)))
  cat("WARNING: no autosomal depth for:", z$sample[is.na(m)], "\n")

ok <- !is.na(z$dp_ratio)
z$sex <- NA
if (sum(ok) >= 2 && length(unique(z$dp_ratio[ok])) >= 2) {
  set.seed(1)
  km <- kmeans(z$dp_ratio[ok], centers = 2, nstart = 50)
  z$sex[ok] <- ifelse(km$cluster == which.min(km$centers), "F", "M")
  cut <- mean(km$centers)
  cat("Cluster centers (Zanc/autosome depth):", round(sort(km$centers), 3),
      " boundary:", round(cut, 3), "\n")
  # k-means always makes two groups; warn if both sit on the same side of 0.75
  if (all(km$centers < 0.75) || all(km$centers > 0.75))
    cat("WARNING: both clusters on one side of 0.75 - may be a single-sex sample\n")
} else {
  cut <- NA; cat("WARNING: too few samples with a ratio to cluster\n")
}

z <- z[order(z$dp_ratio), ]
z$species <- sp
write.table(z[, c("sample", "species", "n_auto", "dp_auto", "n_Zanc", "dp_Zanc",
                  "dp_ratio", "sex")],
            file.path(out, "sex_calls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table(z$sex, useNA = "ifany"))

cols <- c(F = "#D55E00", M = "#0072B2")
nF <- sum(z$sex == "F", na.rm = TRUE); nM <- sum(z$sex == "M", na.rm = TRUE)

## Plot 1: ratio per sample
zs <- z[ok, ]
pdf(file.path(out, "sex_calls.pdf"), width = max(5, 2 + 0.18 * nrow(zs)), height = 4.5)
par(mar = c(7, 4.5, 3, 1))
plot(seq_len(nrow(zs)), zs$dp_ratio, pch = 21, bg = cols[zs$sex], cex = 1.2,
     ylim = c(0, max(1.3, zs$dp_ratio)), xaxt = "n", xlab = "",
     ylab = "Ancestral Z / autosome mean depth",
     main = sprintf("%s (F=%d, M=%d)", sp, nF, nM))
axis(1, at = seq_len(nrow(zs)), labels = zs$sample, las = 2, cex.axis = 0.6)
abline(h = c(0.5, 1), lty = 2, col = "grey60")
if (!is.na(cut)) abline(h = cut, col = "red", lty = 3)
legend("bottomright", c("F (ZW)", "M (ZZ)", "k-means boundary"),
       pt.bg = c(cols, NA), pch = c(21, 21, NA), lty = c(NA, NA, 3),
       col = c("black", "black", "red"), bty = "n", cex = 0.8)
dev.off()

## Plot 2: raw windowed depth along the Z, every individual
w <- read.delim(file.path(out, "Z_window_depth.tsv"))
w <- merge(w, z[, c("sample", "sex")], by = "sample")
w$sex[is.na(w$sex)] <- "unknown"
w$mid_mb <- (w$win_start + w$win_end) / 2 / 1e6
w <- w[order(w$sample, w$mid_mb), ]
write.table(w, file.path(out, "Z_window_depth_sex.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

lcols <- c(cols, unknown = "grey40")
ymax <- quantile(w$mean_dp, 0.995, na.rm = TRUE) * 1.05   # trim spikes
pdf(file.path(out, "Z_window_depth.pdf"), width = 11, height = 4.5)
par(mar = c(4.5, 4.5, 3, 1))
plot(NA, xlim = range(w$mid_mb), ylim = c(0, ymax),
     xlab = "Position on Z (Mb)", ylab = "Mean DP per window",
     main = sprintf("%s (F=%d, M=%d)", sp, nF, nM))
rect(anc[1], -1, anc[2], ymax * 2, col = "grey93", border = NA)
# females drawn last so they are never hidden under males
for (sx in c("unknown", "M", "F")) for (s in unique(w$sample[w$sex == sx])) {
  wi <- w[w$sample == s, ]
  lines(wi$mid_mb, wi$mean_dp, col = adjustcolor(lcols[sx], 0.6), lwd = 0.7)
}
box()
legend("topleft", c("F (ZW)", "M (ZZ)", "ancestral Z (used for ratio)"),
       col = c(cols, "grey85"), lwd = c(1.5, 1.5, 8), bty = "n", cex = 0.8, ncol = 3)
dev.off()
EOF

echo "[$(date)] Done: ${OUTDIR}/sex_calls.tsv, sex_calls.pdf, Z_window_depth.pdf"
