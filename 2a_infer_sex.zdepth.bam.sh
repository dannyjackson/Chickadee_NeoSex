#!/bin/bash
# -----------------------------------------------------------------------------
# Infer sex of each sample in one species from read depth in the BAMs:
#   ratio = mean depth on the ancestral Z (ANC_START-ANC_END)
#         / mean depth across the autosomes        (same BAM, same filters)
#   males (ZZ) ~ 1, females (ZW) ~ 0.5
# Samples are assigned M/F by 2-cluster k-means on the ratio (within species).
# Also plots raw mean depth in windows along the whole Z for every individual.
#
# Samples come from reference_lists/<species>_samples.txt (first column); their
# BAMs are looked up in BAMLIST by read-group SM tag or by file name.
#
# Usage:
#   bash 2a_infer_sex.zdepth.sh <species>      # e.g. moch
#
# Run all four species (from 2_Infer_Sex.md):
#   PROJDIR=/nas/bcch_for_danny
#   SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex
#   for sp in bcch cach moch boch; do
#     mkdir -p ${PROJDIR}/analyses/sex_inference/${sp}
#     nohup bash ${SCRIPTDIR}/2a_infer_sex.zdepth.sh ${sp} \
#       > ${PROJDIR}/analyses/sex_inference/${sp}/infer_sex.log 2>&1 &
#   done
#
# Optional env overrides:
#   BAMLIST=...  ANC_START=65501052 ANC_END=146584261  ZCHR=CM059970.1
#   WCHR=CM059969.1  MINMQ=20  WIN=500000  NJOBS=4 (BAMs processed in parallel)
# Per-BAM depth results are cached in analyses/sex_inference/bam_depth/ and
# reused on reruns; delete that folder to recompute (e.g. after changing MINMQ).
# -----------------------------------------------------------------------------
set -euo pipefail

SP=${1:?Usage: $0 <species>   (e.g. moch)}

PROJDIR=/nas/bcch_for_danny
SAMPLES=${SAMPLES:-${PROJDIR}/reference_lists/${SP}_samples.txt}
BAMLIST=${BAMLIST:-${PROJDIR}/reference_lists/sorted_bams.noW_genome.txt}
OUTBASE=${PROJDIR}/analyses/sex_inference
OUTDIR=${OUTBASE}/${SP}
CACHE=${OUTBASE}/bam_depth
ZCHR=${ZCHR:-CM059970.1}           # Z
WCHR=${WCHR:-CM059969.1}           # W (excluded from autosomes)
ANC_START=${ANC_START:-65501052}    # ancestral Z
ANC_END=${ANC_END:-146584261}
MINMQ=${MINMQ:-20}                  # min mapping quality (both Z and autosomes)
WIN=${WIN:-500000}
NJOBS=${NJOBS:-4}
mkdir -p ${OUTDIR} ${CACHE}

echo "[$(date)] Species: ${SP}"
echo "  Samples from: ${SAMPLES}"
echo "  BAM list:     ${BAMLIST}"
echo "  Z: ${ZCHR}; ancestral Z ${ANC_START}-${ANC_END}; MAPQ>=${MINMQ}; window ${WIN}"

# --- 1. Match this species' samples to BAMs (by SM tag or file-name stem) ---
[[ -s ${SAMPLES} ]] || { echo "ERROR: sample list ${SAMPLES} missing or empty"; exit 1; }
awk 'NF && $1 !~ /^#/ {print $1}' ${SAMPLES} | tr -d '\r' > ${OUTDIR}/samples.txt
: > ${OUTDIR}/bams.tsv
while read -r bam; do
  [[ -z ${bam} ]] && continue
  stem=$(basename ${bam}); stem=${stem%%.*}
  sm=$(samtools view -H ${bam} | awk -F'\t' '/^@RG/{for(i=2;i<=NF;i++) if($i~/^SM:/){print substr($i,4); exit}}')
  hit=$(awk -v a="${sm}" -v b="${stem}" -v c="${bam}" '$1==a || $1==b || $1==c {print $1; exit}' ${OUTDIR}/samples.txt)
  [[ -n ${hit} ]] && printf "%s\t%s\n" "${hit}" "${bam}" >> ${OUTDIR}/bams.tsv
done < ${BAMLIST}

n_s=$(wc -l < ${OUTDIR}/samples.txt); n_b=$(wc -l < ${OUTDIR}/bams.tsv)
echo "  ${n_b} of ${n_s} samples matched to a BAM"
if (( n_b < n_s )); then
  echo "  WARNING: no BAM found for:"; cut -f1 ${OUTDIR}/bams.tsv | grep -vxFf - ${OUTDIR}/samples.txt | sed 's/^/    /' || true
fi
(( n_b > 0 )) || { echo "ERROR: no BAMs matched"; exit 1; }

# --- 2. Depth per BAM (cached): autosomes, ancestral Z, Z windows ---
#   autosomes: samtools coverage over all chromosomes; length-weighted mean over
#              CM* contigs other than Z and W
#   Z:         samtools depth -a over the Z; mean over ancestral Z + windows
depth_one () {
  local sample=$1 bam=$2 tag
  tag=$(basename ${bam}); tag=${tag%%.*}
  local f=${CACHE}/${tag}.MQ${MINMQ}.W${WIN}
  if [[ -s ${f}.summary.tsv && -s ${f}.windows.tsv ]]; then
    echo "  [cached] ${sample}"; return
  fi
  echo "  [$(date +%T)] ${sample}: ${bam}"
  auto=$(samtools coverage --min-MQ ${MINMQ} ${bam} \
    | awk -v z=${ZCHR} -v w=${WCHR} -F'\t' 'NR>1 && $1 ~ /^CM/ && $1 != z && $1 != w {
          L = $3 - $2 + 1; s += $7 * L; t += L }
        END{ printf "%s\t%s", t, (t ? s / t : "NA") }')
  samtools depth -a -Q ${MINMQ} -r ${ZCHR} ${bam} \
  | awk -v a0=${ANC_START} -v a1=${ANC_END} -v win=${WIN} -v wout=${f}.windows.tmp \
        -v auto="${auto}" -v sm="${sample}" 'BEGIN{OFS="\t"}
      { k = int(($2 - 1) / win); n[k]++; d[k] += $3; if (k > mk) mk = k
        if ($2 >= a0 && $2 <= a1) { na++; da += $3 } }
      END{
        split(auto, A, "\t")
        print "sample", "auto_bp", "dp_auto", "Zanc_bp", "dp_Zanc"
        print sm, A[1], A[2], na + 0, (na ? da / na : "NA")
        print "sample", "win_start", "win_end", "n_bp", "mean_dp" > wout
        for (k = 0; k <= mk; k++) if (n[k])
          print sm, k * win + 1, (k + 1) * win, n[k], d[k] / n[k] > wout
      }' > ${f}.summary.tmp
  mv ${f}.windows.tmp ${f}.windows.tsv
  mv ${f}.summary.tmp ${f}.summary.tsv
}
export -f depth_one
export CACHE MINMQ WIN ZCHR WCHR ANC_START ANC_END

echo "[$(date)] Computing depth from BAMs (${NJOBS} at a time)"
xargs -a ${OUTDIR}/bams.tsv -P ${NJOBS} -L 1 bash -c 'depth_one "$0" "$1"'

# Gather this species' results
{ echo -e "sample\tauto_bp\tdp_auto\tZanc_bp\tdp_Zanc"
  cut -f2 ${OUTDIR}/bams.tsv | while read -r bam; do
    t=$(basename ${bam}); tail -n +2 ${CACHE}/${t%%.*}.MQ${MINMQ}.W${WIN}.summary.tsv; done
} > ${OUTDIR}/bam_depth.tsv
{ echo -e "sample\twin_start\twin_end\tn_bp\tmean_dp"
  cut -f2 ${OUTDIR}/bams.tsv | while read -r bam; do
    t=$(basename ${bam}); tail -n +2 ${CACHE}/${t%%.*}.MQ${MINMQ}.W${WIN}.windows.tsv; done
} > ${OUTDIR}/Z_window_depth.tsv

# --- 3. Ratio, k-means sex calls, plots ---
echo "[$(date)] Clustering and plotting"
Rscript - ${OUTDIR} ${SP} ${ANC_START} ${ANC_END} <<'EOF'
args <- commandArgs(TRUE)
out <- args[1]; sp <- toupper(args[2]); anc <- as.numeric(args[3:4]) / 1e6

z <- read.delim(file.path(out, "bam_depth.tsv"))
z$dp_ratio <- z$dp_Zanc / z$dp_auto

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
write.table(z[, c("sample", "species", "dp_auto", "dp_Zanc", "dp_ratio", "sex")],
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
     xlab = "Position on Z (Mb)", ylab = "Mean depth per window (BAM)",
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
