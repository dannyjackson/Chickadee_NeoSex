#!/bin/bash
# -----------------------------------------------------------------------------
# Infer sex of each sample from mean read depth on the Z, normalized by the
# pseudoautosomal region (PAR, first PAR_END bp of the Z).
#   ratio = mean depth (Z-specific, > PAR_END + BUFFER) / mean depth (PAR)
#   males (ZZ) ~ 1.0, females (ZW) ~ 0.5
# Samples are assigned M/F by 2-cluster k-means on the ratio.
#
# Run:
#   nohup bash infer_sex_Zdepth.sh > logs/infer_sex_Zdepth.log 2>&1 &
# Optional overrides: PAR_END=50000000 BUFFER=1000000 MINQ=30
# -----------------------------------------------------------------------------
set -euo pipefail

PROJDIR=/nas/bcch_for_danny/
SCRIPTDIR=${PROJDIR}/scripts/Chickadee_NeoSex
VCF=${PROJDIR}/data/vcfs/Z.all_ind.diploid.vcf.gz
OUTDIR=${PROJDIR}/results/sex_inference
PAR_END=${PAR_END:-50000000}   # end of PAR control
BUFFER=${BUFFER:-1000000}      # sites near the PAR boundary are skipped
MINQ=${MINQ:-30}
mkdir -p ${OUTDIR} logs

echo "[$(date)] Summarizing per-sample depth from ${VCF}"

# Mean DP per sample in PAR and Z-specific regions (biallelic SNPs, missing GT/DP skipped)
bcftools view -m2 -M2 -v snps -i "QUAL>=${MINQ}" ${VCF} \
| bcftools query -f '%POS[\t%SAMPLE\t%GT\t%DP]\n' \
| awk -v par=${PAR_END} -v buf=${BUFFER} 'BEGIN{OFS="\t"}
    {
      if ($1 <= par) r = "PAR"
      else if ($1 > par + buf) r = "Z"
      else next
      for (i = 2; i <= NF; i += 3) {
        if ($(i+1) ~ /\./ || $(i+2) == ".") continue
        s[$i] = 1; n[$i, r]++; d[$i, r] += $(i+2)
      }
    }
    END{
      print "sample", "n_PAR", "dp_PAR", "n_Z", "dp_Z", "dp_ratio"
      for (x in s) {
        p = n[x, "PAR"] ? d[x, "PAR"] / n[x, "PAR"] : "NA"
        z = n[x, "Z"]   ? d[x, "Z"]   / n[x, "Z"]   : "NA"
        print x, n[x, "PAR"] + 0, p, n[x, "Z"] + 0, z, (p != "NA" && z != "NA" && p > 0) ? z / p : "NA"
      }
    }' \
| { IFS= read -r hdr; echo "${hdr}"; sort -k1,1; } > ${OUTDIR}/Z_PAR_depth.tsv

echo "[$(date)] Clustering"

Rscript - ${OUTDIR} <<'EOF'
out <- commandArgs(TRUE)[1]
d <- read.delim(file.path(out, "Z_PAR_depth.tsv"))
ok <- !is.na(d$dp_ratio)

set.seed(1)
km <- kmeans(d$dp_ratio[ok], centers = 2, nstart = 50)
fem <- which.min(km$centers)                      # lower-depth cluster = ZW
d$sex <- NA
d$sex[ok] <- ifelse(km$cluster == fem, "F", "M")

write.table(d, file.path(out, "sex_calls.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("Cluster centers (Z/PAR depth):", round(sort(km$centers), 3), "\n")
print(table(d$sex, useNA = "ifany"))

pdf(file.path(out, "sex_calls.pdf"), width = 6, height = 4.5)
cols <- c(F = "#D55E00", M = "#0072B2")
o <- order(d$dp_ratio)
plot(seq_along(o), d$dp_ratio[o], pch = 21, bg = cols[d$sex[o]],
     xlab = "Samples (ranked)", ylab = "Z-specific / PAR mean depth",
     ylim = c(0, max(1.2, d$dp_ratio, na.rm = TRUE)))
abline(h = c(0.5, 1), lty = 2, col = "grey50")
legend("topleft", c("F (ZW)", "M (ZZ)"), pt.bg = cols, pch = 21, bty = "n")
dev.off()
EOF

echo "[$(date)] Done: ${OUTDIR}/sex_calls.tsv and sex_calls.pdf"