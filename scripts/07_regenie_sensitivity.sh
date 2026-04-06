#!/usr/bin/env bash
###############################################################################
# 07_regenie_sensitivity.sh
#
# Run regenie v2.2.4 (two-step whole-genome regression) for the European-only
# sensitivity analysis.
#
# IMPORTANT: This script runs on Helix Research Network (HRN) infrastructure.
# Paths marked [HELIX INTERNAL] require access to that environment and cannot
# be executed externally.  See Methods for data access details.
#
# Prerequisites:
#   - Run 07_regenie_sensitivity.R first to generate input files in output/regenie/
#   - regenie v2.2.4 available in PATH (or set REGENIE_BIN below)
#
# Outputs (written to output/regenie/):
#   regenie_step1*                — step-1 null model files (LOCO predictions)
#   regenie_step2_glp1r_*.regenie — step-2 results for GLP1R
#   regenie_step2_gipr_*.regenie  — step-2 results for GIPR
#
# Paper: "Rare GLP1R variants selectively attenuate weight loss with
#         semaglutide versus tirzepatide"
###############################################################################

set -euo pipefail

REGENIE_BIN="${REGENIE_BIN:-regenie}"   # override with full path if needed
REGENIE_VERSION="2.2.4"

OUTDIR="output/regenie"
PHENO="${OUTDIR}/regenie_pheno.txt"
COVAR="${OUTDIR}/regenie_covariates.txt"
ANNO="${OUTDIR}/regenie_anno.txt"
SETLIST="${OUTDIR}/regenie_setlist.txt"
MASKS="${OUTDIR}/regenie_masks.txt"

# Verify version
echo "regenie version check:"
"${REGENIE_BIN}" --version | head -1
echo "Expected: ${REGENIE_VERSION}"
echo ""

###############################################################################
# STEP 1: Whole-genome null model (ridge regression on common variants)
#
# Common-variant PLINK files for step 1 are QC-filtered HRN array genotypes:
#   MAF > 1%, genotype missingness < 5%, HWE p > 1e-15, LD-pruned
# [HELIX INTERNAL: paths to QC-filtered common variant PLINK files]
###############################################################################

echo "=== REGENIE STEP 1 ==="

"${REGENIE_BIN}" \
  --step 1 \
  --bed   "[HELIX INTERNAL: path to QC-filtered common variant PLINK .bed file (without extension)]" \
  --phenoFile "${PHENO}" \
  --covarFile "${COVAR}" \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix "${OUTDIR}/tmp_regenie_lowmem" \
  --out "${OUTDIR}/regenie_step1" \
  --threads 8

echo "Step 1 complete."
echo ""

###############################################################################
# STEP 2: Rare variant burden / SKAT tests for GLP1R
#
# Exome variant files (BGEN format) for step 2 are the HRN exome sequencing
# data filtered to the GLP1R locus (chr6:39,000,000–39,200,000, GRCh38).
# [HELIX INTERNAL: paths to GLP1R and GIPR region BGEN files]
###############################################################################

echo "=== REGENIE STEP 2: GLP1R ==="

"${REGENIE_BIN}" \
  --step 2 \
  --bgen  "[HELIX INTERNAL: path to GLP1R region BGEN file]" \
  --sample "[HELIX INTERNAL: path to BGEN .sample file]" \
  --phenoFile "${PHENO}" \
  --covarFile "${COVAR}" \
  --pred "${OUTDIR}/regenie_step1_pred.list" \
  --anno-file  "${ANNO}" \
  --set-list   "${SETLIST}" \
  --mask-def   "${MASKS}" \
  --aaf-bins 0.001 \
  --build-mask "max" \
  --vc-tests skato,acatv \
  --bsize 200 \
  --minMAC 1 \
  --out "${OUTDIR}/regenie_step2_glp1r" \
  --threads 8

echo "GLP1R step 2 complete."
echo ""

###############################################################################
# STEP 2: Rare variant burden / SKAT tests for GIPR
#
# [HELIX INTERNAL: paths to GIPR region BGEN files]
###############################################################################

echo "=== REGENIE STEP 2: GIPR ==="

"${REGENIE_BIN}" \
  --step 2 \
  --bgen  "[HELIX INTERNAL: path to GIPR region BGEN file]" \
  --sample "[HELIX INTERNAL: path to BGEN .sample file]" \
  --phenoFile "${PHENO}" \
  --covarFile "${COVAR}" \
  --pred "${OUTDIR}/regenie_step1_pred.list" \
  --anno-file  "${ANNO}" \
  --set-list   "${SETLIST}" \
  --mask-def   "${MASKS}" \
  --aaf-bins 0.001 \
  --build-mask "max" \
  --vc-tests skato,acatv \
  --bsize 200 \
  --minMAC 1 \
  --out "${OUTDIR}/regenie_step2_gipr" \
  --threads 8

echo "GIPR step 2 complete."
echo ""
echo "All done. Parse results with scripts/07_regenie_sensitivity.R."
