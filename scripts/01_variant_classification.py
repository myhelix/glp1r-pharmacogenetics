###############################################################################
# 01_variant_classification.py
#
# Extract rare coding variant carriers for GLP1R and GIPR from pre-filtered
# exome sequencing data (Helix Research Network, 2025 Q3 release).
#
# Variants are pre-filtered to:
#   - MAF < 0.1% in gnomAD v4.1 (each ancestry group)
#   - MANE Select transcript coding variants
#   - Unrelated individuals per ancestry group
#
# Outputs:
#   - glp1r_rare_coding_carriers.tsv  (one row per carrier-variant pair)
#   - gipr_rare_coding_carriers.tsv
#
# IMPORTANT: This script runs on Helix Research Network infrastructure using
# Apache Spark / Hail 0.2.x. It cannot be executed outside that environment.
# See Methods for data access details.
#
# Paper: "Rare GLP1R variants selectively attenuate weight loss with
#         semaglutide versus tirzepatide"
###############################################################################

import hail as hl
import pandas as pd

# Initialize Hail
# [HELIX INTERNAL: adjust master, tmp_dir for your cluster configuration]
hl.init(
    master='local[8]',
    default_reference='GRCh38',
    tmp_dir='[HELIX INTERNAL: temporary directory path]'
)

###############################################################################
# PARAMETERS
###############################################################################

GENES = ['GLP1R', 'GIPR']

# [HELIX INTERNAL: pre-filtered matrix tables are stored in internal storage.
#  Each matrix table contains all-cohort rare coding variants for the
#  respective gene, filtered as described above. See Methods.]
MT_PATHS = {
    'GLP1R': '[HELIX INTERNAL: path to GLP1R pre-filtered matrix table]',
    'GIPR':  '[HELIX INTERNAL: path to GIPR pre-filtered matrix table]',
}

# Output files (relative to repository root)
OUTPUT_PATHS = {
    'GLP1R': 'output/glp1r_rare_coding_carriers.tsv',
    'GIPR':  'output/gipr_rare_coding_carriers.tsv',
}

###############################################################################
# EXTRACTION FUNCTION
###############################################################################

def extract_carriers(gene: str) -> None:
    """
    Load pre-filtered matrix table, filter to gene of interest,
    retain only carrier entries (n_alt_alleles > 0), and export.

    Fields exported per carrier-variant pair:
      d_id                    - de-identified individual ID
      variant_id              - chr:pos:ref:alt (GRCh38)
      n_alt_alleles           - 1 (het) or 2 (hom alt)
      MANE_CT_gene            - gene symbol
      MANE_CT_consequence     - VEP consequence (MANE Select transcript)
      MANE_CT_HGVSP_custom    - protein-level HGVS notation
      REVEL                   - REVEL pathogenicity score (missense)
      MANE_CT_am_pathogenicity - AlphaMissense pathogenicity score
      MANE_CT_lof             - predicted loss-of-function flag (LOFTEE)
    """
    print(f"\n{'='*60}")
    print(f"Processing {gene}")
    print('='*60)

    # Load pre-filtered matrix table
    mt = hl.read_matrix_table(MT_PATHS[gene])

    n_variants = mt.count_rows()
    n_samples  = mt.count_cols()
    print(f"  Matrix table: {n_variants:,} variants x {n_samples:,} samples")

    # Filter to gene of interest
    mt = mt.filter_rows(hl.literal([gene]).contains(mt.MANE_CT_gene))
    n_gene = mt.count_rows()
    print(f"  {gene} variants: {n_gene:,}")

    if n_gene == 0:
        print(f"  WARNING: No {gene} variants found. Check MT_PATHS['{gene}'].")
        return

    # Construct variant ID: chr:pos:ref:alt
    mt = mt.annotate_rows(
        variant_id = hl.str(mt.locus) + ':' + mt.alleles[0] + ':' + mt.alleles[1]
    )

    # Count alt alleles per entry (0 = hom ref, 1 = het, 2 = hom alt)
    mt = mt.annotate_entries(
        n_alt_alleles = mt.GT.n_alt_alleles()
    )

    # Retain only carriers (het or hom alt)
    mt = mt.filter_entries(mt.n_alt_alleles > 0)

    n_carriers = mt.count_cols()
    print(f"  Carriers: {n_carriers:,}")

    # Export flat table of carrier-variant pairs
    entries = mt.entries()
    entries = entries.key_by()
    entries = entries.select(
        'd_id',
        'variant_id',
        'n_alt_alleles',
        'MANE_CT_gene',
        'MANE_CT_consequence',
        'MANE_CT_HGVSP_custom',
        'REVEL',
        'MANE_CT_am_pathogenicity',
        'MANE_CT_lof'
    )

    output_path = OUTPUT_PATHS[gene]
    entries.export(output_path)
    print(f"  Exported to: {output_path}")


###############################################################################
# RUN
###############################################################################

for gene in GENES:
    extract_carriers(gene)

print("\nDone.")
