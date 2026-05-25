"""
DevVarP Pipeline - Script 03b: Fetal Regulatory Annotation
===========================================================
Re-runs regulatory overlap analysis using human fetal
chromatin accessibility sites (Cao et al. 2020, Science)
instead of generic ENCODE cCREs.
This gives us fetal-specific regulatory context — the key
biological layer missing from existing tool evaluations.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import pyranges as pr

BASE_DIR  = Path.home() / "devvarp_project"
PROC_DIR  = BASE_DIR / "data" / "processed"
RAW_DIR   = BASE_DIR / "data" / "raw"

VARIANTS_FILE = PROC_DIR / "curated_variants.tsv"
FETAL_FILE    = RAW_DIR  / "fetal_enhancers.txt"

print("=" * 60)
print("DevVarP Pipeline - Script 03b: Fetal Annotation")
print("=" * 60)

print("\nSTEP 1: Loading curated variants...")
variants = pd.read_csv(VARIANTS_FILE, sep='\t', low_memory=False)
print(f"  Loaded {len(variants):,} variants")
print("\nSTEP 2: Loading fetal chromatin accessibility sites...")
fetal_raw = pd.read_csv(FETAL_FILE, sep='\t', low_memory=False)
print(f"  Loaded {len(fetal_raw):,} fetal chromatin sites")
print(f"  Columns: {list(fetal_raw.columns)}")
print(f"  Sample peaks: {fetal_raw['peak'].head(3).tolist()}")

# Parse chromosome from peak_id format: chr1_9992_10688
fetal_raw['chrom'] = fetal_raw['peak'].apply(
    lambda x: str(x).split('_')[0])
fetal_raw['fetal_start'] = fetal_raw['start'].astype(int)
fetal_raw['fetal_end']   = fetal_raw['end'].astype(int)

# Keep only standard chromosomes
valid_chroms = [f'chr{i}' for i in list(range(1,23))+['X','Y']]
fetal_raw = fetal_raw[fetal_raw['chrom'].isin(valid_chroms)]
print(f"  After filtering standard chroms: {len(fetal_raw):,} sites")
print(f"  Chromosomes covered: {sorted(fetal_raw['chrom'].unique())[:5]}...")

print("\nSTEP 3: Creating genomic ranges...")
variants_clean = variants.dropna(subset=['chrom','Start','Stop'])
variants_clean = variants_clean[
    variants_clean['chrom'].astype(str).str.startswith('chr')]

var_ranges = pr.PyRanges(pd.DataFrame({
    'Chromosome':  variants_clean['chrom'].astype(str),
    'Start':       variants_clean['Start'].astype(int),
    'End':         variants_clean['Stop'].astype(int),
    'variant_idx': variants_clean.index,
    'onset_group': variants_clean['onset_group'],
    'label':       variants_clean['label']
}))

fetal_ranges = pr.PyRanges(pd.DataFrame({
    'Chromosome': fetal_raw['chrom'].astype(str),
    'Start':      fetal_raw['fetal_start'],
    'End':        fetal_raw['fetal_end']
}))

print(f"  Variant ranges: {len(var_ranges):,}")
print(f"  Fetal ranges:   {len(fetal_ranges):,}")

print("\nSTEP 4: Finding overlaps with fetal elements...")
overlaps = var_ranges.join(fetal_ranges, how='left')
overlaps_df = overlaps.df
print(f"  Total overlap records: {len(overlaps_df):,}")
print("\nSTEP 5: Calculating fetal regulatory enrichment per onset group...")

overlaps_df['in_fetal'] = (
    overlaps_df['Start_b'].notna() &
    (overlaps_df['Start_b'] != -1)
)

dedup = overlaps_df.sort_values('in_fetal', ascending=False)
dedup = dedup.drop_duplicates(subset='variant_idx', keep='first')

print(f"\n  FETAL REGULATORY ENRICHMENT BY ONSET GROUP:")
print(f"  {'Group':<20} {'Path N':>8} {'Path %Fetal':>12} {'Benign N':>10} {'Ben %Fetal':>12}")
print(f"  {'-'*66}")

results = []
for group in ['G1_Congenital','G2_Neonatal','G3_Infantile','G4_Childhood']:
    subset = dedup[dedup['onset_group'] == group]
    path   = subset[subset['label'] == 1]
    benign = subset[subset['label'] == 0]

    path_pct   = path['in_fetal'].sum()   / len(path)   * 100 if len(path)   > 0 else 0
    benign_pct = benign['in_fetal'].sum() / len(benign) * 100 if len(benign) > 0 else 0

    print(f"  {group:<20} {len(path):>8} {path_pct:>11.1f}%  {len(benign):>9} {benign_pct:>11.1f}%")

    results.append({
        'onset_group':         group,
        'n_pathogenic':        len(path),
        'path_fetal_pct':      round(path_pct,  2),
        'n_benign':            len(benign),
        'benign_fetal_pct':    round(benign_pct, 2),
        'enrichment_ratio':    round(path_pct / benign_pct, 3) if benign_pct > 0 else np.nan
    })

results_df = pd.DataFrame(results)

print(f"\n  ENRICHMENT RATIOS (pathogenic % / benign %):")
print(f"  {'Group':<20} {'Ratio':>8} {'Interpretation':>20}")
print(f"  {'-'*52}")
for _, row in results_df.iterrows():
    ratio = row['enrichment_ratio']
    if ratio > 1.2:
        interp = "ENRICHED in fetal"
    elif ratio < 0.8:
        interp = "DEPLETED in fetal"
    else:
        interp = "neutral"
    print(f"  {row['onset_group']:<20} {ratio:>8.3f} {interp:>20}")

out_fetal    = PROC_DIR / "fetal_annotated_variants.tsv"
out_enrichment = PROC_DIR / "fetal_enrichment.tsv"

dedup.to_csv(out_fetal, sep='\t', index=False)
results_df.to_csv(out_enrichment, sep='\t', index=False)

print(f"\n  Saved to: {out_fetal}")
print(f"  Saved to: {out_enrichment}")
print("\nScript 03b complete. Run Script 05 next:")
print("  python 05_analyze_plot.py")
