"""
DevVarP Pipeline - Script 02: Variant Curation
Filters ClinVar to high-confidence non-coding variants
in pediatric-onset rare diseases, stratified by
HPO age-of-onset terms into four developmental groups.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import re
from collections import Counter

BASE_DIR    = Path.home() / "devvarp_project"
RAW_DIR     = BASE_DIR / "data" / "raw"
PROC_DIR    = BASE_DIR / "data" / "processed"
PROC_DIR.mkdir(parents=True, exist_ok=True)

CLINVAR_FILE = RAW_DIR / "variant_summary.txt"
HPO_FILE     = RAW_DIR / "genes_to_phenotype.txt"
HPO_OBO      = RAW_DIR / "hp.obo"

print("=" * 60)
print("DevVarP Pipeline - Script 02: Variant Curation")
print("=" * 60)
ONSET_GROUPS = {
    'HP:0003577': 'G1_Congenital',
    'HP:0003623': 'G2_Neonatal',
    'HP:0003593': 'G3_Infantile',
    'HP:0011463': 'G3_Infantile',
    'HP:0003621': 'G4_Childhood',
    'HP:0030674': 'G1_Congenital',
}

print("\nSTEP 1: Parsing HPO terms...")
onset_term_names = {}
current_id = None
with open(HPO_OBO, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('id: HP:'):
            current_id = line.replace('id: ', '')
        elif line.startswith('name: ') and current_id:
            onset_term_names[current_id] = line.replace('name: ', '')
print(f"  Loaded {len(onset_term_names)} HPO terms")

print("\nSTEP 2: Loading HPO gene-phenotype annotations...")
hpo_df = pd.read_csv(HPO_FILE, sep='\t', comment='#',
    names=['ncbi_gene_id','gene_symbol','hpo_id',
           'hpo_name','frequency','disease_id'],
    low_memory=False)
print(f"  Loaded {len(hpo_df):,} gene-phenotype associations")

gene_onset = {}
for _, row in hpo_df.iterrows():
    if row['hpo_id'] in ONSET_GROUPS:
        gene = row['gene_symbol']
        group = ONSET_GROUPS[row['hpo_id']]
        if gene not in gene_onset:
            gene_onset[gene] = group
        else:
            priority = {'G1_Congenital':1,'G2_Neonatal':2,
                       'G3_Infantile':3,'G4_Childhood':4}
            if priority.get(group,99) < priority.get(gene_onset[gene],99):
                gene_onset[gene] = group

print(f"  Genes with pediatric onset: {len(gene_onset):,}")
group_counts = Counter(gene_onset.values())
for group, count in sorted(group_counts.items()):
    print(f"    {group}: {count} genes")
print("\nSTEP 3: Loading and filtering ClinVar...")
print("  (This may take 2-3 minutes)")

COLS_NEEDED = [
    'Name', 'GeneSymbol', 'ClinicalSignificance',
    'ReviewStatus', 'Type', 'Chromosome',
    'Start', 'Stop', 'Assembly', 'PhenotypeIDS'
]

PATHOGENIC_TERMS = {
    'Pathogenic', 'Likely pathogenic',
    'Pathogenic/Likely pathogenic'
}
BENIGN_TERMS = {
    'Benign', 'Likely benign',
    'Benign/Likely benign'
}
GOOD_REVIEW = {
    'criteria provided, multiple submitters, no conflicts',
    'reviewed by expert panel',
    'practice guideline'
}

chunk_size = 100_000
chunks_kept = []
total_rows = 0
kept_rows = 0

for chunk in pd.read_csv(
    CLINVAR_FILE, sep='\t', low_memory=False,
    chunksize=chunk_size,
    usecols=lambda c: c in COLS_NEEDED,
    on_bad_lines='skip'
):
    total_rows += len(chunk)
    chunk = chunk[chunk['Assembly'] == 'GRCh38']
    chunk = chunk[chunk['ReviewStatus'].isin(GOOD_REVIEW)]
    chunk = chunk[chunk['ClinicalSignificance'].isin(
        PATHOGENIC_TERMS | BENIGN_TERMS)]
    chunk = chunk[chunk['GeneSymbol'].notna()]
    chunk = chunk[chunk['GeneSymbol'] != '-']
    chunk = chunk[chunk['GeneSymbol'].isin(gene_onset.keys())]
    if len(chunk) > 0:
        chunks_kept.append(chunk)
        kept_rows += len(chunk)
    if total_rows % 500_000 == 0:
        print(f"  Processed {total_rows:,} rows, kept {kept_rows:,}...")

print(f"  Total rows processed: {total_rows:,}")
print(f"  Rows after filters: {kept_rows:,}")
clinvar = pd.concat(chunks_kept, ignore_index=True)
print("\nSTEP 4: Identifying non-coding variants...")

def is_noncoding(name):
    if pd.isna(name):
        return False
    name = str(name)
    intronic = re.search(r'c\.\d+[+\-](\d+)', name)
    if intronic:
        offset = int(intronic.group(1))
        if offset > 2:
            return True
    if 'c.*' in name:
        return True
    if re.search(r'c\.-\d+', name):
        return True
    if name.startswith('n.'):
        return True
    if 'upstream' in name.lower():
        return True
    if 'regulatory' in name.lower():
        return True
    return False

clinvar['is_noncoding'] = clinvar['Name'].apply(is_noncoding)
clinvar_nc = clinvar[clinvar['is_noncoding']].copy()
print(f"  All filtered variants: {len(clinvar):,}")
print(f"  Non-coding variants: {len(clinvar_nc):,}")

print("\nSTEP 5: Assigning onset groups and labels...")
clinvar_nc['onset_group'] = clinvar_nc['GeneSymbol'].map(gene_onset)
clinvar_nc['label'] = clinvar_nc['ClinicalSignificance'].apply(
    lambda x: 1 if x in PATHOGENIC_TERMS else 0)
clinvar_nc['chrom'] = clinvar_nc['Chromosome'].astype(str)
clinvar_nc['chrom'] = clinvar_nc['chrom'].apply(
    lambda x: f"chr{x}" if not str(x).startswith('chr') else x)
clinvar_nc = clinvar_nc.dropna(subset=['Start','Stop'])
clinvar_nc['Start'] = clinvar_nc['Start'].astype(int)
clinvar_nc['Stop']  = clinvar_nc['Stop'].astype(int)

print("\nSTEP 6: Summary and saving...")
print(f"\n  FINAL DATASET SUMMARY:")
print(f"  Total non-coding variants: {len(clinvar_nc):,}")
print(f"  Pathogenic: {clinvar_nc['label'].sum():,}")
print(f"  Benign: {(clinvar_nc['label']==0).sum():,}")
print()
print("  By onset group:")
for group in ['G1_Congenital','G2_Neonatal','G3_Infantile','G4_Childhood']:
    subset = clinvar_nc[clinvar_nc['onset_group']==group]
    path = subset['label'].sum()
    benign = (subset['label']==0).sum()
    print(f"    {group}: {len(subset):,} total ({path} pathogenic, {benign} benign)")

out_file = PROC_DIR / "curated_variants.tsv"
clinvar_nc.to_csv(out_file, sep='\t', index=False)
print(f"\n  Saved to: {out_file}")
print("\nScript 02 complete. Run Script 03 next:")
print("  python 03_annotate_regulatory.py")
