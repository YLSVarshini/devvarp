"""
Permutation test to validate G1 Congenital depletion finding.
Shuffles onset group labels 1000 times and recalculates
enrichment ratio to confirm result is not due to chance.
"""
import pandas as pd
import numpy as np
from pathlib import Path

BASE_DIR = Path.home() / "devvarp_project"
PROC_DIR = BASE_DIR / "data" / "processed"

df = pd.read_csv(PROC_DIR / "fetal_annotated_variants.tsv", sep='\t')
print(f"Loaded {len(df):,} variants")

# Observed enrichment ratio for G1
g1 = df[df['onset_group'] == 'G1_Congenital']
path   = g1[g1['label'] == 1]
benign = g1[g1['label'] == 0]
observed_ratio = (path['in_fetal'].mean() / benign['in_fetal'].mean())
print(f"Observed G1 enrichment ratio: {observed_ratio:.4f}")

# Permutation test — shuffle labels 1000 times
print("Running 1000 permutations...")
n_permutations = 1000
perm_ratios = []

g1_labels = g1['label'].values
g1_fetal  = g1['in_fetal'].values
n = len(g1_labels)

for i in range(n_permutations):
    shuffled = np.random.permutation(g1_labels)
    perm_path   = g1_fetal[shuffled == 1]
    perm_benign = g1_fetal[shuffled == 0]
    if len(perm_path) > 0 and len(perm_benign) > 0:
        ratio = perm_path.mean() / perm_benign.mean()
        perm_ratios.append(ratio)

perm_ratios = np.array(perm_ratios)

# P-value = proportion of permutations with ratio <= observed
perm_pvalue = (perm_ratios <= observed_ratio).mean()

print(f"\nPERMUTATION TEST RESULTS:")
print(f"  Observed ratio:     {observed_ratio:.4f}")
print(f"  Permutation mean:   {perm_ratios.mean():.4f}")
print(f"  Permutation std:    {perm_ratios.std():.4f}")
print(f"  Permutation p-val:  {perm_pvalue:.4f}")
print(f"  Significant:        {perm_pvalue < 0.05}")

# Save
results = pd.DataFrame({
    'test': ['permutation'],
    'observed_ratio': [observed_ratio],
    'perm_mean': [perm_ratios.mean()],
    'perm_std':  [perm_ratios.std()],
    'perm_pvalue': [perm_pvalue]
})
results.to_csv(PROC_DIR / "permutation_results.tsv", sep='\t', index=False)
print("\nSaved to permutation_results.tsv")
