"""
DevVarP Pipeline - Script 08: Disease System Analysis
======================================================
Splits G1 Congenital variants by organ system and tests
whether fetal regulatory depletion is consistent across
disease systems or driven by one specific system.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

BASE_DIR = Path.home() / "devvarp_project"
PROC_DIR = BASE_DIR / "data" / "processed"
FIG_DIR  = BASE_DIR / "figures"

print("=" * 60)
print("DevVarP Pipeline - Script 08: Disease System Analysis")
print("=" * 60)

annotated = pd.read_csv(PROC_DIR / "fetal_annotated_variants.tsv", sep='\t')
clinvar   = pd.read_csv(PROC_DIR / "curated_variants.tsv", sep='\t')
clinvar['variant_idx'] = clinvar.index

merged = annotated.merge(
    clinvar[['variant_idx','GeneSymbol','PhenotypeIDS']],
    on='variant_idx', how='left'
)

print(f"Merged variants: {len(merged):,}")
print(f"GeneSymbol non-null: {merged['GeneSymbol'].notna().sum():,}")

g1 = merged[merged['onset_group'] == 'G1_Congenital'].copy()
print(f"G1 Congenital: {len(g1):,}")

HEART_GENES = {
    'GATA4','NKX2-5','TBX5','TBX20','HAND1','HAND2',
    'MYH6','MYH7','ACTC1','TNNT2','TPM1','MYL2','MYL3',
    'SCN5A','KCNQ1','KCNH2','CACNA1C','RYR2','PKP2',
    'DSP','DSG2','DSC2','JUP','TMEM43','LMNA','EMD',
    'FBN1','FBN2','TGFBR1','TGFBR2','SMAD3','ACTA2',
    'MYH11','COL3A1','ELN','NOTCH1','JAG1','PTPN11',
    'RAF1','BRAF','HRAS','KRAS','SOS1','SHOC2','CBL'
}

BRAIN_GENES = {
    'MECP2','CDKL5','ARX','SCN1A','SCN2A','SCN8A',
    'KCNQ2','KCNQ3','STXBP1','SYNGAP1','SHANK3',
    'PTEN','TSC1','TSC2','NF1','NF2','VHL',
    'ASPM','MCPH1','CDK5RAP2','CENPJ','STIL',
    'LIS1','DCX','TUBA1A','TUBB2B','TUBB3',
    'PAX6','SOX2','OTX2','FOXG1','MEF2C',
    'DYRK1A','ANKRD11','KAT6A','SETD5','MED13L',
    'ADGRV1','CLDN14','GJB2','GJB6','MYO7A'
}

METABOLIC_GENES = {
    'CFTR','PAH','PCSK9','LDLR','APOB','APOE',
    'GBA','HEXA','HEXB','GALC','ARSA','NPC1','NPC2',
    'ATP7A','ATP7B','SLC25A13','SLC25A4','POLG',
    'SURF1','SCO2','COX10','BCS1L','LRPPRC',
    'MMUT','PCCA','PCCB','IVD','ACADM','ACADVL',
    'GAA','GLA','GBA2','NAGLU','IDUA','IDS',
    'ABCD1','PEX1','PEX6','PEX10','PEX12',
    'HBB','HBA1','HBA2','G6PD','PKLR'
}

SKELETAL_GENES = {
    'COL1A1','COL1A2','COL2A1','COL9A1','COL9A2',
    'COL10A1','COL11A1','COL11A2','COMP','MATN3',
    'FGFR1','FGFR2','FGFR3','FGFRL1',
    'EXT1','EXT2','GPC3','GPC5','GPC6',
    'RUNX2','SP7','SOX9','SOX5','SOX6',
    'PTCH1','SMO','GLI1','GLI2','GLI3',
    'IHH','SHH','BMP4','BMP5','BMP7','BMPR1A',
    'TRPV4','DTDST','DYMECLIN','SBDS','RMRP'
}

def classify_system(gene):
    if pd.isna(gene): return 'Other Congenital'
    if gene in HEART_GENES:     return 'Congenital Heart'
    if gene in BRAIN_GENES:     return 'Congenital Brain'
    if gene in METABOLIC_GENES: return 'Congenital Metabolic'
    if gene in SKELETAL_GENES:  return 'Congenital Skeletal'
    return 'Other Congenital'

g1['disease_system'] = g1['GeneSymbol'].apply(classify_system)

print("\nG1 variants by disease system:")
print(g1['disease_system'].value_counts())

print("\nFETAL REGULATORY DEPLETION BY DISEASE SYSTEM:")
print(f"  {'System':<25} {'Path N':>7} {'Path%':>7} {'Ben N':>7} {'Ben%':>7} {'OR':>7} {'p':>8}")
print(f"  {'-'*72}")

systems = ['Congenital Heart','Congenital Brain',
           'Congenital Metabolic','Congenital Skeletal','Other Congenital']
results = []

for system in systems:
    subset = g1[g1['disease_system'] == system]
    path   = subset[subset['label'] == 1]
    benign = subset[subset['label'] == 0]

    if len(path) < 5:
        print(f"  {system:<25} insufficient n={len(path)}")
        continue

    path_in    = int(path['in_fetal'].sum())
    path_out   = int(len(path) - path_in)
    benign_in  = int(benign['in_fetal'].sum())
    benign_out = int(len(benign) - benign_in)

    path_pct   = path_in / len(path) * 100
    benign_pct = benign_in / len(benign) * 100 if len(benign) > 0 else 0

    table = [[path_in, path_out],[benign_in, benign_out]]
    or_val, pval = fisher_exact(table)

    print(f"  {system:<25} {len(path):>7} {path_pct:>6.1f}% "
          f"{len(benign):>7} {benign_pct:>6.1f}% "
          f"{or_val:>6.3f}  {pval:>7.4f}")

    results.append({
        'disease_system':   system,
        'n_pathogenic':     len(path),
        'path_fetal_pct':   round(path_pct, 2),
        'n_benign':         len(benign),
        'benign_fetal_pct': round(benign_pct, 2),
        'odds_ratio':       round(or_val, 3),
        'pvalue':           round(pval, 4)
    })

results_df = pd.DataFrame(results)
results_df.to_csv(PROC_DIR / "disease_system_results.tsv", sep='\t', index=False)

if len(results_df) > 0:
    fig, ax = plt.subplots(figsize=(11, 6))
    x     = np.arange(len(results_df))
    width = 0.38

    bars1 = ax.bar(x - width/2, results_df['path_fetal_pct'], width,
                   label='Pathogenic', color='#d73027',
                   edgecolor='black', linewidth=0.9)
    bars2 = ax.bar(x + width/2, results_df['benign_fetal_pct'], width,
                   label='Benign', color='#4575b4',
                   edgecolor='black', linewidth=0.9)

    for bar, n in zip(bars1, results_df['n_pathogenic']):
        ax.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() + 0.3,
                f'n={int(n)}', ha='center', va='bottom',
                fontsize=8, color='#d73027', fontweight='bold')

    for bar, pval in zip(bars1, results_df['pvalue']):
        sig = '***' if pval<0.001 else '**' if pval<0.01 else '*' if pval<0.05 else 'ns'
        ax.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() + 1.8,
                sig, ha='center', va='bottom',
                fontsize=11, fontweight='bold')

    ax.set_xlabel('Congenital Disease System', fontsize=11)
    ax.set_ylabel('% Variants in Fetal Elements', fontsize=11)
    ax.set_title('Fetal Regulatory Depletion Across Congenital Disease Systems\n(G1 Congenital Onset Only)',
                 fontsize=12, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(results_df['disease_system'],
                       rotation=15, ha='right', fontsize=9)
    ax.legend(fontsize=10)
    ax.set_ylim(0, 35)
    ax.yaxis.grid(True, alpha=0.3)

    plt.tight_layout()
    out = FIG_DIR / "figure4_disease_systems.png"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\nFigure saved: {out}")

print("\nScript 08 complete.")
