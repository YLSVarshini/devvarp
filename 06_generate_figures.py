import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
from statsmodels.stats.multitest import multipletests

BASE_DIR = Path.home() / "devvarp_project"
PROC_DIR = BASE_DIR / "data" / "processed"
FIG_DIR  = BASE_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

C_PATH   = '#B2182B'
C_BEN    = '#2166AC'
C_DEPLET = '#B2182B'
C_ENRICH = '#1B7837'
C_NEUT   = '#636363'
C_LINE   = '#252525'
FONT     = 'DejaVu Serif'

plt.rcParams.update({
    'font.family':        FONT,
    'font.size':          11,
    'axes.titlesize':     11,
    'axes.titleweight':   'bold',
    'axes.titlepad':      12,
    'axes.labelsize':     10.5,
    'xtick.labelsize':    9.5,
    'ytick.labelsize':    9.5,
    'legend.fontsize':    9.5,
    'legend.framealpha':  0.9,
    'legend.edgecolor':   '#cccccc',
    'axes.spines.top':    False,
    'axes.spines.right':  False,
    'axes.linewidth':     0.8,
    'figure.dpi':         150,
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.pad_inches': 0.2,
})

def sig_stars(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

enrich  = pd.read_csv(PROC_DIR / "fetal_enrichment.tsv",       sep='\t')
stats   = pd.read_csv(PROC_DIR / "statistical_tests.tsv",      sep='\t')
disease = pd.read_csv(PROC_DIR / "disease_system_results.tsv", sep='\t')

groups  = ['G1_Congenital','G2_Neonatal','G3_Infantile','G4_Childhood']
glabels = ['Congenital\n(G1)','Neonatal\n(G2)','Infantile\n(G3)','Childhood\n(G4)']

ratios      = enrich.set_index('onset_group').loc[groups,'enrichment_ratio'].values
path_pcts   = enrich.set_index('onset_group').loc[groups,'path_fetal_pct'].values
benign_pcts = enrich.set_index('onset_group').loc[groups,'benign_fetal_pct'].values
n_path      = enrich.set_index('onset_group').loc[groups,'n_pathogenic'].values
n_benign    = enrich.set_index('onset_group').loc[groups,'n_benign'].values
pvals_corr  = stats.set_index('onset_group').loc[groups,'pvalue_corrected'].values

_, pvals_dis, _, _ = multipletests(disease['pvalue'].values, method='fdr_bh')
disease['pvalue_corrected'] = pvals_dis
disease['enrichment_ratio'] = disease['path_fetal_pct'] / disease['benign_fetal_pct']

print("Data loaded.")
# FIGURE 1
print("Figure 1...")
fig, ax = plt.subplots(figsize=(7.5, 5.5))
fig.subplots_adjust(top=0.88, bottom=0.16, left=0.13, right=0.96)

bar_colors = [C_DEPLET if r<0.8 else C_ENRICH if r>1.2 else C_NEUT for r in ratios]
bars = ax.bar(glabels, ratios, color=bar_colors, width=0.52,
              edgecolor=C_LINE, linewidth=0.9, zorder=3)

ax.axhline(y=1.0, color=C_LINE, linestyle='--', linewidth=1.4,
           alpha=0.6, zorder=2, xmin=0.02, xmax=0.98)

for bar, n, ratio, p in zip(bars, n_path, ratios, pvals_corr):
    ax.text(bar.get_x() + bar.get_width()/2, ratio/2,
            f'n={int(n)}',
            ha='center', va='center',
            fontsize=9, color='white', fontweight='bold', fontfamily=FONT)
    sig = sig_stars(p)
    ax.text(bar.get_x() + bar.get_width()/2, ratio + 0.05,
            sig, ha='center', va='bottom',
            fontsize=13 if sig!='ns' else 10,
            fontweight='bold' if sig!='ns' else 'normal',
            color=C_LINE, fontfamily=FONT)

ax.set_ylabel('Enrichment Ratio  (Pathogenic % / Benign %)', labelpad=10, fontfamily=FONT)
ax.set_xlabel('Developmental Onset Group', labelpad=8, fontfamily=FONT)
ax.set_title('Fetal Regulatory Element Enrichment in Pathogenic\nNon-Coding Variants by Developmental Onset Group', fontfamily=FONT)
ax.set_ylim(0, 1.95)
ax.set_xlim(-0.5, 3.5)
ax.set_yticks([0, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75])
ax.yaxis.grid(True, alpha=0.2, linewidth=0.6, zorder=0)
ax.set_axisbelow(True)
ax.legend(handles=[
    mpatches.Patch(color=C_DEPLET, label='Depleted  (ratio < 0.8)'),
    mpatches.Patch(color=C_ENRICH, label='Enriched  (ratio > 1.2)'),
    mpatches.Patch(color=C_NEUT,   label='Neutral  (0.8-1.2)'),
], loc='upper left', framealpha=0.9)

fig.text(0.5, 0.01,
         'FDR-corrected  ***p<0.001  **p<0.01  *p<0.05  ns = not significant',
         ha='center', va='bottom', fontsize=8,
         color='#555555', style='italic', fontfamily=FONT)

plt.savefig(FIG_DIR / "figure1_enrichment_ratio.png")
plt.close()
print("  Figure 1 done.")

# FIGURE 2
print("Figure 2...")
fig, ax = plt.subplots(figsize=(6.5, 5.5))
fig.subplots_adjust(top=0.88, bottom=0.12, left=0.20, right=0.90)

heatmap_data = pd.DataFrame({
    'Pathogenic': path_pcts,
    'Benign':     benign_pcts
}, index=glabels)

sns.heatmap(heatmap_data, annot=False, cmap='RdYlGn',
            linewidths=1.5, linecolor='white',
            cbar_kws={'label':'% Variants in Fetal Elements',
                      'shrink':0.82, 'aspect':20},
            ax=ax, vmin=10, vmax=32)

for i in range(len(glabels)):
    for j, (pct, n) in enumerate([
        (path_pcts[i],   int(n_path[i])),
        (benign_pcts[i], int(n_benign[i]))
    ]):
        tc = 'white' if pct > 26 else '#252525'
        ax.text(j+0.5, i+0.38, f'{pct:.1f}%',
                ha='center', va='center',
                fontsize=13, fontweight='bold', fontfamily=FONT, color=tc)
        ax.text(j+0.5, i+0.65, f'(n={n:,})',
                ha='center', va='center',
                fontsize=9, fontfamily=FONT, color=tc)

ax.set_title('Percentage of Variants Overlapping\nFetal Chromatin Accessibility Sites', fontfamily=FONT)
ax.set_xlabel('Variant Class', labelpad=8, fontfamily=FONT)
ax.set_ylabel('Developmental Onset Group', labelpad=8, fontfamily=FONT)
ax.set_xticklabels(['Pathogenic','Benign'], rotation=0, fontsize=11, fontfamily=FONT)
ax.set_yticklabels(glabels, rotation=0, fontsize=10, fontfamily=FONT)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=9)
cbar.set_label('% Variants in Fetal Elements', fontsize=9.5, fontfamily=FONT)

plt.savefig(FIG_DIR / "figure2_heatmap.png")
plt.close()
print("  Figure 2 done.")
# FIGURE 3
print("Figure 3...")
fig = plt.figure(figsize=(13, 6))
gs  = gridspec.GridSpec(1, 2, wspace=0.40,
                         left=0.08, right=0.97,
                         top=0.90, bottom=0.12)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

x     = np.arange(len(glabels))
width = 0.34

b1 = ax1.bar(x - width/2, path_pcts, width,
             color=C_PATH, edgecolor=C_LINE, linewidth=0.9,
             label='Pathogenic', zorder=3)
b2 = ax1.bar(x + width/2, benign_pcts, width,
             color=C_BEN, alpha=0.85, edgecolor=C_LINE, linewidth=0.9,
             label='Benign', zorder=3)

for bar, n in zip(b1, n_path):
    ax1.text(bar.get_x() + bar.get_width()/2,
             bar.get_height() + 0.6,
             f'n={int(n)}',
             ha='center', va='bottom',
             fontsize=8, color=C_PATH, fontweight='bold', fontfamily=FONT)

for bar, n in zip(b2, n_benign):
    ax1.text(bar.get_x() + bar.get_width()/2,
             bar.get_height() + 0.6,
             f'n={int(n)}',
             ha='center', va='bottom',
             fontsize=7.5, color='#1A4A80', fontfamily=FONT)

ax1.set_xlabel('Developmental Onset Group', labelpad=8, fontfamily=FONT)
ax1.set_ylabel('% Variants in Fetal Elements', labelpad=8, fontfamily=FONT)
ax1.set_title('Pathogenic vs Benign Overlap\nwith Fetal Regulatory Elements', fontfamily=FONT)
ax1.set_xticks(x)
ax1.set_xticklabels(glabels, fontsize=9.5, fontfamily=FONT)
ax1.legend(fontsize=10, framealpha=0.9, loc='upper left')
ax1.set_ylim(0, 42)
ax1.set_xlim(-0.5, 3.5)
ax1.yaxis.grid(True, alpha=0.2, linewidth=0.6, zorder=0)
ax1.set_axisbelow(True)

ax2.plot(glabels, ratios, 'o-', color=C_DEPLET,
         linewidth=2.8, markersize=10,
         markerfacecolor='white', markeredgewidth=2.5, zorder=4)
ax2.axhline(y=1.0, color=C_LINE, linestyle='--',
            linewidth=1.4, alpha=0.6, zorder=2,
            xmin=0.02, xmax=0.98)
ax2.fill_between(glabels, ratios, 1.0,
                 where=[r < 1 for r in ratios],
                 alpha=0.10, color=C_DEPLET, label='Depleted region', zorder=1)
ax2.fill_between(glabels, ratios, 1.0,
                 where=[r >= 1 for r in ratios],
                 alpha=0.10, color=C_ENRICH, label='Enriched region', zorder=1)

for i, (ratio, p) in enumerate(zip(ratios, pvals_corr)):
    sig = sig_stars(p)
    col = C_DEPLET if ratio < 1 else C_ENRICH
    ax2.annotate(sig, xy=(i, ratio), xytext=(0, 14),
                 textcoords='offset points', ha='center',
                 fontsize=13 if sig!='ns' else 10,
                 fontweight='bold' if sig!='ns' else 'normal',
                 color=col, fontfamily=FONT)
    ax2.text(i, ratio - 0.055, f'{ratio:.2f}',
             ha='center', va='top', fontsize=9,
             fontfamily=FONT, color=C_LINE)

ax2.set_xlabel('Developmental Onset Group', labelpad=8, fontfamily=FONT)
ax2.set_ylabel('Enrichment Ratio', labelpad=8, fontfamily=FONT)
ax2.set_title('Depletion-to-Enrichment Gradient\nAcross Developmental Onset', fontfamily=FONT)
ax2.legend(fontsize=10, loc='upper left', framealpha=0.9)
ax2.set_ylim(0.55, 1.78)
ax2.set_xticks(range(len(glabels)))
ax2.set_xticklabels(glabels, fontsize=9.5, fontfamily=FONT)
ax2.yaxis.grid(True, alpha=0.2, linewidth=0.6, zorder=0)
ax2.set_axisbelow(True)

plt.savefig(FIG_DIR / "figure3_gradient.png")
plt.close()
print("  Figure 3 done.")

# FIGURE 4
print("Figure 4...")
fig = plt.figure(figsize=(13, 6.5))
gs  = gridspec.GridSpec(1, 2, wspace=0.40,
                         left=0.08, right=0.97,
                         top=0.90, bottom=0.18)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

xd = np.arange(len(disease))
dlabels_short = ['Congenital\nHeart','Congenital\nBrain',
                 'Congenital\nSkeletal','Other\nCongenital']

b1 = ax1.bar(xd - width/2, disease['path_fetal_pct'], width,
             color=C_PATH, edgecolor=C_LINE, linewidth=0.9,
             label='Pathogenic', zorder=3)
b2 = ax1.bar(xd + width/2, disease['benign_fetal_pct'], width,
             color=C_BEN, alpha=0.85, edgecolor=C_LINE, linewidth=0.9,
             label='Benign', zorder=3)

for bar, n in zip(b1, disease['n_pathogenic']):
    ax1.text(bar.get_x() + bar.get_width()/2,
             bar.get_height() + 0.5,
             f'n={int(n)}',
             ha='center', va='bottom',
             fontsize=8, color=C_PATH, fontweight='bold', fontfamily=FONT)

for bar, p in zip(b1, disease['pvalue_corrected']):
    sig = sig_stars(p)
    ax1.text(bar.get_x() + bar.get_width()/2,
             bar.get_height() + 2.2,
             sig, ha='center', va='bottom',
             fontsize=12, fontweight='bold', color=C_LINE, fontfamily=FONT)

ax1.set_xlabel('Congenital Disease System', labelpad=8, fontfamily=FONT)
ax1.set_ylabel('% Variants in Fetal Elements', labelpad=8, fontfamily=FONT)
ax1.set_title('Pathogenic vs Benign Overlap\nby Disease System', fontfamily=FONT)
ax1.set_xticks(xd)
ax1.set_xticklabels(dlabels_short, rotation=0, ha='center',
                    fontsize=9.5, fontfamily=FONT)
ax1.legend(fontsize=10, framealpha=0.9)
ax1.set_ylim(0, 36)
ax1.set_xlim(-0.5, 3.5)
ax1.yaxis.grid(True, alpha=0.2, linewidth=0.6, zorder=0)
ax1.set_axisbelow(True)

ratio_colors = [C_DEPLET if r<0.8 else C_ENRICH if r>1.2 else C_NEUT
                for r in disease['enrichment_ratio']]
rb = ax2.bar(xd, disease['enrichment_ratio'],
             color=ratio_colors, edgecolor=C_LINE,
             linewidth=0.9, width=0.50, zorder=3)

ax2.axhline(y=1.0, color=C_LINE, linestyle='--',
            linewidth=1.4, alpha=0.6, zorder=2,
            xmin=0.02, xmax=0.98)
ax2.axhspan(0.8, 1.2, alpha=0.04, color=C_NEUT, zorder=1)

for bar, ratio, p in zip(rb, disease['enrichment_ratio'],
                          disease['pvalue_corrected']):
    ax2.text(bar.get_x() + bar.get_width()/2,
             ratio + 0.025, f'{ratio:.2f}',
             ha='center', va='bottom',
             fontsize=10, fontweight='bold', color=C_LINE, fontfamily=FONT)
    sig = sig_stars(p)
    ax2.text(bar.get_x() + bar.get_width()/2,
             ratio + 0.10, sig,
             ha='center', va='bottom',
             fontsize=12, fontweight='bold', color=C_LINE, fontfamily=FONT)

ax2.set_xlabel('Congenital Disease System', labelpad=8, fontfamily=FONT)
ax2.set_ylabel('Enrichment Ratio (Pathogenic % / Benign %)', labelpad=8, fontfamily=FONT)
ax2.set_title('Fetal Regulatory Enrichment Ratio\nby Disease System (FDR Corrected)', fontfamily=FONT)
ax2.set_xticks(xd)
ax2.set_xticklabels(dlabels_short, rotation=0, ha='center',
                    fontsize=9.5, fontfamily=FONT)
ax2.set_ylim(0, 1.65)
ax2.set_xlim(-0.5, 3.5)
ax2.yaxis.grid(True, alpha=0.2, linewidth=0.6, zorder=0)
ax2.set_axisbelow(True)
ax2.legend(handles=[
    mpatches.Patch(color=C_DEPLET, label='Depleted  (ratio < 0.8)'),
    mpatches.Patch(color=C_ENRICH, label='Enriched  (ratio > 1.2)'),
    mpatches.Patch(color=C_NEUT,   label='Neutral  (0.8-1.2)'),
], fontsize=9, loc='upper right', framealpha=0.9)

plt.savefig(FIG_DIR / "figure4_disease_systems.png")
plt.close()
print("  Figure 4 done.")
print("\nAll figures saved to:", FIG_DIR)
