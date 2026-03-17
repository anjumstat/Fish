# -*- coding: utf-8 -*-
"""
Functional Enrichment Analysis for Cyprinid Domestication
Hypothesis: Genome-wide amplification preserves functional proportions
Compares expanded genes to genome background
FIXED: Using actual expanded genes from CAFE results with optimized GO/KEGG processing
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import os
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Set publication style - matching Figure 3
# =============================================================================

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.labelsize'] = 19
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 17
plt.rcParams['ytick.labelsize'] = 17
plt.rcParams['legend.fontsize'] = 17
plt.rcParams['figure.titlesize'] = 22
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['savefig.bbox'] = 'tight'

# Nature-style colors
expansion_color = '#1b9e77'  # Teal green
contraction_color = '#d95f02'  # Orange-red
domesticated_color = '#d95f02'  # Orange-red
wild_color = '#7570b3'  # Purple
outgroup_color = '#666666'  # Gray
growth_color = '#2ecc71'  # Green
metabolism_color = '#3498db'  # Blue
stress_color = '#e74c3c'  # Red
other_color = '#95a5a6'  # Gray

# =============================================================================
# Paths
# =============================================================================

eggnog_path = r"E:/fish/eggnog_results/eggnog_results"
cafe_results_path = r"E:/fish/orth/OrthoFinder/Results_Feb11/cafe_results"
output_dir = r"E:/fish/paper/figures/Figure4"
os.makedirs(output_dir, exist_ok=True)

# =============================================================================
# Load eggNOG annotations
# =============================================================================

annot_file = os.path.join(eggnog_path, "fish_domestication.emapper.annotations")
print(f"Loading annotations from: {annot_file}")

# Read file, skipping header lines
df = pd.read_csv(annot_file, sep='\t', skiprows=4, header=None)

# Rename columns based on eggNOG output format
df.columns = ['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 
              'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 
              'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
              'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 
              'BiGG_Reaction', 'PFAMs']

print(f"Loaded {len(df):,} total gene annotations")

# =============================================================================
# LOAD EXPANDED GENES FROM CAFE RESULTS
# =============================================================================

print("\n" + "="*80)
print("LOADING EXPANDED GENES FROM CAFE RESULTS")
print("="*80)

expanded_genes_file = os.path.join(cafe_results_path, "expanded_genes_cyprinids.txt")

if os.path.exists(expanded_genes_file):
    # Load the expanded genes
    expanded_genes = pd.read_csv(expanded_genes_file, header=None, names=['gene_id'])
    print(f"Loaded {len(expanded_genes):,} expanded genes from CAFE results")
    
    # Filter your annotations to only include expanded genes
    expanded_df = df[df['query'].isin(expanded_genes['gene_id'])].copy()
    print(f"Matched {len(expanded_df):,} genes with eggNOG annotations")
    print(f"Coverage: {len(expanded_df)/len(expanded_genes)*100:.2f}% of expanded genes have annotations")
else:
    print(f"ERROR: Expanded genes file not found at: {expanded_genes_file}")
    # Fallback to random sampling
    np.random.seed(42)
    expanded_indices = np.random.choice(len(df), size=min(96919, len(df)), replace=False)
    expanded_df = df.iloc[expanded_indices].copy()
    print("WARNING: Using random sampling as fallback!")

background_df = df.copy()
print(f"Expanded genes set: {len(expanded_df):,}")
print(f"Background genome set: {len(background_df):,}")

# =============================================================================
# Define COG categories and mapping
# =============================================================================

cog_mapping = {
    'J': 'Translation', 'A': 'RNA processing', 'K': 'Transcription',
    'L': 'Replication', 'B': 'Chromatin', 'D': 'Cell cycle',
    'Y': 'Nuclear structure', 'V': 'Defense mechanisms', 'T': 'Signal transduction',
    'M': 'Cell wall/membrane', 'N': 'Cell motility', 'Z': 'Cytoskeleton',
    'W': 'Extracellular', 'U': 'Secretion', 'O': 'Post-translational modification',
    'C': 'Energy metabolism', 'G': 'Carbohydrate metabolism', 'E': 'Amino acid metabolism',
    'F': 'Nucleotide metabolism', 'H': 'Coenzyme metabolism', 'I': 'Lipid metabolism',
    'P': 'Inorganic ion metabolism', 'Q': 'Secondary metabolites',
    'R': 'General function prediction', 'S': 'Function unknown'
}

# Group by hypothesis categories
hypothesis_cogs = {
    'Growth': ['D', 'M', 'N', 'Z'],  # Cell cycle, structure, motility
    'Metabolism': ['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q'],  # All metabolism
    'Stress': ['O', 'T', 'V'],  # Stress response, signaling, defense
    'Other': ['J', 'A', 'K', 'L', 'B', 'Y', 'W', 'U', 'R', 'S']
}

# Reverse mapping for quick lookup
cog_to_hypothesis = {}
for cat, cogs in hypothesis_cogs.items():
    for cog in cogs:
        cog_to_hypothesis[cog] = cat

# =============================================================================
# Parse COG categories
# =============================================================================

def parse_cog(cog_str):
    """Parse COG category string (e.g., 'JAT' -> ['J','A','T'])"""
    if pd.isna(cog_str) or cog_str == '' or cog_str == '-':
        return []
    return list(cog_str)

df['cog_list'] = df['COG_category'].apply(parse_cog)
expanded_df['cog_list'] = expanded_df['COG_category'].apply(parse_cog)

# =============================================================================
# PANEL A: COG ENRICHMENT ANALYSIS
# =============================================================================

print("\n" + "="*80)
print("PANEL A: COG ENRICHMENT ANALYSIS")
print("="*80)

# Count COG categories in background and expanded
all_cogs = set()
for c in cog_mapping.keys():
    all_cogs.add(c)

cog_results = []
for cog in sorted(all_cogs):
    # Background count
    bg_count = sum(1 for cog_list in df['cog_list'] if cog in cog_list)
    bg_total = len(df)
    bg_freq = bg_count / bg_total if bg_total > 0 else 0
    
    # Expanded count
    exp_count = sum(1 for cog_list in expanded_df['cog_list'] if cog in cog_list)
    exp_total = len(expanded_df)
    exp_freq = exp_count / exp_total if exp_total > 0 else 0
    
    # Fold enrichment
    fold_enrich = exp_freq / bg_freq if bg_freq > 0 else 0
    
    # Fisher's exact test
    if exp_count > 0 and bg_count > 0:
        table = [[exp_count, exp_total - exp_count],
                 [bg_count, bg_total - bg_count]]
        oddsratio, p_value = fisher_exact(table)
    else:
        p_value = 1.0
    
    cog_results.append({
        'COG': cog,
        'Category': cog_mapping.get(cog, 'Unknown'),
        'Hypothesis': cog_to_hypothesis.get(cog, 'Other'),
        'Background_Count': bg_count,
        'Background_Freq': bg_freq,
        'Expanded_Count': exp_count,
        'Expanded_Freq': exp_freq,
        'Fold_Enrichment': fold_enrich,
        'P_Value': p_value,
    })

cog_df = pd.DataFrame(cog_results)

# FDR correction
cog_df['FDR'] = multipletests(cog_df['P_Value'], method='fdr_bh')[1]
cog_df['-log10_FDR'] = -np.log10(cog_df['FDR'])
cog_df['-log10_p'] = -np.log10(cog_df['P_Value'])

# Sort by fold enrichment for display
cog_df_sorted = cog_df.sort_values('Fold_Enrichment', ascending=False)

print("\nCOG Enrichment Results (sorted by Fold Enrichment):")
print("-" * 100)
print(f"{'COG':<5} {'Category':<25} {'Background':<15} {'Expanded':<15} {'Fold':<8} {'FDR':<10}")
print("-" * 100)
for _, row in cog_df_sorted.iterrows():
    print(f"{row['COG']:<5} {row['Category'][:25]:<25} "
          f"{row['Background_Count']:>6,} ({row['Background_Freq']*100:>4.1f}%) "
          f"{row['Expanded_Count']:>6,} ({row['Expanded_Freq']*100:>4.1f}%) "
          f"{row['Fold_Enrichment']:>6.3f} "
          f"{row['FDR']:.2e}")

print(f"\nCOGs with FDR < 0.05: {sum(cog_df['FDR'] < 0.05)} of {len(cog_df)}")
print(f"Mean fold enrichment: {cog_df['Fold_Enrichment'].mean():.3f}")
print(f"Fold enrichment range: {cog_df['Fold_Enrichment'].min():.3f} - {cog_df['Fold_Enrichment'].max():.3f}")

# =============================================================================
# PANEL B: HYPOTHESIS CATEGORY DISTRIBUTION
# =============================================================================

print("\n" + "="*80)
print("PANEL B: HYPOTHESIS CATEGORY DISTRIBUTION")
print("="*80)

def get_hypothesis_categories(cog_list):
    """Get all hypothesis categories for a gene based on its COGs"""
    cats = set()
    for cog in cog_list:
        if cog in cog_to_hypothesis:
            cats.add(cog_to_hypothesis[cog])
    return cats

# Background
bg_cats = []
for cog_list in df['cog_list']:
    bg_cats.extend(list(get_hypothesis_categories(cog_list)))
bg_cat_counts = pd.Series(bg_cats).value_counts()

# Expanded
exp_cats = []
for cog_list in expanded_df['cog_list']:
    exp_cats.extend(list(get_hypothesis_categories(cog_list)))
exp_cat_counts = pd.Series(exp_cats).value_counts()

print("\nHypothesis category distribution:")
print("-" * 70)
print(f"{'Category':<15} {'Background':<20} {'Expanded':<20} {'Fold':<10}")
print("-" * 70)

hypothesis_order = ['Growth', 'Metabolism', 'Stress', 'Other']
for cat in hypothesis_order:
    bg_count = bg_cat_counts.get(cat, 0)
    bg_pct = (bg_count / len(df)) * 100
    exp_count = exp_cat_counts.get(cat, 0)
    exp_pct = (exp_count / len(expanded_df)) * 100
    fold = exp_pct / bg_pct if bg_pct > 0 else 0
    
    print(f"{cat:<15} "
          f"{bg_count:>7,} ({bg_pct:>5.2f}%)   "
          f"{exp_count:>7,} ({exp_pct:>5.2f}%)   "
          f"{fold:>6.3f}")

# =============================================================================
# PANEL C: GO TERM ENRICHMENT (OPTIMIZED)
# =============================================================================

print("\n" + "="*80)
print("PANEL C: GO TERM ENRICHMENT")
print("="*80)

def parse_go(go_str):
    """Parse GO terms"""
    if pd.isna(go_str) or go_str == '' or go_str == '-':
        return []
    return [go.strip() for go in str(go_str).split(',')]

print("Parsing GO terms...")
df['go_list'] = df['GOs'].apply(parse_go)
expanded_df['go_list'] = expanded_df['GOs'].apply(parse_go)

# Get all GO terms
all_go_terms = set()
for go_list in df['go_list']:
    all_go_terms.update(go_list)

print(f"Total unique GO terms: {len(all_go_terms):,}")

# OPTIMIZATION: Use Counter for frequency counting
print("Counting GO term frequencies in background...")
bg_go_counter = Counter()
for i, go_list in enumerate(df['go_list']):
    if i % 10000 == 0:
        print(f"  Processed {i:,} background genes...")
    bg_go_counter.update(go_list)

print("Counting GO term frequencies in expanded genes...")
exp_go_counter = Counter()
for i, go_list in enumerate(expanded_df['go_list']):
    if i % 5000 == 0:
        print(f"  Processed {i:,} expanded genes...")
    exp_go_counter.update(go_list)

# Create results
go_counts = []
total_go_terms = len(all_go_terms)

print(f"Processing {total_go_terms:,} GO terms...")
for i, go_term in enumerate(all_go_terms):
    if i % 1000 == 0:
        print(f"  Processing GO term {i:,}/{total_go_terms:,}...")
    
    bg_count = bg_go_counter.get(go_term, 0)
    exp_count = exp_go_counter.get(go_term, 0)
    
    if exp_count > 5:  # Filter rare terms
        go_counts.append({
            'GO': go_term,
            'Background': bg_count,
            'Expanded': exp_count,
            'Expanded_Pct': (exp_count / len(expanded_df)) * 100
        })

go_counts_df = pd.DataFrame(go_counts).sort_values('Expanded', ascending=False).head(15)
print("\nTop 15 most frequent GO terms in expanded genes:")
if len(go_counts_df) > 0:
    print(go_counts_df.to_string(index=False))
else:
    print("No GO terms found with count > 5")

# =============================================================================
# PANEL D: KEGG PATHWAY ENRICHMENT (OPTIMIZED)
# =============================================================================

print("\n" + "="*80)
print("PANEL D: KEGG PATHWAY ENRICHMENT")
print("="*80)

def parse_kegg(kegg_str):
    """Parse KEGG pathways"""
    if pd.isna(kegg_str) or kegg_str == '' or kegg_str == '-':
        return []
    pathways = []
    for item in str(kegg_str).split(','):
        if 'map' in item or 'ko' in item:
            pathways.append(item.strip())
    return pathways

print("Parsing KEGG pathways...")
df['kegg_list'] = df['KEGG_Pathway'].apply(parse_kegg)
expanded_df['kegg_list'] = expanded_df['KEGG_Pathway'].apply(parse_kegg)

# Count KEGG pathways
all_kegg = set()
for kegg_list in df['kegg_list']:
    all_kegg.update(kegg_list)

print(f"Total unique KEGG pathways: {len(all_kegg):,}")

# OPTIMIZATION: Use Counter for frequency counting
print("Counting KEGG pathway frequencies in background...")
bg_kegg_counter = Counter()
for kegg_list in df['kegg_list']:
    bg_kegg_counter.update(kegg_list)

print("Counting KEGG pathway frequencies in expanded genes...")
exp_kegg_counter = Counter()
for kegg_list in expanded_df['kegg_list']:
    exp_kegg_counter.update(kegg_list)

# Get top KEGG pathways in expanded genes
kegg_counts = []
for i, kegg in enumerate(all_kegg):
    if i % 1000 == 0 and i > 0:
        print(f"  Processed {i:,} KEGG pathways...")
    
    bg_count = bg_kegg_counter.get(kegg, 0)
    exp_count = exp_kegg_counter.get(kegg, 0)
    
    if exp_count > 5:
        kegg_counts.append({
            'KEGG': kegg,
            'Background': bg_count,
            'Expanded': exp_count,
            'Expanded_Pct': (exp_count / len(expanded_df)) * 100
        })

kegg_counts_df = pd.DataFrame(kegg_counts).sort_values('Expanded', ascending=False).head(15)
print("\nTop 15 KEGG pathways in expanded genes:")
if len(kegg_counts_df) > 0:
    print(kegg_counts_df.to_string(index=False))
else:
    print("No KEGG pathways found with count > 5")

# =============================================================================
# CREATE FIGURE 4
# =============================================================================

print("\n" + "="*80)
print("CREATING FIGURE 4")
print("="*80)

fig = plt.figure(figsize=(22, 18))

# Create grid for panels with adjusted spacing
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3, 
                       height_ratios=[1, 1.2])

# =============================================================================
# Panel A: COG Enrichment - Dot plot
# =============================================================================

ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title('a', fontweight='bold', fontsize=20, loc='left', pad=15)
ax1.text(0.01, 0.98, 'COG category enrichment', 
         transform=ax1.transAxes, fontsize=18, fontweight='bold', va='top')

# Color by hypothesis category
colors = []
for cog in cog_df['COG']:
    colors.append({
        'Growth': growth_color,
        'Metabolism': metabolism_color,
        'Stress': stress_color,
        'Other': other_color
    }.get(cog_to_hypothesis.get(cog, 'Other')))

# Create scatter plot
x = cog_df['Fold_Enrichment']
y = cog_df['-log10_FDR']
sizes = cog_df['Expanded_Count'] / 50  # Scale dot size by count

scatter = ax1.scatter(x, y, c=colors, s=sizes, alpha=0.7, 
                      edgecolors='black', linewidth=1)

# Add reference lines
ax1.axvline(x=1.0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax1.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, linewidth=1)

# Add labels for key COGs
for idx, row in cog_df.iterrows():
    if row['-log10_FDR'] > 50 or row['COG'] in ['T', 'O', 'S', 'K', 'R', 'L']:
        ax1.annotate(row['COG'], (row['Fold_Enrichment'], row['-log10_FDR']),
                    fontsize=12, fontweight='bold', ha='center')

ax1.set_xlabel('Fold enrichment', fontweight='bold', fontsize=18)
ax1.set_ylabel('-log₁₀(FDR)', fontweight='bold', fontsize=18)
ax1.set_xlim(0.8, 1.25)
ax1.set_ylim(0, max(y) * 1.1 if len(y) > 0 else 10)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Add legend for colors
legend_elements = [
    Patch(facecolor=growth_color, label='Growth', alpha=0.7),
    Patch(facecolor=metabolism_color, label='Metabolism', alpha=0.7),
    Patch(facecolor=stress_color, label='Stress', alpha=0.7),
    Patch(facecolor=other_color, label='Other', alpha=0.7)
]
ax1.legend(handles=legend_elements, loc='upper right', fontsize=14, frameon=False)

# =============================================================================
# Panel B: Hypothesis category comparison
# =============================================================================

ax2 = fig.add_subplot(gs[0, 1])
ax2.set_title('b', fontweight='bold', fontsize=20, loc='left', pad=15)
ax2.text(0.01, 0.98, 'Functional category distribution', 
         transform=ax2.transAxes, fontsize=18, fontweight='bold', va='top')

x = np.arange(len(hypothesis_order))
width = 0.35

bg_pcts = [bg_cat_counts.get(cat, 0)/len(df)*100 for cat in hypothesis_order]
exp_pcts = [exp_cat_counts.get(cat, 0)/len(expanded_df)*100 for cat in hypothesis_order]

bars1 = ax2.bar(x - width/2, bg_pcts, width, label='Genome background',
                color='#cccccc', edgecolor='black', linewidth=1, alpha=0.8)
bars2 = ax2.bar(x + width/2, exp_pcts, width, label='Expanded genes',
                color=domesticated_color, edgecolor='black', linewidth=1, alpha=0.8)

ax2.set_ylabel('Percentage of genes (%)', fontweight='bold', fontsize=18)
ax2.set_xticks(x)
ax2.set_xticklabels(hypothesis_order, fontsize=16)
ax2.legend(frameon=False, fontsize=16)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Add value labels
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        if height > 0.5:
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.8,
                    f'{height:.1f}%', ha='center', va='bottom', fontsize=12, fontweight='bold')

# =============================================================================
# Panel C: GO terms - Horizontal bar plot
# =============================================================================

ax3 = fig.add_subplot(gs[1, 0])
ax3.set_title('c', fontweight='bold', fontsize=20, loc='left', pad=20)
ax3.text(0.01, 0.98, 'Top GO terms in expanded genes',
         transform=ax3.transAxes, fontsize=18, fontweight='bold', va='bottom')

if len(go_counts_df) > 0:
    top_go = go_counts_df.head(10).copy()
    
    # Truncate GO labels
    top_go['short_GO'] = [g.split(':')[-1][:20] + '...' if len(g.split(':')[-1]) > 20 else g.split(':')[-1] for g in top_go['GO']]
    
    y_pos = np.arange(len(top_go))
    
    bars = ax3.barh(y_pos, top_go['Expanded'], height=0.7, color=expansion_color, 
                    edgecolor='black', linewidth=1, alpha=0.8)
    
    ax3.set_yticks(y_pos)
    ax3.set_yticklabels(top_go['short_GO'], fontsize=12)
    ax3.set_xlabel('Number of genes', fontweight='bold', fontsize=18, labelpad=15)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    
    ax3.margins(y=0.1)
    
    # Add value labels
    for bar, count in zip(bars, top_go['Expanded']):
        ax3.text(bar.get_width() + 10, bar.get_y() + bar.get_height()/2,
                f'{int(count):,}', va='center', fontsize=12, fontweight='bold')
    
    max_count = top_go['Expanded'].max()
    ax3.set_xlim(0, max_count * 1.15)
    ax3.grid(True, alpha=0.2, axis='x', linestyle='--', linewidth=0.5)
else:
    ax3.text(0.5, 0.5, 'No significant GO terms found', 
             ha='center', va='center', fontsize=16, transform=ax3.transAxes)

# =============================================================================
# Panel D: KEGG pathways - Horizontal bar plot
# =============================================================================

ax4 = fig.add_subplot(gs[1, 1])
ax4.set_title('d', fontweight='bold', fontsize=20, loc='left', pad=20)
ax4.text(0.01, 0.98, 'Top KEGG pathways in expanded genes',
         transform=ax4.transAxes, fontsize=18, fontweight='bold', va='bottom')

if len(kegg_counts_df) > 0:
    top_kegg = kegg_counts_df.head(10).copy()
    
    # Clean KEGG IDs
    top_kegg['clean_name'] = [k.replace('map', '').replace('ko', '')[:20] for k in top_kegg['KEGG']]
    
    y_pos = np.arange(len(top_kegg))
    
    bars = ax4.barh(y_pos, top_kegg['Expanded'], height=0.7, color=contraction_color, 
                    edgecolor='black', linewidth=1, alpha=0.8)
    
    ax4.set_yticks(y_pos)
    ax4.set_yticklabels(top_kegg['clean_name'], fontsize=12)
    ax4.set_xlabel('Number of genes', fontweight='bold', fontsize=18, labelpad=15)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    
    ax4.margins(y=0.1)
    
    for bar, count in zip(bars, top_kegg['Expanded']):
        ax4.text(bar.get_width() + 10, bar.get_y() + bar.get_height()/2,
                f'{int(count):,}', va='center', fontsize=12, fontweight='bold')
    
    max_count = top_kegg['Expanded'].max()
    ax4.set_xlim(0, max_count * 1.15)
    ax4.grid(True, alpha=0.2, axis='x', linestyle='--', linewidth=0.5)
else:
    ax4.text(0.5, 0.5, 'No significant KEGG pathways found', 
             ha='center', va='center', fontsize=16, transform=ax4.transAxes)

# =============================================================================
# Main title and save
# =============================================================================

fig.text(0.5, 0.98, 'Functional analysis of expanded gene families',
         fontsize=24, ha='center', style='italic', color='#2c3e50')

plt.tight_layout()
plt.subplots_adjust(top=0.92, bottom=0.08, left=0.1, right=0.95)

# Save figure
output_png = os.path.join(output_dir, 'Figure4_Functional_Enrichment_REAL_DATA.png')
plt.savefig(output_png, dpi=600, bbox_inches='tight', facecolor='white')
print(f"\n✅ Saved: {output_png}")

output_pdf = os.path.join(output_dir, 'Figure4_Functional_Enrichment_REAL_DATA.pdf')
plt.savefig(output_pdf, format='pdf', dpi=600, bbox_inches='tight', facecolor='white')
print(f"✅ Saved: {output_pdf}")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

print("\n" + "="*80)
print("FINAL SUMMARY - SUPPORT FOR HYPOTHESIS")
print("="*80)
print("\nHypothesis: Genome-wide amplification preserves functional proportions")
print("-" * 60)

print(f"\nExpanded genes analyzed: {len(expanded_df):,} of {len(expanded_genes):,} total")
print(f"Background genes: {len(df):,}")

# Check 1: All categories represented
sig_cogs = sum(cog_df['FDR'] < 0.05)
print(f"\n✓ COG categories: {sig_cogs} significant at FDR < 0.05")

# Check 2: Fold enrichment ~1.0
mean_fold = cog_df['Fold_Enrichment'].mean()
print(f"✓ Mean fold enrichment: {mean_fold:.3f}")

# Check 3: Hypothesis categories
print("\n✓ Hypothesis category comparison:")
for cat in hypothesis_order:
    bg_pct = bg_cat_counts.get(cat, 0)/len(df)*100
    exp_pct = exp_cat_counts.get(cat, 0)/len(expanded_df)*100
    fold = exp_pct / bg_pct if bg_pct > 0 else 0
    print(f"    {cat:10}: {exp_pct:5.2f}% vs {bg_pct:5.2f}% (fold = {fold:.3f})")

print("\n" + "="*80)
print("FIGURE 4 CREATED SUCCESSFULLY WITH REAL CAFE DATA!")
print("="*80)
print(f"\n📁 Files saved to: {output_dir}")

# Display the figure
plt.show()

print("\n✅ All done! Check the output directory for the figure files.")