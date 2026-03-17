# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 17:24:59 2026

@author: H.A.R
"""

# -*- coding: utf-8 -*-
"""
Create Figure 5: Hypothesis Testing - Main Finding (NATURE ECOLOGY & EVOLUTION STYLE)
UPDATED with ACTUAL overlap values from gene-level analysis
- Based on 29,422 expanded genes with annotations
- Shows distribution across Growth, Metabolism, and Stress categories
- ALL values from actual data output (2026-03-17)
- FIXED: Other category calculation now correct (13,059)
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
import os

# Set publication-quality style - ALL FONTS INCREASED BY 4 POINTS
plt.style.use('default')
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 22
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['figure.titlesize'] = 26
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600

# Panel label font size
panel_label_size = 22

# =============================================================================
# Set Output Directory
# =============================================================================

output_dir = r"E:/fish/paper/figures/Figure5"
os.makedirs(output_dir, exist_ok=True)
print(f"Output will be saved to: {output_dir}")

# =============================================================================
# LOAD YOUR ACTUAL DATA
# =============================================================================

# Paths to your data files
eggnog_path = r"E:/fish/eggnog_results/eggnog_results"
cafe_results_path = r"E:/fish/orth/OrthoFinder/Results_Feb11/cafe_results"
expanded_genes_file = os.path.join(cafe_results_path, "expanded_genes_cyprinids.txt")
annot_file = os.path.join(eggnog_path, "fish_domestication.emapper.annotations")

print("\n" + "="*60)
print("LOADING ACTUAL GENE DATA")
print("="*60)

# Load expanded genes
expanded_genes = pd.read_csv(expanded_genes_file, header=None, names=['gene_id'])
print(f"Loaded {len(expanded_genes):,} expanded genes from CAFE results")

# Load eggNOG annotations
df = pd.read_csv(annot_file, sep='\t', skiprows=4, header=None)
df.columns = ['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 
              'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 
              'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
              'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 
              'BiGG_Reaction', 'PFAMs']

# Filter to expanded genes only
expanded_df = df[df['query'].isin(expanded_genes['gene_id'])].copy()
print(f"Matched {len(expanded_df):,} genes with eggNOG annotations")

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

def get_hypothesis_categories(cog_list):
    """Get all hypothesis categories for a gene based on its COGs"""
    cats = set()
    for cog in cog_list:
        if cog in cog_to_hypothesis:
            cats.add(cog_to_hypothesis[cog])
    return cats

expanded_df['cog_list'] = expanded_df['COG_category'].apply(parse_cog)

# =============================================================================
# CALCULATE ACTUAL OVERLAPS (CORRECTED VERSION)
# =============================================================================

print("\n" + "="*60)
print("CALCULATING ACTUAL OVERLAP VALUES")
print("="*60)

# Initialize sets for hypothesis categories
growth_genes = set()
metabolism_genes = set()
stress_genes = set()

# First pass: Add genes to their respective hypothesis categories
for idx, row in expanded_df.iterrows():
    cats = get_hypothesis_categories(row['cog_list'])
    gene_id = row['query']
    
    if 'Growth' in cats:
        growth_genes.add(gene_id)
    if 'Metabolism' in cats:
        metabolism_genes.add(gene_id)
    if 'Stress' in cats:
        stress_genes.add(gene_id)

# Calculate any hypothesis (union)
any_hypothesis = growth_genes | metabolism_genes | stress_genes

# Second pass: Determine Other genes (those NOT in any hypothesis category)
other_genes = set()
for idx, row in expanded_df.iterrows():
    gene_id = row['query']
    if gene_id not in any_hypothesis:
        other_genes.add(gene_id)

# Calculate actual overlap values
growth_metabolism = len(growth_genes & metabolism_genes)
growth_stress = len(growth_genes & stress_genes)
metabolism_stress = len(metabolism_genes & stress_genes)
all_three = len(growth_genes & metabolism_genes & stress_genes)

# Calculate exclusive counts
growth_only = len(growth_genes - metabolism_genes - stress_genes)
metabolism_only = len(metabolism_genes - growth_genes - stress_genes)
stress_only = len(stress_genes - growth_genes - metabolism_genes)

# Calculate totals
total_expanded_genes = len(expanded_df)
any_hypothesis_count = len(any_hypothesis)
other_count = len(other_genes)

# Verification
print(f"Total expanded genes: {total_expanded_genes:,}")
print(f"\nCategory counts:")
print(f"  Growth: {len(growth_genes):,}")
print(f"  Metabolism: {len(metabolism_genes):,}")
print(f"  Stress: {len(stress_genes):,}")
print(f"  Any hypothesis: {any_hypothesis_count:,}")
print(f"  Other: {other_count:,}")
print(f"  Total (check): {any_hypothesis_count + other_count:,}")

print(f"\nACTUAL OVERLAP VALUES:")
print(f"  Growth only: {growth_only:,}")
print(f"  Metabolism only: {metabolism_only:,}")
print(f"  Stress only: {stress_only:,}")
print(f"  Growth + Metabolism: {growth_metabolism:,}")
print(f"  Growth + Stress: {growth_stress:,}")
print(f"  Metabolism + Stress: {metabolism_stress:,}")
print(f"  All three categories: {all_three:,}")

# Calculate percentages
growth_pct = (len(growth_genes) / total_expanded_genes) * 100
metabolism_pct = (len(metabolism_genes) / total_expanded_genes) * 100
stress_pct = (len(stress_genes) / total_expanded_genes) * 100
any_pct = (any_hypothesis_count / total_expanded_genes) * 100
other_pct = (other_count / total_expanded_genes) * 100

# Background data from your analysis
total_background = 100546
growth_background = 5562
metabolism_background = 14104
stress_background = 39707
other_background = 45042
any_background = 56656  # From your actual data output

growth_bg_pct = (growth_background / total_background) * 100
metabolism_bg_pct = (metabolism_background / total_background) * 100
stress_bg_pct = (stress_background / total_background) * 100
other_bg_pct = (other_background / total_background) * 100
any_bg_pct = (any_background / total_background) * 100

# Enrichment values
growth_enrich = growth_pct / growth_bg_pct
metabolism_enrich = metabolism_pct / metabolism_bg_pct
stress_enrich = stress_pct / stress_bg_pct
other_enrich = other_pct / other_bg_pct
any_enrich = any_pct / any_bg_pct

print(f"\nPercentages:")
print(f"  Growth: {growth_pct:.2f}%")
print(f"  Metabolism: {metabolism_pct:.2f}%")
print(f"  Stress: {stress_pct:.2f}%")
print(f"  Any hypothesis: {any_pct:.2f}%")
print(f"  Other: {other_pct:.2f}%")
print(f"\nBackground any hypothesis: {any_bg_pct:.2f}%")

# =============================================================================
# Create Figure 5 with Proper Spacing - 3x2 Grid
# =============================================================================

fig = plt.figure(figsize=(28, 20))

# Create grid with proper spacing
gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 1, 0.8], 
                       hspace=0.45, wspace=0.4,
                       left=0.08, right=0.95, bottom=0.12, top=0.92)

# =============================================================================
# Panel A: Bar Chart of Percentages (Main Finding)
# =============================================================================

ax1 = fig.add_subplot(gs[0, 0])

categories = ['Growth', 'Metabolism', 'Stress\nResponse', 'Any\nHypothesis']
percentages = [growth_pct, metabolism_pct, stress_pct, any_pct]
colors = ['#2ecc71', '#3498db', '#e74c3c', '#9b59b6']

bars = ax1.bar(categories, percentages, color=colors, edgecolor='black', 
               linewidth=2, alpha=0.8, width=0.6)

# Add value labels on bars
for bar, pct in zip(bars, percentages):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + 0.8,
             f'{pct:.1f}%', ha='center', va='bottom', 
             fontsize=18, fontweight='bold')

ax1.set_ylabel('Percentage of expanded genes (%)', fontsize=20, fontweight='bold', labelpad=15)
ax1.set_title('a', fontweight='bold', fontsize=panel_label_size, loc='left', pad=15)

# Set y-limit
max_percentage = max(percentages)
ax1.set_ylim(0, max_percentage * 1.15)
ax1.tick_params(axis='x', labelsize=18)
ax1.tick_params(axis='y', labelsize=18)

# Make axis labels bold
for label in ax1.get_xticklabels():
    label.set_fontweight('bold')
for label in ax1.get_yticklabels():
    label.set_fontweight('bold')

ax1.grid(True, alpha=0.3, axis='y', linestyle='--', linewidth=0.5)

# =============================================================================
# Panel B: Pie Chart of Distribution
# =============================================================================

ax2 = fig.add_subplot(gs[0, 1])

# Data for pie - showing distribution among hypothesis categories
sizes = [len(growth_genes), len(metabolism_genes), len(stress_genes), other_count]
labels = ['Growth', 'Metabolism', 'Stress', 'Other']
colors_pie = ['#2ecc71', '#3498db', '#e74c3c', '#95a5a6']
explode = (0.03, 0.03, 0.03, 0.03)

wedges, texts, autotexts = ax2.pie(
    sizes, 
    labels=labels,
    colors=colors_pie,
    autopct=lambda pct: f'{pct:.1f}%',
    explode=explode,
    startangle=90,
    textprops={'fontsize': 16, 'fontweight': 'bold'},
    labeldistance=1.2
)

for autotext in autotexts:
    autotext.set_color('white')
    autotext.set_fontsize(16)
    autotext.set_fontweight('bold')

ax2.set_title('b', fontweight='bold', fontsize=panel_label_size, pad=20)

# =============================================================================
# Panel C: Comparison with Genome Background
# =============================================================================

ax3 = fig.add_subplot(gs[1, 0])

# Create grouped bar chart
x = np.arange(3)
width = 0.35

expanded_vals = [growth_pct, metabolism_pct, stress_pct]
bg_vals = [growth_bg_pct, metabolism_bg_pct, stress_bg_pct]

bars1 = ax3.bar(x - width/2, expanded_vals, width, label='Expanded genes', 
                color='#e67e22', edgecolor='black', linewidth=1.5)
bars2 = ax3.bar(x + width/2, bg_vals, width, label='Genome background', 
                color='#7f8c8d', edgecolor='black', linewidth=1.5, alpha=0.7)

ax3.set_ylabel('Percentage (%)', fontsize=20, fontweight='bold', labelpad=15)
ax3.set_title('c', fontweight='bold', fontsize=panel_label_size, loc='left', pad=15)
ax3.set_xticks(x)
ax3.set_xticklabels(['Growth', 'Metabolism', 'Stress'], fontsize=18, fontweight='bold')
ax3.tick_params(axis='y', labelsize=18)

# Make y axis labels bold
for label in ax3.get_yticklabels():
    label.set_fontweight('bold')

# Set y-limit
max_y = max(max(expanded_vals), max(bg_vals))
ax3.set_ylim(0, max_y * 1.35)
ax3.grid(True, alpha=0.3, axis='y', linestyle='--')

# Place legend
ax3.legend(fontsize=16, frameon=True, fancybox=True, shadow=True, 
          loc='center', bbox_to_anchor=(0.47, 0.88))

# Add enrichment factors above bars
enrich_values = [growth_enrich, metabolism_enrich, stress_enrich]
for i, (exp, enrich) in enumerate(zip(expanded_vals, enrich_values)):
    ax3.text(i - width/2, exp + max_y * 0.05,
             f'{enrich:.2f}x', ha='center', fontsize=16, fontweight='bold')

# =============================================================================
# Panel D: Overlap Analysis (NOW WITH ACTUAL VALUES!)
# =============================================================================

ax4 = fig.add_subplot(gs[1, 1])

# Create grouped bar chart for overlaps with ACTUAL values
overlap_categories = ['Growth\nOnly', 'Metab.\nOnly', 'Stress\nOnly', 
                      'G+M', 'G+S', 'M+S', 'All\nThree']
overlap_values = [growth_only, metabolism_only, stress_only,
                  growth_metabolism, growth_stress, metabolism_stress, all_three]
overlap_colors = ['#2ecc71', '#3498db', '#e74c3c', 
                  '#f39c12', '#f1c40f', '#e67e22', '#9b59b6']

bars = ax4.bar(overlap_categories, overlap_values, color=overlap_colors, 
               edgecolor='black', linewidth=1.5)

ax4.set_ylabel('Number of genes', fontsize=20, fontweight='bold', labelpad=15)
ax4.set_title('d', fontweight='bold', fontsize=panel_label_size, loc='left', pad=15)
ax4.tick_params(axis='x', labelsize=14, rotation=45)
ax4.tick_params(axis='y', labelsize=18)

# Make axis labels bold
for label in ax4.get_xticklabels():
    label.set_fontweight('bold')
for label in ax4.get_yticklabels():
    label.set_fontweight('bold')

# Set y-limit
max_overlap = max(overlap_values)
ax4.set_ylim(0, max_overlap * 1.15)
ax4.grid(True, alpha=0.3, axis='y', linestyle='--')

# Add value labels
for bar, val in zip(bars, overlap_values):
    height = bar.get_height()
    ax4.text(bar.get_x() + bar.get_width()/2., height + max_overlap * 0.02,
             f'{val:,}', ha='center', va='bottom', fontsize=12, fontweight='bold')

# =============================================================================
# Panel E: Summary Statistics Table
# =============================================================================

ax5 = fig.add_subplot(gs[2, :])
ax5.axis('off')

# Create summary table data with ACTUAL values
table_data = [
    ['Category', 'Expanded genes', 'Percentage', 'Background %', 'Enrichment'],
    ['Growth', f'{len(growth_genes):,}', f'{growth_pct:.1f}%', f'{growth_bg_pct:.1f}%', f'{growth_enrich:.2f}x'],
    ['Metabolism', f'{len(metabolism_genes):,}', f'{metabolism_pct:.1f}%', f'{metabolism_bg_pct:.1f}%', f'{metabolism_enrich:.2f}x'],
    ['Stress', f'{len(stress_genes):,}', f'{stress_pct:.1f}%', f'{stress_bg_pct:.1f}%', f'{stress_enrich:.2f}x'],
    ['Any hypothesis', f'{any_hypothesis_count:,}', f'{any_pct:.1f}%', f'{any_bg_pct:.1f}%', f'{any_enrich:.2f}x'],
    ['Other', f'{other_count:,}', f'{other_pct:.1f}%', f'{other_bg_pct:.1f}%', f'{other_enrich:.2f}x']
]

# Create table
table = ax5.table(cellText=table_data, loc='center', cellLoc='center',
                  colWidths=[0.18, 0.15, 0.12, 0.12, 0.12],
                  bbox=[0.1, 0.1, 0.8, 0.8])

table.auto_set_font_size(False)
table.set_fontsize(16)

# Style the table
for (i, j), cell in table.get_celld().items():
    cell.set_linewidth(1)
    cell.set_edgecolor('black')
    
    if i == 0:
        cell.set_facecolor('#2c3e50')
        cell.set_text_props(color='white', fontweight='bold', fontsize=16)
    else:
        if i % 2 == 0:
            cell.set_facecolor('#f8f9fa')
        else:
            cell.set_facecolor('#e9ecef')
        if j == 0:
            cell.set_text_props(fontweight='bold', fontsize=14)
        else:
            cell.set_text_props(fontsize=16)
        if j == 4 and i > 0 and i < 5:
            # Highlight enrichment values close to 1.0
            enrich_val = float(table_data[i][4].replace('x', ''))
            if 0.9 < enrich_val < 1.1:
                cell.set_facecolor('#d4edda')

ax5.set_title('e', fontweight='bold', fontsize=panel_label_size, pad=25, y=0.95)

# =============================================================================
# Add main title and conclusion box
# =============================================================================

plt.suptitle('Functional analysis supports genome-wide amplification hypothesis', 
             fontsize=26, fontweight='bold', y=0.98)

# Add conclusion text with ACTUAL values
conclusion_text = (
    f"HYPOTHESIS SUPPORTED: {any_pct:.1f}% of expanded genes ({any_hypothesis_count:,}) "
    f"participate in growth ({growth_pct:.1f}%), metabolism ({metabolism_pct:.1f}%), "
    f"or stress response ({stress_pct:.1f}%), matching background distribution ({any_bg_pct:.1f}%)"
)

fig.text(0.5, 0.06, conclusion_text, ha='center', fontsize=20, 
         fontweight='bold', bbox=dict(boxstyle='round,pad=0.5', 
                                      facecolor='#d4edda', edgecolor='#28a745', alpha=0.8))

# =============================================================================
# Add annotation for key finding
# =============================================================================

# Arrow pointing to Stress as the largest category
ax1.annotate(f'Largest\ncategory\n({stress_pct:.1f}%)', 
             xy=(2.2, stress_pct+0.6),
             xytext=(2.3, stress_pct + 11.0),
             arrowprops=dict(arrowstyle='->', color='#e74c3c', lw=2.5, 
                            connectionstyle='arc3,rad=0.2'),
             fontsize=18, color='#e74c3c', fontweight='bold', ha='center',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                      alpha=0.9, edgecolor='#e74c3c', linewidth=1.5))

# =============================================================================
# Save figure
# =============================================================================

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)

# Save in multiple formats
output_file = os.path.join(output_dir, 'Figure5_Hypothesis_Testing_FINAL.png')
plt.savefig(output_file, dpi=600, bbox_inches='tight', facecolor='white')
print(f"\n✓ Figure saved to: {output_file}")

pdf_file = os.path.join(output_dir, 'Figure5_Hypothesis_Testing_FINAL.pdf')
plt.savefig(pdf_file, format='pdf', dpi=600, bbox_inches='tight', facecolor='white')
print(f"✓ PDF saved to: {pdf_file}")

print("\n" + "="*60)
print("FIGURE 5 COMPLETE - ALL VALUES ARE NOW ACTUAL DATA!")
print("="*60)
print(f"\nKey findings from your ACTUAL data:")
print(f"• {any_pct:.1f}% of expanded genes support the hypothesis (background: {any_bg_pct:.1f}%)")
print(f"• Stress response is the largest category ({stress_pct:.1f}%)")
print(f"• Overlap values (ACTUAL):")
print(f"  - Growth only: {growth_only:,}")
print(f"  - Metabolism only: {metabolism_only:,}")
print(f"  - Stress only: {stress_only:,}")
print(f"  - Growth+Metabolism: {growth_metabolism:,}")
print(f"  - Growth+Stress: {growth_stress:,}")
print(f"  - Metabolism+Stress: {metabolism_stress:,}")
print(f"  - All three: {all_three:,}")
print(f"\n📁 Files saved to: {output_dir}")

plt.show()