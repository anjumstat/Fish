# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 10:29:35 2026

@author: H.A.R
"""

# -*- coding: utf-8 -*-
"""
Create Figure 2: Orthogroup Overview with Increased Row Spacing
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality style with larger fonts
plt.style.use('default')
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['figure.dpi'] = 500
plt.rcParams['savefig.dpi'] = 500
plt.rcParams['savefig.bbox'] = 'tight'

# =============================================================================
# Set Output Directory
# =============================================================================

output_dir = r"E:/fish/paper/figures/Figure2"
os.makedirs(output_dir, exist_ok=True)
print(f"Output will be saved to: {output_dir}")

# =============================================================================
# Load OrthoFinder Results
# =============================================================================

# Path to your OrthoFinder results
ortho_path = r"E:/fish/orth/OrthoFinder/Results_Feb11"
orthogroups_file = os.path.join(ortho_path, "Orthogroups", "Orthogroups.tsv")
gene_counts_file = os.path.join(ortho_path, "Orthogroups", "Orthogroups.GeneCount.tsv")

# Check if files exist
if not os.path.exists(orthogroups_file):
    print(f"ERROR: File not found: {orthogroups_file}")
    exit()

# Load orthogroups
print("\nLoading Orthogroups.tsv...")
og = pd.read_csv(orthogroups_file, sep='\t')
print(f"  Found {len(og)} orthogroups")

# Load gene counts
print("Loading Orthogroups.GeneCount.tsv...")
if os.path.exists(gene_counts_file):
    counts = pd.read_csv(gene_counts_file, sep='\t')
    print(f"  Found {len(counts)} orthogroups with counts")
else:
    print(f"  Warning: {gene_counts_file} not found")
    counts = None

# Get species columns
species_cols = [col for col in og.columns if col not in ['Orthogroup', 'Total']]
print(f"  Species ({len(species_cols)} total):")

# Simplify species names for display
species_short = []
for sp in species_cols:
    if 'Carassius' in sp:
        species_short.append('C. auratus')
    elif 'Cyprinus' in sp:
        species_short.append('C. carpio')
    elif 'Larimichthys' in sp:
        species_short.append('L. crocea')
    elif 'Oryzias' in sp:
        species_short.append('O. sinensis')
    elif 'Sinocyclocheilus' in sp and 'anshuiensis' in sp:
        species_short.append('S. anshuiensis')
    elif 'Sinocyclocheilus' in sp and 'grahami' in sp:
        species_short.append('S. grahami')
    else:
        species_short.append(sp[:15])

print(f"\nSimplified species names: {species_short}")

# =============================================================================
# Create main figure with subplots - INCREASED ROW SPACING
# =============================================================================

# Create figure with larger size to accommodate spacing
fig = plt.figure(figsize=(18, 16))  # Increased height from 15 to 16

# Create a grid for all panels with increased vertical spacing
# height_ratios: first row taller, middle row medium, bottom row compact
# hspace: controls space between rows (INCREASED from 0.4 to 0.6)
gs = fig.add_gridspec(3, 3, height_ratios=[1.3, 1.1, 0.7], hspace=0.6, wspace=0.3)

# =============================================================================
# Panel A: Presence/Absence Heatmap (Row 1)
# =============================================================================

ax1 = fig.add_subplot(gs[0, :])
ax1.set_title('A. Orthogroup Presence/Absence Across Species', fontweight='bold', fontsize=20, loc='left', pad=20)

print("Creating presence/absence matrix...")

# Create presence matrix for top 100 orthogroups
n_orthogroups = min(100, len(og))
presence_matrix = []
og_indices = []

for idx in range(n_orthogroups):
    row = og.iloc[idx]
    presence_row = []
    for sp in species_cols:
        if pd.notna(row[sp]) and str(row[sp]).strip() != '' and str(row[sp]).strip().lower() != 'nan':
            presence_row.append(1)
        else:
            presence_row.append(0)
    presence_matrix.append(presence_row)
    og_indices.append(f"OG{idx+1}")

presence_matrix = np.array(presence_matrix)

# Plot heatmap
im = ax1.imshow(presence_matrix, cmap='YlGnBu', aspect='auto', interpolation='nearest')

# Set y-ticks
y_tick_positions = range(0, n_orthogroups, 10)
y_tick_labels = [og_indices[i] for i in y_tick_positions]
ax1.set_yticks(y_tick_positions)
ax1.set_yticklabels(y_tick_labels, fontsize=14, fontweight='bold')

# Set x-ticks
ax1.set_xticks(range(len(species_short)))
ax1.set_xticklabels(species_short, rotation=45, ha='right', fontsize=14, fontweight='bold')

ax1.set_ylabel('Orthogroups (first 100 shown)', fontweight='bold', fontsize=16, labelpad=15)
plt.colorbar(im, ax=ax1, label='Present (1) / Absent (0)')

# Add grid lines
ax1.set_xticks(np.arange(-0.5, len(species_short), 1), minor=True)
ax1.set_yticks(np.arange(-0.5, n_orthogroups, 10), minor=True)
ax1.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# =============================================================================
# Panel B: Orthogroup Size Distribution (Row 2, Col 1)
# =============================================================================

ax2 = fig.add_subplot(gs[1, 0])
ax2.set_title('B. Orthogroup Size Distribution', fontweight='bold', fontsize=18, loc='left', pad=15)

# Calculate orthogroup sizes
og_sizes = []
for idx, row in og.iterrows():
    size = 0
    for sp in species_cols:
        if pd.notna(row[sp]) and str(row[sp]).strip() != '' and str(row[sp]).strip().lower() != 'nan':
            genes = str(row[sp]).split(', ')
            size += len(genes)
    og_sizes.append(size)

# Plot histogram
ax2.hist(og_sizes, bins=50, color='steelblue', edgecolor='black', alpha=0.7, log=True)
ax2.set_xlabel('Orthogroup Size (genes)', fontweight='bold', fontsize=14, labelpad=10)
ax2.set_ylabel('Number of Orthogroups (log scale)', fontweight='bold', fontsize=13, labelpad=10)
ax2.set_xscale('log')
ax2.grid(True, alpha=0.3, linestyle='--')

ax2.tick_params(axis='both', which='major', labelsize=14, width=1.5)
for label in ax2.get_xticklabels() + ax2.get_yticklabels():
    label.set_fontweight('bold')

# Add statistics
mean_size = np.mean(og_sizes)
median_size = np.median(og_sizes)
max_size = np.max(og_sizes)
ax2.text(0.7, 0.95, f'Mean: {mean_size:.1f}\nMedian: {median_size:.0f}\nMax: {max_size}',
         transform=ax2.transAxes, fontsize=14, fontweight='bold', verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5, edgecolor='black'))

# =============================================================================
# Panel C: Single vs Multi-copy Orthogroups (Row 2, Col 2)
# =============================================================================

ax3 = fig.add_subplot(gs[1, 1])
ax3.set_title('C. Orthogroups Present in All Species', fontweight='bold', fontsize=18, loc='center', pad=15)

# Classify orthogroups
single_copy_all = 0
multi_copy_all = 0
single_copy_any = 0
multi_copy_any = 0

for idx, row in og.iterrows():
    present_in_all = True
    has_multi = False
    total_genes = 0
    
    for sp in species_cols:
        if pd.notna(row[sp]) and str(row[sp]).strip() != '' and str(row[sp]).strip().lower() != 'nan':
            genes = str(row[sp]).split(', ')
            gene_count = len(genes)
            total_genes += gene_count
            
            if gene_count > 1:
                has_multi = True
        else:
            present_in_all = False
    
    if total_genes > 0:
        if has_multi:
            multi_copy_any += 1
        else:
            single_copy_any += 1
    
    if present_in_all:
        if has_multi:
            multi_copy_all += 1
        else:
            single_copy_all += 1

# Create pie chart
labels = ['Single-copy\n(all 6 species)', 'Multi-copy\n(all 6 species)']
sizes = [single_copy_all, multi_copy_all]

if sum(sizes) > 0:
    colors = ['#66c2a5', '#fc8d62']
    explode = (0.05, 0)
    
    wedges, texts, autotexts = ax3.pie(sizes, explode=explode, labels=labels, colors=colors,
                                        autopct='%1.1f%%', startangle=90,
                                        textprops={'fontsize': 14, 'fontweight': 'bold'})
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(14)
else:
    ax3.text(0.5, 0.5, 'No orthogroups\npresent in all species', 
             ha='center', va='center', transform=ax3.transAxes, 
             fontweight='bold', fontsize=14)

# =============================================================================
# Panel D: Per-Species Gene Distribution (Row 2, Col 3)
# =============================================================================

ax4 = fig.add_subplot(gs[1, 2])
ax4.set_title('D. Genes per Species', fontweight='bold', fontsize=16, loc='left', pad=15)

# Calculate gene counts per species
species_gene_counts = []
for sp in species_cols:
    count = 0
    for idx, row in og.iterrows():
        if pd.notna(row[sp]) and str(row[sp]).strip() != '' and str(row[sp]).strip().lower() != 'nan':
            genes = str(row[sp]).split(', ')
            count += len(genes)
    species_gene_counts.append(count)

# Create bar plot
x_pos = np.arange(len(species_short))
bars = ax4.bar(x_pos, species_gene_counts, color='#8da0cb', edgecolor='black', alpha=0.8, linewidth=1.5)
ax4.set_xticks(x_pos)
ax4.set_xticklabels(species_short, rotation=45, ha='right', fontsize=14, fontweight='bold')
ax4.set_ylabel('Number of Genes', fontweight='bold', fontsize=14, labelpad=10)

ax4.tick_params(axis='y', which='major', labelsize=14, width=1.5)
for label in ax4.get_yticklabels():
    label.set_fontweight('bold')

# Add value labels
for bar, count in zip(bars, species_gene_counts):
    height = bar.get_height()
    ax4.text(bar.get_x() + bar.get_width()/2., height,
             f'{count:,}', ha='center', va='bottom', fontsize=10, 
             fontweight='bold', rotation=0)

ax4.grid(True, alpha=0.3, axis='y', linestyle='--', linewidth=1)

# =============================================================================
# Panel E: Summary Statistics Table (Row 3)
# =============================================================================

ax5 = fig.add_subplot(gs[2, :])
ax5.axis('off')
ax5.set_title('E. Summary Statistics', fontweight='bold', fontsize=18, loc='left', y=0.98, pad=20)

# Calculate statistics
total_genes = sum(species_gene_counts)

# Count orthogroups in domesticated vs wild
domesticated_cols = [col for col in species_cols if 'Carassius' in col or 'Cyprinus' in col]
wild_cols = [col for col in species_cols if 'Sinocyclocheilus' in col]

domesticated_only = 0
wild_only = 0
shared_dom_wild = 0

for idx, row in og.iterrows():
    in_dom = False
    in_wild = False
    
    for sp in domesticated_cols:
        if pd.notna(row[sp]) and str(row[sp]).strip() != '' and str(row[sp]).strip().lower() != 'nan':
            in_dom = True
            break
    
    for sp in wild_cols:
        if pd.notna(row[sp]) and str(row[sp]).strip() != '' and str(row[sp]).strip().lower() != 'nan':
            in_wild = True
            break
    
    if in_dom and not in_wild:
        domesticated_only += 1
    elif not in_dom and in_wild:
        wild_only += 1
    elif in_dom and in_wild:
        shared_dom_wild += 1

stats_data = [
    ['Total Orthogroups', f'{len(og):,}'],
    ['Total Genes', f'{total_genes:,}'],
    ['Single-copy OGs (any species)', f'{single_copy_any:,}'],
    ['Multi-copy OGs (any species)', f'{multi_copy_any:,}'],
    ['Single-copy (all 6 species)', f'{single_copy_all:,}'],
    ['Multi-copy (all 6 species)', f'{multi_copy_all:,}'],
    ['Domesticated-specific OGs', f'{domesticated_only:,}'],
    ['Wild-specific OGs', f'{wild_only:,}'],
    ['Shared OGs (domesticated & wild)', f'{shared_dom_wild:,}'],
    ['Mean Orthogroup Size', f'{mean_size:.2f}'],
    ['Median Orthogroup Size', f'{median_size:.0f}'],
]

# Create table
table = ax5.table(cellText=stats_data, colLabels=['Metric', 'Value'],
                  loc='center', cellLoc='left', colWidths=[0.4, 0.15])
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1, 1.8)

# Style the table
for (i, j), cell in table.get_celld().items():
    if i == 0:
        cell.set_facecolor('#2c3e50')
        cell.set_text_props(color='white', fontweight='bold', fontsize=12)
    else:
        cell.set_facecolor('#f8f9fa' if i % 2 == 0 else '#e9ecef')
        cell.set_text_props(fontweight='bold' if j == 0 else 'normal')
    if j == 1 and i > 0:
        cell.set_text_props(ha='right', fontweight='bold')
    cell.set_edgecolor('black')
    cell.set_linewidth(1)

# =============================================================================
# Add visual separators between rows (optional)
# =============================================================================

# Get the positions of the subplots to add horizontal lines
bbox_a = ax1.get_position()
bbox_b = ax2.get_position()
bbox_e = ax5.get_position()

# Add horizontal line between Row 1 and Row 2
line1 = plt.Line2D([0.02, 0.98], [bbox_a.y0 - 0.02, bbox_a.y0 - 0.02], 
                   transform=fig.transFigure, color='gray', linestyle='--', linewidth=1, alpha=0.5)
fig.add_artist(line1)

# Add horizontal line between Row 2 and Row 3
line2 = plt.Line2D([0.02, 0.98], [bbox_b.y0 - 0.02, bbox_b.y0 - 0.02], 
                   transform=fig.transFigure, color='gray', linestyle='--', linewidth=1, alpha=0.5)
fig.add_artist(line2)

# =============================================================================
# Final formatting and save
# =============================================================================

plt.tight_layout()

# Save figure
output_file = os.path.join(output_dir, 'figure2_orthogroup_overview.png')
plt.savefig(output_file, dpi=500, bbox_inches='tight', facecolor='white')
print(f"\n✓ Figure saved to: {output_file}")

pdf_file = os.path.join(output_dir, 'figure2_orthogroup_overview.pdf')
plt.savefig(pdf_file, format='pdf', bbox_inches='tight', facecolor='white')
print(f"✓ PDF saved to: {pdf_file}")

svg_file = os.path.join(output_dir, 'figure2_orthogroup_overview.svg')
plt.savefig(svg_file, format='svg', bbox_inches='tight', facecolor='white')
print(f"✓ SVG saved to: {svg_file}")

print("\n" + "="*60)
print("FIGURE 2 CREATED SUCCESSFULLY")
print("="*60)

plt.show()