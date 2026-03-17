# -*- coding: utf-8 -*-
"""
Create publication-quality figure from CAFE results
Visualizing domesticated vs wild expansion/contraction patterns
REDESIGNED LAYOUT - Panels A,B in row1; C,E in row2; D in row3
FONT SIZES INCREASED BY 8 POINTS TOTAL
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# =============================================================================
# Your data from Table 5
# =============================================================================

# Domesticated vs Wild comparisons
dom_expansions = 35786
dom_contractions = 4384
wild_expansions = 12730
wild_contractions = 30213

# Significance data
total_families = 24965
sig_p05 = 6388
sig_p01 = 4305
sig_p001 = 2049

# Species-specific data (from your clade results)
species_data = {
    'C. auratus': {'expansions': 18215, 'contractions': 2163, 'group': 'Domesticated'},
    'C. carpio': {'expansions': 17571, 'contractions': 2221, 'group': 'Domesticated'},
    'S. anshuiensis': {'expansions': 5516, 'contractions': 16097, 'group': 'Wild'},
    'S. grahami': {'expansions': 7214, 'contractions': 14116, 'group': 'Wild'},
    'L. crocea': {'expansions': 14770, 'contractions': 5863, 'group': 'Outgroup'},
    'O. sinensis': {'expansions': 13883, 'contractions': 6001, 'group': 'Outgroup'}
}

# =============================================================================
# Create output directory
# =============================================================================

output_dir = "E:/fish/paper/figures/Figure3"
os.makedirs(output_dir, exist_ok=True)

# =============================================================================
# Set publication style - Nature Ecology & Evolution (EXTRA LARGE FONTS)
# =============================================================================

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18  # Increased from 14 to 18 (+4 more)
plt.rcParams['axes.labelsize'] = 19  # Increased from 15 to 19 (+4 more)
plt.rcParams['axes.titlesize'] = 20  # Increased from 16 to 20 (+4 more)
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['legend.fontsize'] = 17  # Increased from 13 to 17 (+4 more)
plt.rcParams['figure.titlesize'] = 22  # Increased from 18 to 22 (+4 more)
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['xtick.labelsize'] = 17  # Increased from 13 to 17 (+4 more)
plt.rcParams['ytick.labelsize'] = 17  # Increased from 13 to 17 (+4 more)

# Nature-style colors - Colorblind friendly
expansion_color = '#1b9e77'  # Teal green
contraction_color = '#d95f02'  # Orange-red
domesticated_color = '#d95f02'  # Orange-red
wild_color = '#7570b3'  # Purple
outgroup_color = '#666666'  # Gray

# =============================================================================
# Create figure with redesigned layout
# 3 rows: Row1 (A,B), Row2 (C,E), Row3 (D full width)
# =============================================================================

fig = plt.figure(figsize=(14, 18))  # Increased figure size to accommodate even larger fonts

# Create grid with 3 rows
# height_ratios: row1 (A,B), row2 (C,E), row3 (D full width)
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.5], hspace=0.45, wspace=0.5)

# =============================================================================
# Row 1: Panel A (left) and Panel B (right)
# =============================================================================

# Panel A: Domesticated vs Wild bar chart (top left)
ax1 = fig.add_subplot(gs[0, 0])

x = np.arange(2)
width = 0.35

bars1 = ax1.bar(x - width/2, [dom_expansions, wild_expansions], width, 
                label='Expansions', color=expansion_color, edgecolor='black', linewidth=1, alpha=0.9)
bars2 = ax1.bar(x + width/2, [dom_contractions, wild_contractions], width, 
                label='Contractions', color=contraction_color, edgecolor='black', linewidth=1, alpha=0.9)

ax1.set_ylabel('Number of gene families', fontweight='bold', labelpad=20)
ax1.set_title('a', fontweight='bold', loc='left', pad=20)
ax1.set_xticks(x)
ax1.set_xticklabels(['Domesticated', 'Wild'], fontsize=18)
ax1.legend(frameon=False, loc='upper right', fontsize=17)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Set y-limit
max_height = max(dom_expansions, wild_expansions, dom_contractions, wild_contractions)
ax1.set_ylim(0, max_height * 1.2)

# Add value labels with increased font size
for bar in bars1:
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + max_height*0.03,
             f'{int(height/1000)}k', ha='center', va='bottom', fontweight='bold', fontsize=16)

for bar in bars2:
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + max_height*0.03,
             f'{int(height/1000)}k', ha='center', va='bottom', fontweight='bold', fontsize=16)


# Panel B: Net change (top right)
ax2 = fig.add_subplot(gs[0, 1])

net_dom = dom_expansions - dom_contractions
net_wild = wild_expansions - wild_contractions

bars = ax2.bar(['Domesticated', 'Wild'], [net_dom, net_wild], 
               color=[domesticated_color, wild_color], edgecolor='black', linewidth=1, alpha=0.9)
ax2.axhline(y=0, color='black', linewidth=1, linestyle='-', alpha=0.3)
ax2.set_ylabel('Net change', fontweight='bold', labelpad=20)
ax2.set_title('b', fontweight='bold', loc='left', pad=20)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.tick_params(axis='x', labelsize=18)

# Set y-limit
max_net = max(abs(net_dom), abs(net_wild))
ax2.set_ylim(-max_net * 1.3, max_net * 1.3)

# Add value labels inside bars with increased font size
for bar in bars:
    height = bar.get_height()
    if height > 0:
        ax2.text(bar.get_x() + bar.get_width()/2., height * 0.3,
                 f'+{int(height/1000)}k', ha='center', va='center', 
                 fontweight='bold', color='white', fontsize=16,
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.6))
    else:
        ax2.text(bar.get_x() + bar.get_width()/2., height * 0.7,
                 f'{int(height/1000)}k', ha='center', va='center', 
                 fontweight='bold', color='white', fontsize=16,
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.6))


# =============================================================================
# Row 2: Panel C (left) and Panel E (right)
# =============================================================================

# Panel C: Ratio comparison (middle left)
ax3 = fig.add_subplot(gs[1, 0])

exp_ratio = dom_expansions / wild_expansions
cont_ratio = dom_contractions / wild_contractions

x = np.arange(2)
width = 0.5

# Create bars
bars = ax3.bar(x, [exp_ratio, cont_ratio], width, 
               color=[expansion_color, contraction_color], 
               edgecolor='black', linewidth=1, alpha=0.9)

# Add reference line at y=1 (equal ratio)
ax3.axhline(y=1, color='#2c3e50', linewidth=1.5, linestyle='--', alpha=0.7, zorder=0)

ax3.set_ylabel('Ratio (domesticated/wild)', fontweight='bold', labelpad=20)
ax3.set_title('c', fontweight='bold', loc='left', pad=20)
ax3.set_xticks(x)
ax3.set_xticklabels(['Expansions', 'Contractions'], fontsize=18)

# Apply log scale
ax3.set_yscale('log')
ax3.set_ylim(0.08, 4)
ax3.set_yticks([0.1, 0.2, 0.5, 1, 2, 4])
ax3.set_yticklabels(['0.1', '0.2', '0.5', '1', '2', '4'], fontsize=17)

ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.grid(axis='y', alpha=0.2, linestyle='-', linewidth=0.5, color='#cccccc')

# Add value labels above bars with increased font size
for bar, ratio in zip(bars, [exp_ratio, cont_ratio]):
    height = bar.get_height()
    ax3.text(bar.get_x() + bar.get_width()/2., height * 1.2,
             f'{ratio:.2f}', ha='center', va='bottom', 
             fontweight='bold', fontsize=17, color='black')

# Add fold-difference text inside bars with increased font size
if exp_ratio > 1:
    ax3.text(0, exp_ratio * 0.5, f'{exp_ratio:.1f}×', 
             ha='center', va='center', fontsize=17, color='white', fontweight='bold')
if cont_ratio < 1:
    ax3.text(1, cont_ratio * 0.5, f'{(1/cont_ratio):.1f}×', 
             ha='center', va='center', fontsize=17, color='white', fontweight='bold')


# Panel E: Statistical significance (middle right)
ax5 = fig.add_subplot(gs[1, 1])

categories = ['p < 0.05', 'p < 0.01', 'p < 0.001']
counts = [sig_p05, sig_p01, sig_p001]
percentages = [sig_p05/total_families*100, sig_p01/total_families*100, sig_p001/total_families*100]

colors = plt.cm.Reds(np.linspace(0.4, 0.8, 3))

bars = ax5.bar(categories, counts, color=colors, edgecolor='black', linewidth=1, alpha=0.9)
ax5.set_ylabel('Number of families', fontweight='bold', labelpad=20)
ax5.set_title('d', fontweight='bold', loc='left', pad=20)
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)

max_count = max(counts)
ax5.set_ylim(0, max_count * 1.3)

# Rotate x-axis labels slightly with increased font size
ax5.set_xticklabels(categories, rotation=10, ha='right', fontsize=17)

# Add value labels with increased font size
for bar, pct in zip(bars, percentages):
    height = bar.get_height()
    ax5.text(bar.get_x() + bar.get_width()/2., height + max_count*0.04,
             f'{int(height/1000)}k\n({pct:.1f}%)', ha='center', va='bottom', 
             fontsize=15, fontweight='bold', linespacing=1.3)

# Add total families text with increased font size
ax5.text(0.98, 0.98, f'Total: {total_families:,}', 
         transform=ax5.transAxes, ha='right', va='top', fontsize=15, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8f9fa', alpha=0.8, edgecolor='none'))


# =============================================================================
# Row 3: Panel D - Species breakdown (full width)
# =============================================================================

ax4 = fig.add_subplot(gs[2, :])

species_list = list(species_data.keys())
x = np.arange(len(species_list))
width = 0.35

# Create bars
exp_bars = ax4.bar(x - width/2, [species_data[s]['expansions'] for s in species_list], 
                   width, label='Expansions', color=expansion_color, edgecolor='black', linewidth=1, alpha=0.9)
cont_bars = ax4.bar(x + width/2, [species_data[s]['contractions'] for s in species_list], 
                    width, label='Contractions', color=contraction_color, edgecolor='black', linewidth=1, alpha=0.9)

ax4.set_ylabel('Number of gene families', fontweight='bold', labelpad=20)
ax4.set_title('e', fontweight='bold', loc='left', pad=25)
ax4.set_xticks(x)
ax4.set_xticklabels(species_list, rotation=45, ha='right', fontsize=17)
ax4.legend(frameon=False, loc='upper right', fontsize=17)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.tick_params(axis='y', labelsize=17)

# Set y-limit
max_species = max([max(species_data[s]['expansions'], species_data[s]['contractions']) for s in species_list])
ax4.set_ylim(0, max_species * 1.2)

# Color-code species labels
for i, (label, species) in enumerate(zip(ax4.get_xticklabels(), species_list)):
    if species_data[species]['group'] == 'Domesticated':
        label.set_color(domesticated_color)
        label.set_fontweight('bold')
    elif species_data[species]['group'] == 'Wild':
        label.set_color(wild_color)
        label.set_fontweight('bold')

# Add value labels with increased font size
for i, bar in enumerate(exp_bars):
    height = bar.get_height()
    if height > 5000:
        ax4.text(bar.get_x() + bar.get_width()/2., height + max_species*0.03,
                 f'{int(height/1000)}k', ha='center', va='bottom', fontsize=15, fontweight='bold')

for i, bar in enumerate(cont_bars):
    height = bar.get_height()
    if height > 5000:
        ax4.text(bar.get_x() + bar.get_width()/2., height + max_species*0.03,
                 f'{int(height/1000)}k', ha='center', va='bottom', fontsize=15, fontweight='bold')


# =============================================================================
# Main title
# =============================================================================

fig.text(0.5, 0.98, 'Gene family evolution during domestication',
         fontsize=24, ha='center', style='italic', color='#2c3e50')  # Increased title size

# =============================================================================
# Adjust layout
# =============================================================================

plt.tight_layout()
plt.subplots_adjust(top=0.95, left=0.12, right=0.95, bottom=0.05)

# =============================================================================
# Save figure
# =============================================================================

output_png = os.path.join(output_dir, 'Figure3_Redesigned_Layout_ExtraLargeFonts.png')
plt.savefig(output_png, dpi=600, bbox_inches='tight', facecolor='white')
print(f"✅ Saved: {output_png}")

output_pdf = os.path.join(output_dir, 'Figure3_Redesigned_Layout_ExtraLargeFonts.pdf')
plt.savefig(output_pdf, format='pdf', dpi=600, bbox_inches='tight', facecolor='white')
print(f"✅ Saved: {output_pdf}")

output_tiff = os.path.join(output_dir, 'Figure3_Redesigned_Layout_ExtraLargeFonts.tiff')
plt.savefig(output_tiff, format='tiff', dpi=600, bbox_inches='tight', facecolor='white')
print(f"✅ Saved: {output_tiff}")

print("\n" + "="*60)
print("FIGURE 3 CREATED SUCCESSFULLY - REDESIGNED LAYOUT WITH EXTRA LARGE FONTS")
print("="*60)
print(f"\n📁 Files saved to: {output_dir}")
print("\n✅ Font sizes increased by 8 points total:")
print("   • Base font: 18 (was 10)")
print("   • Axis labels: 19 (was 11)")
print("   • Axis titles: 20 (was 12)")
print("   • Legend: 17 (was 9)")
print("   • Tick labels: 17 (was 9)")
print("   • Bar value labels: 15-17")
print("   • Figure title: 24 (was 14)")
print("\n✅ New layout:")
print("  Row 1: Panel a (Domesticated vs Wild) and Panel b (Net change)")
print("  Row 2: Panel c (Ratio comparison) and Panel d (Statistical significance)")
print("  Row 3: Panel e (Species breakdown) - Full width")

plt.show()