# Fish
# Cyprinid Domestication Genomics - Data and Scripts Repository

## 📋 Overview
This repository contains all data and scripts required to reproduce the analyses and figures from our study on gene family evolution during cyprinid domestication. The study compares domesticated species (Carassius auratus - goldfish, Cyprinus carpio - common carp) with wild cavefish relatives (Sinocyclocheilus anshuiensis, S. grahami) and outgroup species (Larimichthys crocea, Oryzias sinensis).

---

## 📁 Repository Structure
├─
│ ├── 📁 orthofinder_results/ # OrthoFinder output (Results_Feb11)
│ │ ├── Orthogroups.tsv # 34,707 orthogroups across 6 species
│ │ ├── Orthogroups.GeneCount.tsv # Gene counts per orthogroup
│ │ └── SpeciesTree_rooted.txt # Rooted species tree for CAFE5
│ ├── 📁 cafe_results/ # CAFE5 analysis results
│ │ ├── Gamma_change.tab # Expansion/contraction values (24,965 families)
│ │ └── expanded_genes_cyprinids.txt # 117,761 expanded gene IDs
│ └── 📁 eggnog_annotations/ # Functional annotations
│ └── fish_domestication.emapper.annotations # eggNOG annotations (100,546 genes)
│
├── 📁 scripts/ # All Python scripts
│ ├── 01_clean_proteomes.py # Protein sequence cleaning
│ ├── 02_figure1_phylogeny.py # Phylogenetic tree (Figure 1)
│ ├── 03_figure2_orthogroup_overview.py # Orthogroup analysis (Figure 2)
│ ├── 04_figure3_CAFE_results.py # Expansion/contraction patterns (Figure 3)
│ ├── 05_figure4_functional_enrichment.py # Functional enrichment (Figure 4)
│ └── 06_figure5_hypothesis_testing.py # Hypothesis testing (Figure 5)

---

## 📊 Data Files Description

### OrthoFinder Results (`/data/orthofinder_results/`)
| File | Description |
|------|-------------|
| `Orthogroups.tsv` | 34,707 orthogroups across all 6 species. Each row represents an orthogroup with gene IDs for each species. |
| `Orthogroups.GeneCount.tsv` | Gene counts per orthogroup for each species. Used as input for CAFE5. |
| `SpeciesTree_rooted.txt` | Rooted phylogenetic tree in Newick format. Used as input for CAFE5. |

### CAFE5 Results (`/data/cafe_results/`)
| File | Description |
|------|-------------|
| `Gamma_change.tab` | CAFE5 output showing expansion/contraction values for 24,965 gene families across all species. Positive values = expansions, negative values = contractions. |
| `expanded_genes_cyprinids.txt` | List of 117,761 expanded gene IDs from domesticated lineages (goldfish and common carp). One gene ID per line. |

### eggNOG Annotations (`/data/eggnog_annotations/`)
| File | Description |
|------|-------------|
| `fish_domestication.emapper.annotations` | eggNOG-mapper functional annotations for 100,546 genes. Includes COG categories, GO terms, KEGG pathways, and descriptions. |

---

## 🐍 Scripts Description

### Script 1: `01_clean_proteomes.py`
**Purpose:** Cleans raw protein FASTA files before OrthoFinder analysis.
**Functions:**
- Removes invalid amino acid characters
- Filters sequences <20 amino acids
- Filters sequences with >50% gaps
- Cleans sequence headers
- Generates cleaning summary CSV

**Input:** Raw protein FASTA files (.faa, .fasta, .fa, .pep)
**Output:** Cleaned FASTA files (prefixed with "cleaned_") and cleaning_summary.csv

---

### Script 2: `02_figure1_phylogeny.py`
**Purpose:** Generates Figure 1 - Colored phylogenetic tree.
**Input:** `data/orthofinder_results/Species_Tree/tree/SpeciesTree_rooted.txt`
**Features:**
- Colors domesticated species (goldfish, common carp) in darkred with pink background
- Colors wild cavefish species in darkblue with light blue background
- Colors outgroup species in darkgreen (Oryzias) and darkorange (Larimichthys)
- Adds legends with CAFE results summary
- Includes scale bar and bootstrap values

---

### Script 3: `03_figure2_orthogroup_overview.py`
**Purpose:** Generates Figure 2 - Orthogroup overview across species.
**Input:** 
- `data/orthofinder_results/Orthogroups/Orthogroups.tsv`
- `data/orthofinder_results/Orthogroups/Orthogroups.GeneCount.tsv`

**Panels:**
- A: Presence/absence heatmap (first 100 orthogroups)
- B: Orthogroup size distribution (log scale)
- C: Single-copy vs multi-copy orthogroups
- D: Gene counts per species
- E: Summary statistics table

**Key Statistics Generated:**
- Total orthogroups: 34,707
- Total genes: 575,162
- Single-copy orthogroups (all species): 79
- Multi-copy orthogroups (all species): 14,937
- Domesticated-specific orthogroups: 4,076
- Wild-specific orthogroups: 4,482
- Shared orthogroups: 24,140

---

### Script 4: `04_figure3_CAFE_results.py`
**Purpose:** Generates Figure 3 - Expansion/contraction patterns from CAFE5.

**Data (hard-coded from CAFE results):**
| Metric | Value |
|--------|-------|
| Domesticated expansions | 35,786 |
| Domesticated contractions | 4,384 |
| Wild expansions | 12,730 |
| Wild contractions | 30,213 |
| Total families | 24,965 |
| Significant (p<0.05) | 6,388 (25.6%) |

**Species-specific Data:**
| Species | Expansions | Contractions | Group |
|---------|------------|--------------|-------|
| C. auratus | 18,215 | 2,163 | Domesticated |
| C. carpio | 17,571 | 2,221 | Domesticated |
| S. anshuiensis | 5,516 | 16,097 | Wild |
| S. grahami | 7,214 | 14,116 | Wild |
| L. crocea | 14,770 | 5,863 | Outgroup |
| O. sinensis | 13,883 | 6,001 | Outgroup |

**Panels:**
- a: Domesticated vs Wild expansions/contractions
- b: Net change (+31,402 domesticated; -17,483 wild)
- c: Ratio comparison (2.8× more expansions, 6.9× fewer contractions)
- d: Statistical significance
- e: Species-specific breakdown

---

### Script 5: `05_figure4_functional_enrichment.py`
**Purpose:** Generates Figure 4 - Functional enrichment analysis of expanded genes.
**Input:** 
- `data/eggnog_annotations/fish_domestication.emapper.annotations`
- `data/cafe_results/expanded_genes_cyprinids.txt`

**Key Findings:**
| Metric | Value |
|--------|-------|
| Expanded genes with annotations | 29,422 |
| Total background genes | 100,546 |
| Significant COG categories | Only COG L (Replication) |
| Mean fold enrichment | 0.925 |
| Total GO terms | 14,830 |
| Total KEGG pathways | 720 |

**Hypothesis Categories:**
| Category | Expanded % | Background % | Enrichment |
|----------|------------|--------------|------------|
| Growth | 5.20% | 5.53% | 0.94× |
| Metabolism | 13.60% | 14.03% | 0.97× |
| Stress | 39.28% | 39.49% | 0.99× |

**Top GO Terms:**
| GO Term | Description | Expanded Count |
|---------|-------------|----------------|
| GO:0008150 | biological_process | 21,174 |
| GO:0005575 | cellular_component | 20,926 |
| GO:0005623 | cell | 20,381 |

**Top KEGG Pathways:**
| KEGG ID | Pathway | Expanded Count |
|---------|---------|----------------|
| ko01100 | Metabolic pathways | 2,120 |
| ko05200 | Pathways in cancer | 1,280 |
| ko04010 | MAPK signaling | 869 |

---

### Script 6: `06_figure5_hypothesis_testing.py`
**Purpose:** Generates Figure 5 - Analysis of growth, metabolism, and stress categories.
**Input:** Same as Script 5

**Key Findings:**
| Category | Count | Percentage |
|----------|-------|------------|
| Growth | 1,531 | 5.20% |
| Metabolism | 4,001 | 13.60% |
| Stress | 11,557 | 39.28% |
| Any hypothesis | 16,363 | 55.61% |
| Other | 13,059 | 44.39% |

**Overlap Analysis:**
| Category | Count |
|----------|-------|
| Growth only | 1,242 |
| Metabolism only | 3,555 |
| Stress only | 10,840 |
| Growth + Metabolism | 9 |
| Growth + Stress | 280 |
| Metabolism + Stress | 437 |
| All three | 0 |

**Background Comparison:**
| Category | Background % | Enrichment |
|----------|--------------|------------|
| Growth | 5.53% | 0.94× |
| Metabolism | 14.03% | 0.97× |
| Stress | 39.49% | 0.99× |
| Any hypothesis | 56.35% | 0.99× |

---

## 🐍 Software Requirements

### Python Packages Required
pandas
numpy
matplotlib
seaborn
scipy
statsmodels
biopython
ete3

### Install with pip:
```bash
pip install pandas numpy matplotlib seaborn scipy statsmodels biopython ete3
External Tools (not included in scripts)
Tool	Purpose
OrthoFinder	Orthogroup inference
CAFE5	Gene family evolution analysis
eggNOG-mapper	Functional annotation
🔄 Workflow Order
Run 01_clean_proteomes.py on raw protein sequences

Run OrthoFinder on cleaned proteomes

Run CAFE5 using OrthoFinder outputs

Run eggNOG-mapper on proteomes

Run figure scripts in order (02-06)

📝 Figure Outputs
All figure scripts generate high-resolution images (600 DPI) in PNG and PDF formats suitable for publication.

Script	Output Files
02_figure1_phylogeny.py	Figure1_phylogeny.png, Figure1_phylogeny.pdf
03_figure2_orthogroup_overview.py	figure2_orthogroup_overview.png, figure2_orthogroup_overview.pdf
04_figure3_CAFE_results.py	Figure3_CAFE_results.png, Figure3_CAFE_results.pdf
05_figure4_functional_enrichment.py	Figure4_functional_enrichment.png, Figure4_functional_enrichment.pdf
06_figure5_hypothesis_testing.py	Figure5_hypothesis_testing.png, Figure5_hypothesis_testing.pdf
Contact
For questions about these data or scripts, please contact the corresponding author.
