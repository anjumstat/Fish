# Method: Using ETE Toolkit - Enhanced Fish Phylogeny with Clean Styling
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import re
import os

# Path to your species tree
tree_file = "E:/fish/orth/OrthoFinder/Results_Feb11/Species_Tree/tree/SpeciesTree_rooted.txt"

# Output directory
output_dir = "E:/fish/paper/figures/Figure1/"
os.makedirs(output_dir, exist_ok=True)

# Read tree
t = Tree(tree_file)

# Print original names for debugging
print("Original leaf names:")
for leaf in t.iter_leaves():
    print(f"  - {leaf.name}")

# Create a dictionary for clean species names with common names
clean_names = {
    "cleaned_Oryzias_sinensis.ASM858656v1.pep.all": "Oryzias sinensis (Japanese rice fish)",
    "cleaned_Larimichthys_crocea.L_crocea_2.0.pep.all": "Larimichthys crocea (Large yellow croaker)",
    "cleaned_Sinocyclocheilus_grahami.SAMN03320097.WGS_v1.1.pep.all": "Sinocyclocheilus grahami (Golden-line fish)",
    "cleaned_Sinocyclocheilus_anshuiensis.SAMN03320099.WGS_v1.1.pep.all": "Sinocyclocheilus anshuiensis (Anshui cavefish)",
    "cleaned_Carassius_auratus.ASM336829v1.pep.all": "Carassius auratus (Goldfish)",
    "cleaned_Cyprinus_carpio_carpio.Cypcar_WagV4.0.pep.all": "Cyprinus carpio (Common carp)"  # Fixed: removed _carpio_carpio issue
}

# Define species groups with multiple possible name formats for matching
domesticated_species_patterns = [
    "Carassius auratus (Goldfish)",
    "Cyprinus carpio (Common carp)",
    "Carassius auratus",
    "Cyprinus carpio",
    "goldfish",
    "carp"
]

sinocyclocheilus_species_patterns = [
    "Sinocyclocheilus grahami (Golden-line fish)",
    "Sinocyclocheilus anshuiensis (Anshui cavefish)",
    "Sinocyclocheilus grahami",
    "Sinocyclocheilus anshuiensis",
    "grahami",
    "anshuiensis"
]

outgroup_species_patterns = [
    "Oryzias sinensis (Japanese rice fish)",
    "Larimichthys crocea (Large yellow croaker)",
    "Oryzias sinensis",
    "Larimichthys crocea",
    "sinensis",
    "crocea"
]

# Store original names for node finding
original_to_clean = {}
cleaned_names_list = []  # Store cleaned names for debugging

for leaf in t.iter_leaves():
    original_name = leaf.name
    if original_name in clean_names:
        clean_name = clean_names[original_name]
        original_to_clean[original_name] = clean_name
        leaf.name = clean_name  # Set clean name for display
        cleaned_names_list.append(clean_name)
    else:
        # Fallback cleaning if not in dictionary
        clean = re.sub(r'cleaned_', '', original_name)
        clean = re.sub(r'\..*', '', clean)
        leaf.name = clean
        cleaned_names_list.append(clean)

print("\nCleaned species names:")
for name in cleaned_names_list:
    print(f"  - {name}")

# Find nodes for each species with flexible matching
def find_node_by_name(tree, name_pattern):
    """Find node by partial name matching"""
    for leaf in tree.iter_leaves():
        if name_pattern.lower() in leaf.name.lower():
            print(f"  ✓ Matched '{name_pattern}' to '{leaf.name}'")
            return leaf
    return None

# Get nodes for each species group
domesticated_nodes = []
for pattern in domesticated_species_patterns:
    node = find_node_by_name(t, pattern)
    if node and node not in domesticated_nodes:
        domesticated_nodes.append(node)
        print(f"  Added domesticated node: {node.name}")

sinocyclocheilus_nodes = []
for pattern in sinocyclocheilus_species_patterns:
    node = find_node_by_name(t, pattern)
    if node and node not in sinocyclocheilus_nodes:
        sinocyclocheilus_nodes.append(node)
        print(f"  Added Sinocyclocheilus node: {node.name}")

outgroup_nodes = []
for pattern in outgroup_species_patterns:
    node = find_node_by_name(t, pattern)
    if node and node not in outgroup_nodes:
        outgroup_nodes.append(node)
        print(f"  Added outgroup node: {node.name}")

print(f"\nFound {len(domesticated_nodes)} domesticated nodes:")
for node in domesticated_nodes:
    print(f"  - {node.name}")

print(f"Found {len(sinocyclocheilus_nodes)} Sinocyclocheilus nodes:")
for node in sinocyclocheilus_nodes:
    print(f"  - {node.name}")

print(f"Found {len(outgroup_nodes)} outgroup nodes:")
for node in outgroup_nodes:
    print(f"  - {node.name}")

# Find MRCA for each group
def find_mrca(tree, nodes):
    if len(nodes) >= 2:
        return tree.get_common_ancestor(nodes)
    elif len(nodes) == 1:
        return nodes[0]
    else:
        return None

domesticated_node = find_mrca(t, domesticated_nodes)
sinocyclocheilus_node = find_mrca(t, sinocyclocheilus_nodes)
outgroup_node = find_mrca(t, outgroup_nodes)

print(f"\nMRCA nodes found:")
print(f"  Domesticated MRCA: {domesticated_node is not None}")
print(f"  Sinocyclocheilus MRCA: {sinocyclocheilus_node is not None}")
print(f"  Outgroup MRCA: {outgroup_node is not None}")

# Create TreeStyle with adjusted margins to bring scale bar closer
ts = TreeStyle()
ts.show_leaf_name = False
ts.show_branch_length = True
ts.show_scale = True
ts.scale_length = 0.09  # Adjust scale bar size
ts.branch_vertical_margin = 25
ts.margin_top = 50  # Reduced from 70
ts.margin_bottom = 100  # Reduced from 200
ts.margin_left = 50  # Reduced from 70
ts.margin_right = 450  # Keep space for labels
ts.mode = "r"

# Scale bar settings (removed invalid scale_font_size)
# ETE3 automatically handles scale bar font

# Add title with smaller margin
title = TextFace("Phylogenetic Tree of Cyprinid Fishes: Domesticated vs Wild Cave Species", 
                 fsize=16, fgcolor="black", fstyle="bold")  # Reduced font size
title.margin_bottom = 10  # Reduced from 20
ts.title.add_face(title, column=0)

# Color the clades
def color_clade(node, color_name, light_color, line_width=3):
    if node:
        print(f"  Coloring clade with {color_name}")
        for n in node.traverse():
            style = NodeStyle()
            style["bgcolor"] = light_color
            style["vt_line_color"] = color_name
            style["hz_line_color"] = color_name
            style["vt_line_width"] = line_width
            style["hz_line_width"] = line_width
            style["size"] = 0
            n.set_style(style)

# Color domesticated clade (red/pink)
if domesticated_node:
    print("Applying domesticated clade color...")
    color_clade(domesticated_node, "darkred", "#FFB6C1", line_width=3)
else:
    print("WARNING: No domesticated node found - coloring individual species instead")
    for node in domesticated_nodes:
        if node:
            color_clade(node, "darkred", "#FFB6C1", line_width=2)

# Color Sinocyclocheilus clade (blue)
if sinocyclocheilus_node:
    print("Applying Sinocyclocheilus clade color...")
    color_clade(sinocyclocheilus_node, "darkblue", "#B0E0E6", line_width=3)
else:
    print("WARNING: No Sinocyclocheilus node found - coloring individual species instead")
    for node in sinocyclocheilus_nodes:
        if node:
            color_clade(node, "darkblue", "#B0E0E6", line_width=2)

# Color outgroup species individually
for node in outgroup_nodes:
    if node:
        if "Oryzias" in node.name or "sinensis" in node.name.lower():
            color_clade(node, "darkgreen", "#C1E1C1", line_width=2)
        elif "Larimichthys" in node.name or "crocea" in node.name.lower():
            color_clade(node, "darkorange", "#FFDAB9", line_width=2)

# Layout function for labels
def my_layout(node):
    if node.is_leaf():
        # Determine color based on species group
        text_color = "black"
        font_style = "normal"
        
        # Check against patterns
        name_lower = node.name.lower()
        
        if any(pattern.lower() in name_lower for pattern in ["carassius", "goldfish", "cyprinus", "carp"]):
            text_color = "darkred"
            font_style = "bold"
        elif any(pattern.lower() in name_lower for pattern in ["sinocyclocheilus", "grahami", "anshuiensis"]):
            text_color = "darkblue"
        elif "oryzias" in name_lower or "sinensis" in name_lower:
            text_color = "darkgreen"
        elif "larimichthys" in name_lower or "crocea" in name_lower:
            text_color = "darkorange"
        
        # Add species name with common name
        name_face = TextFace(f"  {node.name}", fsize=11,  # Slightly smaller font
                           fgcolor=text_color, fstyle=font_style)
        name_face.margin_left = 3  # Reduced margin
        node.add_face(name_face, column=0, position="branch-right")
        
    elif node.name and node.name.replace('.', '').replace('-', '').isdigit():
        try:
            if float(node.name) >= 70:
                boot_face = TextFace(f"{node.name}", fsize=9, fgcolor="gray")  # Removed space
                node.add_face(boot_face, column=0, position="branch-top")
        except:
            pass

ts.layout_fn = my_layout

# Create legend with reduced spacing
legend_y_position = 15  # Reduced from 30

# Domesticated species legend
legend1 = TextFace("  ● Domesticated: Goldfish & Common Carp", 
                   fsize=11, fgcolor="darkred", fstyle="bold")  # Smaller font
legend1.margin_top = legend_y_position
legend1.margin_left = 40  # Reduced from 50
legend1.margin_bottom = 3  # Reduced

# Wild cavefish legend
legend2 = TextFace("  ● Wild Cavefish: S. grahami & S. anshuiensis", 
                   fsize=11, fgcolor="darkblue")
legend2.margin_top = 3
legend2.margin_left = 40
legend2.margin_bottom = 3

# Outgroup legend
legend3 = TextFace("  ● Outgroup 1: Oryzias sinensis", 
                   fsize=11, fgcolor="darkgreen")
legend3.margin_top = 3
legend3.margin_left = 40
legend3.margin_bottom = 3

legend4 = TextFace("  ● Outgroup 2: Larimichthys crocea", 
                   fsize=11, fgcolor="darkorange")
legend4.margin_top = 3
legend4.margin_left = 40
legend4.margin_bottom = 3

# Add CAFE results summary note with reduced spacing
note1 = TextFace("  * Domesticated: 35,786 expansions (+31,402)", 
                fsize=11, fgcolor="darkred", fstyle="bold")
note1.margin_top = 15  # Reduced from 20
note1.margin_left = 40
note1.margin_bottom = 2

note2 = TextFace("  * Wild cavefish: 30,213 contractions (-17,483)", 
                fsize=11, fgcolor="darkblue", fstyle="bold")
note2.margin_top = 2
note2.margin_left = 40
note2.margin_bottom = 2

note3 = TextFace("  * All 24 COG categories represented (fold ~1.0)", 
                fsize=11, fgcolor="black")
note3.margin_top = 2
note3.margin_left = 40

# Add all legends to title
ts.title.add_face(legend1, column=0)
ts.title.add_face(legend2, column=0)
ts.title.add_face(legend3, column=0)
ts.title.add_face(legend4, column=0)
ts.title.add_face(note1, column=0)
ts.title.add_face(note2, column=0)
ts.title.add_face(note3, column=0)

# Save files with slightly smaller dimensions to bring elements closer
png_path = os.path.join(output_dir, "fish_phylogeny_colored_final.png")
pdf_path = os.path.join(output_dir, "fish_phylogeny_colored_final.pdf")

t.render(png_path, w=1800, h=1000, tree_style=ts)  # Reduced from 2000x1200
t.render(pdf_path, tree_style=ts)

print(f"\n✅ Tree saved to: {png_path}")
print(f"✅ PDF saved to: {pdf_path}")

# Open the image (Windows only)
if os.name == 'nt':
    os.startfile(png_path)

print("\n✅ Enhanced fish phylogeny complete with scale bar positioned closer!")