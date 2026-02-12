from Bio import Phylo
from io import StringIO

def get_character_table_biopython(newick_str):
    """Get character table using BioPython's Phylo module."""
    # Parse tree
    tree = Phylo.read(StringIO(newick_str), "newick")
    
    # Get all terminal (leaf) names
    all_leaves = [clade.name for clade in tree.get_terminals()]
    all_leaves.sort()  # Lexicographic order
    leaf_to_idx = {leaf: i for i, leaf in enumerate(all_leaves)}
    
    table = []
    
    # Function to process clades
    def process_clade(clade, parent_leaves=None):
        """Process clade and add splits."""
        if parent_leaves is None:
            parent_leaves = set(all_leaves)
        
        # Get leaves in this clade
        clade_leaves = {leaf.name for leaf in clade.get_terminals()}
        
        # Check if split is non-trivial
        if len(clade_leaves) > 1 and len(clade_leaves) < len(all_leaves):
            # Create character row
            row = ['0'] * len(all_leaves)
            for leaf in clade_leaves:
                row[leaf_to_idx[leaf]] = '1'
            table.append(''.join(row))
        
        # Process children
        for child in clade.clades:
            if not child.is_terminal():  # Only internal nodes
                process_clade(child, clade_leaves)
    
    # Start processing from root
    process_clade(tree.root)
    
    return table

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        newick = f.readline().strip()
    
    # Get character table
    table = get_character_table_biopython(newick)
    
    # Output
    for row in table:
        print(row)

if __name__ == "__main__":
    main()