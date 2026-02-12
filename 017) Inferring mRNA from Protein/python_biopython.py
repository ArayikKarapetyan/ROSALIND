from Bio.Data import CodonTable

def count_rna_strings_biopython(protein):
    """Count using BioPython's codon tables."""
    # Get the standard genetic code table
    standard_table = CodonTable.unambiguous_rna_by_name["Standard"]
    
    # Count codons for each amino acid
    codon_counts = {}
    for codon, aa in standard_table.forward_table.items():
        codon_counts[aa] = codon_counts.get(aa, 0) + 1
    
    # Add stop codons count (not in forward_table)
    stop_codon_count = len(standard_table.stop_codons)
    
    mod = 1000000
    total = 1
    
    for aa in protein:
        if aa not in codon_counts:
            return 0
        total = (total * codon_counts[aa]) % mod
    
    # Multiply by stop codon possibilities
    total = (total * stop_codon_count) % mod
    
    return total


with open("Dataset.txt", "r") as file:
    protein = file.read().strip()
    
print(f"Python BioPython: {count_rna_strings_biopython(protein)}")