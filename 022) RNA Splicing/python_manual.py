def parse_fasta_manual(fasta_string):
    """Parse FASTA format manually."""
    sequences = {}
    current_id = ""
    
    for line in fasta_string.strip().split('\n'):
        if line.startswith('>'):
            current_id = line[1:].strip()
            sequences[current_id] = ""
        else:
            sequences[current_id] += line.strip()
    
    return sequences

def remove_introns_manual(dna, introns):
    """Remove introns from DNA sequence."""
    # Start with the full DNA
    result = dna
    
    # Remove each intron
    for intron in introns:
        result = result.replace(intron, "")
    
    return result

def dna_to_rna_manual(dna):
    """Transcribe DNA to RNA by replacing T with U."""
    return dna.replace('T', 'U')

def translate_rna_to_protein_manual(rna):
    """Translate RNA to protein using genetic code."""
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    protein = []
    
    # Process in steps of 3 (codons)
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        amino_acid = genetic_code.get(codon, '')
        
        # Stop codon - end translation
        if amino_acid == '*':
            break
        
        if amino_acid:
            protein.append(amino_acid)
    
    return ''.join(protein)

def rna_splicing_manual(fasta_string):
    """Main function for RNA splicing."""
    sequences = parse_fasta_manual(fasta_string)
    
    # First sequence is the DNA, rest are introns
    seq_ids = list(sequences.keys())
    dna = sequences[seq_ids[0]]
    introns = [sequences[seq_id] for seq_id in seq_ids[1:]]
    
    # Step 1: Remove introns
    exons_only = remove_introns_manual(dna, introns)
    
    # Step 2: Transcribe to RNA
    rna = dna_to_rna_manual(exons_only)
    
    # Step 3: Translate to protein
    protein = translate_rna_to_protein_manual(rna)
    
    return protein


with open("Dataset.txt", "r") as file:
    fasta_data = file.read()
protein = rna_splicing_manual(fasta_data)

print("Python Manual Solution:")
print(protein)


