def parse_fasta_manual(fasta_string):
    """Parse FASTA format manually."""
    lines = fasta_string.strip().split('\n')
    sequences = {}
    current_id = ""
    
    for line in lines:
        if line.startswith('>'):
            current_id = line[1:].strip()
            sequences[current_id] = ""
        else:
            sequences[current_id] += line.strip()
    
    return sequences

def reverse_complement_manual(dna):
    """Generate reverse complement of DNA string."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna))

def translate_codon(codon):
    """Translate a DNA codon to amino acid."""
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    return genetic_code.get(codon, '')

def translate_orf(dna_sequence, start_pos):
    """Translate an ORF starting from start_pos."""
    protein = ""
    i = start_pos
    
    while i + 3 <= len(dna_sequence):
        codon = dna_sequence[i:i+3]
        amino_acid = translate_codon(codon)
        
        if amino_acid == '*':  # Stop codon
            return protein
        
        protein += amino_acid
        i += 3
    
    return ""  # No stop codon found

def find_all_orfs_manual(dna):
    """Find all ORFs in a DNA sequence and its reverse complement."""
    proteins = set()
    
    # Check all 6 reading frames: 3 forward, 3 reverse complement
    sequences_to_check = [
        dna,  # Forward strand
        reverse_complement_manual(dna)  # Reverse complement
    ]
    
    for seq in sequences_to_check:
        # Check 3 reading frames for each sequence
        for frame in range(3):
            i = frame
            while i + 3 <= len(seq):
                codon = seq[i:i+3]
                if codon == 'ATG':  # Start codon
                    protein = translate_orf(seq, i)
                    if protein:  # Valid ORF (has stop codon)
                        proteins.add(protein)
                    i += 3  # Move to next codon
                else:
                    i += 1  # Check next position
    
    return proteins

def process_dna_orfs_manual(fasta_string):
    """Process FASTA DNA and find all ORFs."""
    sequences = parse_fasta_manual(fasta_string)
    all_proteins = set()
    
    for dna in sequences.values():
        proteins = find_all_orfs_manual(dna)
        all_proteins.update(proteins)
    
    return all_proteins



    
with open("Dataset.txt", "r") as file:
    fasta_data = file.read()

proteins = process_dna_orfs_manual(fasta_data)

for protein in proteins:
    print(protein)
