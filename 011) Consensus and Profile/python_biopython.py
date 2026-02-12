from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from io import StringIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

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

def compute_profile_and_consensus_biopython(fasta_string):
    """Compute profile and consensus using BioPython."""
    # Parse sequences
    fasta_io = StringIO(fasta_string)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    
    # Extract sequences
    seq_list = [str(record.seq) for record in records]
    n = len(seq_list[0])
    
    # Initialize profile
    profile = {
        'A': [0] * n,
        'C': [0] * n,
        'G': [0] * n,
        'T': [0] * n
    }
    
    # Fill profile
    for seq in seq_list:
        for j, nucleotide in enumerate(seq):
            profile[nucleotide][j] += 1
    
    # Compute consensus
    consensus_chars = []
    for j in range(n):
        counts = [profile[nuc][j] for nuc in ['A', 'C', 'G', 'T']]
        max_idx = counts.index(max(counts))
        consensus_chars.append(['A', 'C', 'G', 'T'][max_idx])
    
    consensus = ''.join(consensus_chars)
    
    return consensus, profile

# Alternative using BioPython's SummaryInfo
def compute_profile_and_consensus_summary(fasta_string):
    """Compute using BioPython's SummaryInfo."""

    
    fasta_io = StringIO(fasta_string)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    
    # Create alignment
    alignment = MultipleSeqAlignment(records)
    
    # Get summary information
    summary = alignment.__dict__.get('_per_column_annotations', {})
    
    # Manual calculation if summary not available
    n = alignment.get_alignment_length()
    m = len(alignment)
    
    profile = {'A': [0]*n, 'C': [0]*n, 'G': [0]*n, 'T': [0]*n}
    
    for i in range(m):
        seq = str(alignment[i].seq)
        for j in range(n):
            profile[seq[j]][j] += 1
    
    # Consensus
    consensus_chars = []
    for j in range(n):
        max_count = -1
        max_nuc = 'A'
        for nuc in ['A', 'C', 'G', 'T']:
            if profile[nuc][j] > max_count:
                max_count = profile[nuc][j]
                max_nuc = nuc
        consensus_chars.append(max_nuc)
    
    consensus = ''.join(consensus_chars)
    
    return consensus, profile



with open("Dataset.txt", "r") as file:
    fasta_data = file.read()

sequences = parse_fasta_manual(fasta_data)
consensus, profile = compute_profile_and_consensus_biopython(sequences)

print("Consensus:", consensus)
print("\nProfile Matrix:")
for nuc in ['A', 'C', 'G', 'T']:
    print(f"{nuc}: {' '.join(map(str, profile[nuc]))}")