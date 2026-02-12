from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO

def find_orfs_biopython(fasta_string):
    """Find ORFs using BioPython."""
    proteins = set()
    
    # Parse FASTA
    fasta_io = StringIO(fasta_string)
    record = SeqIO.read(fasta_io, "fasta")
    dna_seq = record.seq
    
    # Check both strands
    for strand, nuc in [(+1, dna_seq), (-1, dna_seq.reverse_complement())]:
        for frame in range(3):
            # Translate in this frame
            length = 3 * ((len(nuc) - frame) // 3)
            
            # Translate this frame
            trans = str(nuc[frame:frame+length].translate(to_stop=True))
            
            # Find all ATG positions in this frame
            trans_len = len(trans)
            seq = str(nuc[frame:frame+length])
            
            # Look for start codons (ATG) and translate from there
            for pos in range(0, len(seq) - 2, 3):
                if seq[pos:pos+3] == 'ATG':
                    # Translate from this ATG
                    subseq = seq[pos:]
                    prot = str(Seq(subseq).translate(to_stop=True))
                    
                    # Only add if we actually hit a stop codon
                    if len(prot) * 3 <= len(subseq):
                        proteins.add(prot)
    
    return proteins

# Alternative simpler approach
def find_orfs_biopython_simple(dna_seq):
    """Simpler BioPython ORF finder."""
    proteins = set()
    
    # Convert to Seq object if needed
    if not isinstance(dna_seq, Seq):
        dna_seq = Seq(dna_seq)
    
    # Check all 6 reading frames
    for strand in [dna_seq, dna_seq.reverse_complement()]:
        for frame in range(3):
            # Get translation for this frame
            trans = str(strand[frame:].translate(table=1, to_stop=False))
            
            # Split on stop codons (*)
            parts = trans.split('*')
            
            # For each part, find proteins starting with M
            for part in parts:
                # Find all M positions
                for i in range(len(part)):
                    if part[i] == 'M':
                        protein = part[i:]
                        # Check if this is a valid ORF (not interrupted by earlier stop)
                        # Since we split on stops, all proteins here start with M
                        # and go to the end of the part (which would be next stop)
                        proteins.add(protein)
    
    return proteins


with open("Dataset.txt", "r") as file:
    fasta_data = file.read()

proteins = find_orfs_biopython(fasta_data)

for protein in proteins:
    print(protein)