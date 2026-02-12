from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO


def find_reverse_palindromes_biopython(fasta_string, min_len=4, max_len=12):
    """Efficient BioPython implementation."""
    fasta_io = StringIO(fasta_string)
    record = SeqIO.read(fasta_io, "fasta")
    dna_seq = str(record.seq)
    
    results = []
    n = len(dna_seq)
    
    # Precompute reverse complement of entire sequence for reference
    seq_obj = Seq(dna_seq)
    rev_comp = str(seq_obj.reverse_complement())
    
    for i in range(n):
        for length in range(min_len, max_len + 1):
            if i + length > n:
                break
            
            # Check palindrome by comparing with reverse complement region
            substring = dna_seq[i:i + length]
            # The reverse complement of this substring is at position n-(i+length) in rev_comp
            rev_comp_sub = rev_comp[n-(i+length):n-i]
            
            if substring == rev_comp_sub:
                results.append((i + 1, length))
    
    return results






with open("Dataset.txt", "r") as file:
    fasta_data = file.read()

palindromes = find_reverse_palindromes_biopython(fasta_data)

for pos, length in palindromes:
    print(f"{pos} {length}")
