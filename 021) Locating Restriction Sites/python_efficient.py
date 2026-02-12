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

# More efficient version
def find_reverse_palindromes_efficient(dna, min_len=4, max_len=12):
    """More efficient palindrome finding."""
    results = []
    n = len(dna)
    
    # Precompute complement for faster lookup
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    for i in range(n):
        # For each position, check possible palindrome lengths
        for length in range(min_len, max_len + 1):
            if i + length > n:
                break
            
            # Check palindrome property
            is_palindrome = True
            half_length = length // 2
            
            for k in range(half_length):
                left = dna[i + k]
                right = dna[i + length - 1 - k]
                
                if complement[left] != right:
                    is_palindrome = False
                    break
            
            if is_palindrome:
                results.append((i + 1, length))
    
    return results






with open("Dataset.txt", "r") as file:
    fasta_data = file.read()
sequences = parse_fasta_manual(fasta_data)

for dna in sequences.values():
    palindromes = find_reverse_palindromes_efficient(dna)
    for pos, length in palindromes:
        print(f"{pos} {length}")