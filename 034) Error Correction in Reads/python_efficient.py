from collections import Counter, defaultdict

def reverse_complement(dna):
    """Fast reverse complement using translation table."""
    trans = str.maketrans('ACGT', 'TGCA')
    return dna.translate(trans)[::-1]

def correct_errors_efficient(reads):
    """Efficient error correction using frequency counting and neighbor checking."""
    # Count frequencies
    freq = Counter()
    for read in reads:
        freq[read] += 1
        freq[reverse_complement(read)] += 1
    
    # Identify correct reads (total count >= 2)
    correct = {read for read, count in freq.items() if count >= 2}
    
    # Create a set for O(1) lookups
    correct_set = set(correct)
    
    # Precompute all possible single mutations for correct reads
    mutation_map = defaultdict(list)
    for correct_read in correct:
        seq = correct_read
        # Generate all single mutations
        for i in range(len(seq)):
            for base in 'ACGT':
                if seq[i] != base:
                    mutation = seq[:i] + base + seq[i+1:]
                    mutation_map[mutation].append(correct_read)
    
    # Find corrections
    corrections = []
    seen = set()
    
    for read in reads:
        if read in correct_set:
            continue
        
        # Check if this read has exactly one correction
        if read in mutation_map:
            possible_corrections = mutation_map[read]
            
            # Filter only those where reverse complement is also considered
            valid_corrections = []
            for corr in possible_corrections:
                if corr in correct_set or reverse_complement(corr) in correct_set:
                    valid_corrections.append(corr)
            
            # Should have exactly one valid correction
            if len(valid_corrections) == 1:
                correction = f"{read}->{valid_corrections[0]}"
                if correction not in seen:
                    corrections.append(correction)
                    seen.add(correction)
    
    return corrections

def main():
    # Read FASTA file
    reads = []
    with open("Dataset.txt", "r") as f:
        current_read = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_read:
                    reads.append(current_read)
                current_read = ""
            else:
                current_read = line  # Assuming single line per sequence
        if current_read:
            reads.append(current_read)
    
    # Correct errors
    corrections = correct_errors_efficient(reads)
    
    # Output corrections
    for correction in corrections:
        print(correction)

if __name__ == "__main__":
    main()