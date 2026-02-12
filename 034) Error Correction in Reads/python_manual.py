def reverse_complement(dna):
    comp = {'A':'T','T':'A','C':'G','G':'C'}
    return ''.join(comp[b] for b in reversed(dna))

def hamming_distance(s1, s2):
    return sum(a!=b for a,b in zip(s1,s2))

def correct_errors(reads):
    from collections import Counter
    counts = Counter(reads)
    
    correct = set([r for r,c in counts.items() if c > 1])
    
    corrections = []
    
    for read, c in counts.items():
        if c > 1:
            continue  # already correct
        # check against correct reads
        for ref in correct:
            if hamming_distance(read, ref) == 1 or hamming_distance(read, reverse_complement(ref)) == 1:
                corrections.append(f"{read}->{ref}")
                break
    
    return corrections

# Example main
def main():
    reads = []
    with open("Dataset.txt","r") as f:
        current = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    reads.append(current)
                    current = ""
            else:
                current += line
        if current:
            reads.append(current)
    
    corrections = correct_errors(reads)
    for corr in corrections:
        print(corr)

if __name__=="__main__":
    main()
