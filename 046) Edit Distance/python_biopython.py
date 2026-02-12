from Bio import SeqIO

def edit_distance_biopython(seq1, seq2):
    """
    Edit distance calculation using BioPython sequences.
    """
    s = str(seq1)
    t = str(seq2)
    
    m, n = len(s), len(t)
    
    # Space-optimized DP
    if m > n:
        s, t = t, s
        m, n = n, m
    
    prev = list(range(m + 1))
    
    for j in range(1, n + 1):
        curr = [j]
        for i in range(1, m + 1):
            if s[i-1] == t[j-1]:
                sub_cost = prev[i-1]
            else:
                sub_cost = prev[i-1] + 1
            
            del_cost = prev[i] + 1
            ins_cost = curr[i-1] + 1
            
            curr.append(min(sub_cost, del_cost, ins_cost))
        
        prev = curr
    
    return prev[m]

def main():
    # Read sequences using BioPython
    records = list(SeqIO.parse("Dataset.txt", "fasta"))
    seq1 = records[0].seq
    seq2 = records[1].seq
    
    # Compute edit distance
    distance = edit_distance_biopython(seq1, seq2)
    
    # Output
    print(distance)

if __name__ == "__main__":
    main()