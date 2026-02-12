def edit_distance(s, t):
    """
    Compute edit distance (Levenshtein distance) between two strings.
    Uses dynamic programming with O(mn) time and O(min(m,n)) space.
    """
    m, n = len(s), len(t)
    
    # Ensure s is the shorter string for space optimization
    if m > n:
        s, t = t, s
        m, n = n, m
    
    # Previous row of distances (initially for empty string)
    prev = list(range(m + 1))
    
    # Compute current row distances
    for j in range(1, n + 1):
        # Current row starts with distance for empty s
        curr = [j]
        
        for i in range(1, m + 1):
            # Cost for substitution
            if s[i-1] == t[j-1]:
                sub_cost = prev[i-1]
            else:
                sub_cost = prev[i-1] + 1
            
            # Costs for deletion and insertion
            del_cost = prev[i] + 1      # delete from s
            ins_cost = curr[i-1] + 1    # insert into s
            
            # Take minimum of three operations
            curr.append(min(sub_cost, del_cost, ins_cost))
        
        # Move to next row
        prev = curr
    
    return prev[m]

def main():
    # Read FASTA file
    sequences = []
    with open("Dataset.txt", "r") as f:
        current_seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line
        if current_seq:
            sequences.append(current_seq)
    
    # Get the two protein strings
    s, t = sequences[0], sequences[1]
    
    # Compute edit distance
    distance = edit_distance(s, t)
    
    # Output result
    print(distance)

if __name__ == "__main__":
    main()