def parse_input_manual(input_string):
    """Parse input: first line n, second line permutation."""
    lines = input_string.strip().split('\n')
    n = int(lines[0].strip())
    permutation = list(map(int, lines[1].strip().split()))
    return n, permutation

# More efficient LIS using patience sorting (O(n log n))
def longest_increasing_subsequence_efficient(seq):
    """Find LIS using patience sorting algorithm."""
    if not seq:
        return []
    
    # piles will store the top card of each pile
    piles = []
    # parent pointers for reconstruction
    parent = [-1] * len(seq)
    pile_tops = []  # indices of pile tops
    
    for i, num in enumerate(seq):
        # Binary search for leftmost pile where top >= num
        left, right = 0, len(piles)
        while left < right:
            mid = (left + right) // 2
            if seq[piles[mid]] < num:
                left = mid + 1
            else:
                right = mid
        
        if left < len(piles):
            piles[left] = i
        else:
            piles.append(i)
        
        # Store parent pointer (previous pile's top)
        if left > 0:
            parent[i] = piles[left - 1]
        else:
            parent[i] = -1
        
        pile_tops = piles[:]
    
    # Reconstruct LIS
    lis = []
    curr = piles[-1]
    while curr != -1:
        lis.append(seq[curr])
        curr = parent[curr]
    
    return lis[::-1]

# Efficient LDS by finding LIS on reversed or negated sequence
def longest_decreasing_subsequence_efficient(seq):
    """Find LDS using efficient method."""
    # LDS is LIS on reversed sequence, or LIS on negated values
    # Method 1: Find LIS on reversed sequence
    rev_seq = seq[::-1]
    # But we need decreasing in original order
    # Better: Find LIS on negated values
    neg_seq = [-x for x in seq]
    lis_neg = longest_increasing_subsequence_efficient(neg_seq)
    lds = [-x for x in lis_neg]
    return lds

def process_permutation_efficient(input_string):
    """Main processing function."""
    n, permutation = parse_input_manual(input_string)
    
    lis = longest_increasing_subsequence_efficient(permutation)
    lds = longest_decreasing_subsequence_efficient(permutation)
    
    return lis, lds

with open("Dataset.txt", "r") as file:
    input_data = file.read()

lis, lds = process_permutation_efficient(input_data)

print(' '.join(map(str, lis)))
print(' '.join(map(str, lds)))
