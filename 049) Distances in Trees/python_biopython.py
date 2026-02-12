from Bio import Phylo
from io import StringIO

def distance_in_tree_biopython(newick_str, node1_name, node2_name):
    tree = Phylo.read(StringIO(newick_str), "newick")

    # IMPORTANT: set all branch lengths to 1
    for clade in tree.find_clades():
        if clade.branch_length is None:
            clade.branch_length = 1

    node1 = next(tree.find_clades(name=node1_name), None)
    node2 = next(tree.find_clades(name=node2_name), None)

    if node1 is None or node2 is None:
        return 0

    return int(tree.distance(node1, node2))



def main():
    with open("Dataset.txt", "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    results = []
    i = 0

    while i < len(lines):
        tree_str = lines[i]
        i += 1

        node1_name, node2_name = lines[i].split()
        i += 1

        dist = distance_in_tree_biopython(tree_str, node1_name, node2_name)
        results.append(str(dist))

    print(" ".join(results))


if __name__ == "__main__":
    main()
