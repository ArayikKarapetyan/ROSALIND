library(ape)

parse_newick_r <- function(newick_str) {
  # Parse Newick string using ape package
  tree <- read.tree(text = newick_str)
  return(tree)
}

get_character_table_r <- function(tree) {
  # Get character table from phylogenetic tree
  
  # Get all leaves (taxa)
  all_leaves <- tree$tip.label
  all_leaves <- sort(all_leaves)  # Lexicographic order
  n_leaves <- length(all_leaves)
  
  # Initialize character table
  char_table <- character(0)
  
  # Function to process nodes
  process_node <- function(node_id, parent_leaves = NULL) {
    # node_id can be tip number or internal node number
    
    if (is.null(parent_leaves)) {
      parent_leaves <- all_leaves
    }
    
    # Get leaves in this subtree
    if (node_id <= length(tree$tip.label)) {
      # Tip node
      node_leaves <- tree$tip.label[node_id]
    } else {
      # Internal node
      # Get all tips under this node
      node_leaves <- extract.clade(tree, node_id)$tip.label
    }
    
    # Check if split is non-trivial
    n_node <- length(node_leaves)
    n_parent <- length(parent_leaves)
    
    if (n_node > 1 && n_node < n_parent) {
      # Create binary character row
      row <- rep("0", n_leaves)
      for (leaf in node_leaves) {
        idx <- which(all_leaves == leaf)
        row[idx] <- "1"
      }
      char_table <<- c(char_table, paste(row, collapse = ""))
    }
    
    # Get children of this node
    # In ape, edge matrix has parent in column 1, child in column 2
    children <- tree$edge[tree$edge[, 1] == node_id, 2]
    
    # Process children
    for (child in children) {
      process_node(child, node_leaves)
    }
  }
  
  # Start from root (node with no parent)
  # Root is the node not appearing in child column
  all_nodes <- unique(c(tree$edge))
  child_nodes <- unique(tree$edge[, 2])
  root <- setdiff(all_nodes, child_nodes)
  
  process_node(root)
  
  return(char_table)
}

main <- function() {
  # Read input from Dataset.txt
  newick <- readLines("55) Creating a Character Table\\Dataset.txt", warn = FALSE)[1]
  
  # Parse tree
  tree <- parse_newick_r(newick)
  
  # Get character table
  char_table <- get_character_table_r(tree)
  
  # Output table
  for (row in char_table) {
    cat(row, "\n")
  }
}

# Run the program
main()