count_internal_nodes_r <- function(n_leaves) {
  # Calculate number of internal nodes in unrooted binary tree
  # Formula: internal_nodes = leaves - 2
  
  return(n_leaves - 2)
}

main <- function() {
  # Read input from Dataset.txt
  input <- readLines("Dataset.txt", warn = FALSE)
  n <- as.integer(input[1])
  
  # Calculate result
  result <- count_internal_nodes_r(n)
  
  # Output result
  cat(result, "\n")
}

# Run the program
main()