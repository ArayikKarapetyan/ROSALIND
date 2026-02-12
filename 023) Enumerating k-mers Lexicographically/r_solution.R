generate_lexicographic_r <- function(alphabet, n) {
  # Use expand.grid to generate all combinations
  # Create list of n copies of the alphabet
  args <- rep(list(alphabet), n)
  
  # Generate all combinations
  combos <- do.call(expand.grid, args)
  
  # Convert to strings
  strings <- apply(combos, 1, paste, collapse = "")
  
  # Sort (expand.grid produces in specific order but we sort to be sure)
  strings <- sort(strings)
  
  return(strings)
}


lines <- readLines("23) Enumerating k-mers Lexicographically\\Dataset.txt", warn = FALSE)
alphabet <- strsplit(lines[1], " ")[[1]]
n <- as.integer(lines[2])

strings <- generate_lexicographic_r(alphabet, n)
for (s in strings) {
  cat(s, "\n")
}
