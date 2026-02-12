generate_strings_varying_length_r <- function(alphabet, n) {
  # Generate all strings of length at most n
  results <- character(0)
  
  # Recursive function
  generate <- function(current) {
    if (nchar(current) > 0) {
      results <<- c(results, current)
    }
    if (nchar(current) < n) {
      for (char in alphabet) {
        generate(paste0(current, char))
      }
    }
  }
  
  generate("")
  return(results)
}

main <- function() {
  # Read input from Dataset.txt
  lines <- readLines("39) Ordering Strings of Varying Length Lexicographically\\Dataset.txt", warn = FALSE)
  
  # Parse alphabet
  alphabet <- strsplit(lines[1], " ")[[1]]
  
  # Parse n
  n <- as.integer(lines[2])
  
  # Generate strings
  strings <- generate_strings_varying_length_r(alphabet, n)
  
  # Output results
  for (s in strings) {
    cat(s, "\n")
  }
}

# Run the program
main()