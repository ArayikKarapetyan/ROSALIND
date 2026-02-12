count_subsets_r <- function(n) {
  MOD <- 1000000
  
  # Use modular exponentiation
  result <- 1
  base <- 2
  
  while (n > 0) {
    if (n %% 2 == 1) {
      result <- (result * base) %% MOD
    }
    base <- (base * base) %% MOD
    n <- n %/% 2
  }
  
  return(result)
}

main <- function() {
  # Read input from Dataset.txt
  n <- as.integer(readLines("44) Counting Subsets\\Dataset.txt", warn = FALSE)[1])
  
  # Calculate number of subsets
  result <- count_subsets_r(n)
  
  # Output result
  cat(result, "\n")
}

# Run the program
main()