sum_combinations_r <- function(n, m) {
  MOD <- 1000000
  
  # Initialize Pascal's triangle row
  row <- integer(n + 1)
  row[1] <- 1  # Index 1 corresponds to C(0,0) in R (1-based indexing)
  
  # Build Pascal's triangle up to row n
  for (i in 1:n) {
    # Update from right to left to avoid overwriting values
    for (j in (i + 1):2) {
      row[j] <- (row[j - 1] + row[j]) %% MOD
    }
    row[1] <- 1
  }
  
  # Sum from m to n (adjust for 1-based indexing)
  total <- 0
  for (k in (m + 1):(n + 1)) {
    total <- (total + row[k]) %% MOD
  }
  
  return(total)
}

main <- function() {
  # Read input from Dataset.txt
  input <- readLines("45) Introduction to Alternative Splicing\\Dataset.txt", warn = FALSE)
  numbers <- strsplit(input[1], " ")[[1]]
  n <- as.integer(numbers[1])
  m <- as.integer(numbers[2])
  
  # Calculate result
  result <- sum_combinations_r(n, m)
  
  # Output
  cat(result, "\n")
}

# Run the program
main()