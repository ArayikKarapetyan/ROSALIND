compute_failure_array_r <- function(s) {
  # Compute failure array (prefix function) for string s
  n <- nchar(s)
  P <- integer(n)
  P[1] <- 0  # First element is 0 by convention
  
  for (i in 2:n) {
    j <- P[i-1] + 1  # Convert to 1-based indexing for R
    
    # Get characters (R is 1-based)
    char_i <- substr(s, i, i)
    
    # Find longest proper prefix that is also suffix
    while (j > 1) {
      char_j <- substr(s, j, j)
      if (char_i == char_j) {
        break
      }
      j <- P[j-1] + 1
    }
    
    # Check if we found a match
    if (j == 1) {
      char_j <- substr(s, 1, 1)
      if (char_i == char_j) {
        P[i] <- 1
      } else {
        P[i] <- 0
      }
    } else {
      P[i] <- j
    }
  }
  
  return(P)
}

main <- function() {
  # Read FASTA file
  lines <- readLines("Dataset.txt", warn = FALSE)
  
  # Extract DNA sequence
  dna <- ""
  in_seq <- FALSE
  for (line in lines) {
    if (startsWith(line, ">")) {
      in_seq <- TRUE
      next
    }
    if (in_seq) {
      dna <- paste0(dna, line)
    }
  }
  
  # Compute failure array
  failure_array <- compute_failure_array_r(dna)
  
  # Output result
  cat(paste(failure_array, collapse = " "), "\n")
}

# Run the program
main()