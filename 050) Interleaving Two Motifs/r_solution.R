shortest_common_supersequence_r <- function(s, t) {
  # Convert strings to character vectors
  s_chars <- strsplit(s, "")[[1]]
  t_chars <- strsplit(t, "")[[1]]
  m <- length(s_chars)
  n <- length(t_chars)
  
  # Compute LCS length table
  dp <- matrix(0, nrow = m + 1, ncol = n + 1)
  
  for (i in 1:m) {
    for (j in 1:n) {
      if (s_chars[i] == t_chars[j]) {
        dp[i + 1, j + 1] <- dp[i, j] + 1
      } else {
        dp[i + 1, j + 1] <- max(dp[i, j + 1], dp[i + 1, j])
      }
    }
  }
  
  # Reconstruct SCS
  scs_chars <- character(0)
  i <- m
  j <- n
  
  while (i > 0 && j > 0) {
    if (s_chars[i] == t_chars[j]) {
      scs_chars <- c(s_chars[i], scs_chars)
      i <- i - 1
      j <- j - 1
    } else if (dp[i, j + 1] > dp[i + 1, j]) {
      scs_chars <- c(s_chars[i], scs_chars)
      i <- i - 1
    } else {
      scs_chars <- c(t_chars[j], scs_chars)
      j <- j - 1
    }
  }
  
  # Add remaining characters
  while (i > 0) {
    scs_chars <- c(s_chars[i], scs_chars)
    i <- i - 1
  }
  while (j > 0) {
    scs_chars <- c(t_chars[j], scs_chars)
    j <- j - 1
  }
  
  return(paste(scs_chars, collapse = ""))
}

main <- function() {
  # Read input from Dataset.txt
  lines <- readLines("50) Interleaving Two Motifs\\Dataset.txt", warn = FALSE)
  
  # Get two DNA strings
  s <- lines[1]
  t <- lines[2]
  
  # Find shortest common supersequence
  scs <- shortest_common_supersequence_r(s, t)
  
  # Output result
  cat(scs, "\n")
}

# Run the program
main()