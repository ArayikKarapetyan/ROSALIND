library(memoise)

longest_common_subsequence_r <- function(s, t) {
  # Convert strings to character vectors
  s_chars <- strsplit(s, "")[[1]]
  t_chars <- strsplit(t, "")[[1]]
  m <- length(s_chars)
  n <- length(t_chars)
  
  # Create DP table
  dp <- matrix(0, nrow = m + 1, ncol = n + 1)
  
  # Fill DP table
  for (i in 1:m) {
    for (j in 1:n) {
      if (s_chars[i] == t_chars[j]) {
        dp[i + 1, j + 1] <- dp[i, j] + 1
      } else {
        dp[i + 1, j + 1] <- max(dp[i, j + 1], dp[i + 1, j])
      }
    }
  }
  
  # Backtrack to reconstruct LCS
  lcs_chars <- character(0)
  i <- m
  j <- n
  
  while (i > 0 && j > 0) {
    if (s_chars[i] == t_chars[j]) {
      lcs_chars <- c(s_chars[i], lcs_chars)
      i <- i - 1
      j <- j - 1
    } else if (dp[i, j + 1] > dp[i + 1, j]) {
      i <- i - 1
    } else {
      j <- j - 1
    }
  }
  
  return(paste(lcs_chars, collapse = ""))
}

main <- function() {
  # Read FASTA file
  lines <- readLines("Dataset.txt", warn = FALSE)
  
  # Extract sequences
  sequences <- character(0)
  current_seq <- ""
  
  for (line in lines) {
    if (startsWith(line, ">")) {
      if (nchar(current_seq) > 0) {
        sequences <- c(sequences, current_seq)
        current_seq <- ""
      }
    } else {
      current_seq <- paste0(current_seq, line)
    }
  }
  
  if (nchar(current_seq) > 0) {
    sequences <- c(sequences, current_seq)
  }
  
  # Get two sequences
  s <- sequences[1]
  t <- sequences[2]
  
  # Find LCS
  lcs <- longest_common_subsequence_r(s, t)
  
  # Output result
  cat(lcs, "\n")
}

# Run the program
main()