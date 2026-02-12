edit_distance_r <- function(s, t) {
  # Convert strings to character vectors
  s_chars <- strsplit(s, "")[[1]]
  t_chars <- strsplit(t, "")[[1]]
  m <- length(s_chars)
  n <- length(t_chars)
  
  # Create DP table
  dp <- matrix(0, nrow = m + 1, ncol = n + 1)
  
  # Initialize
  for (i in 1:(m + 1)) {
    dp[i, 1] <- i - 1
  }
  for (j in 1:(n + 1)) {
    dp[1, j] <- j - 1
  }
  
  # Fill DP table
  for (i in 2:(m + 1)) {
    for (j in 2:(n + 1)) {
      if (s_chars[i - 1] == t_chars[j - 1]) {
        sub_cost <- dp[i - 1, j - 1]
      } else {
        sub_cost <- dp[i - 1, j - 1] + 1
      }
      
      del_cost <- dp[i - 1, j] + 1
      ins_cost <- dp[i, j - 1] + 1
      
      dp[i, j] <- min(sub_cost, del_cost, ins_cost)
    }
  }
  
  return(dp[m + 1, n + 1])
}

main <- function() {
  # Read FASTA file
  lines <- readLines("Dataset.txt", warn = FALSE)
  
  # Parse sequences
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
  
  # Get two protein strings
  s <- sequences[1]
  t <- sequences[2]
  
  # Compute edit distance
  distance <- edit_distance_r(s, t)
  
  # Output result
  cat(distance, "\n")
}

# Run the program
main()