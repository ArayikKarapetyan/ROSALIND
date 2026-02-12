# ===============================
# Parse input
# ===============================
parse_input_r <- function(input_string) {
  lines <- strsplit(input_string, "\n")[[1]]
  n <- as.integer(lines[1])
  permutation <- as.integer(strsplit(lines[2], " ")[[1]])
  list(n = n, permutation = permutation)
}

# ===============================
# Longest Increasing Subsequence (O(n^2))
# ===============================
longest_increasing_subsequence_r <- function(seq) {
  n <- length(seq)
  if (n == 0) return(integer(0))
  
  dp <- rep(1, n)
  prev <- rep(-1, n)
  
  for (i in seq_len(n)) {
    for (j in seq_len(i - 1)) {
      if (seq[j] < seq[i] && dp[j] + 1 > dp[i]) {
        dp[i] <- dp[j] + 1
        prev[i] <- j
      }
    }
  }
  
  # Reconstruct LIS
  lis <- integer(0)
  curr <- which.max(dp)
  while (curr != -1) {
    lis <- c(seq[curr], lis)
    curr <- prev[curr]
  }
  
  lis
}

# ===============================
# Longest Decreasing Subsequence (O(n^2))
# ===============================
longest_decreasing_subsequence_r <- function(seq) {
  n <- length(seq)
  if (n == 0) return(integer(0))
  
  dp <- rep(1, n)
  prev <- rep(-1, n)
  
  for (i in seq_len(n)) {
    for (j in seq_len(i - 1)) {
      if (seq[j] > seq[i] && dp[j] + 1 > dp[i]) {
        dp[i] <- dp[j] + 1
        prev[i] <- j
      }
    }
  }
  
  # Reconstruct LDS
  lds <- integer(0)
  curr <- which.max(dp)
  while (curr != -1) {
    lds <- c(seq[curr], lds)
    curr <- prev[curr]
  }
  
  lds
}

# ===============================
# Main processing
# ===============================
process_permutation_r <- function(input_string) {
  data <- parse_input_r(input_string)
  permutation <- data$permutation
  
  lis <- longest_increasing_subsequence_r(permutation)
  lds <- longest_decreasing_subsequence_r(permutation)
  
  list(lis = lis, lds = lds)
}

lines <- readLines("24) Longest Increasing Subsequence\\Dataset.txt", warn = FALSE)

input_string <- paste(lines, collapse = "\n")

result <- process_permutation_r(input_string)

cat(paste(result$lis, collapse = " "), "\n")
cat(paste(result$lds, collapse = " "), "\n")
