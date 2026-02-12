count_rna_matchings_dp <- function(rna) {
  MOD <- 1000000
  n <- nchar(rna)
  bases <- strsplit(rna, "")[[1]]
  
  # Helper function
  can_pair <- function(b1, b2) {
    (b1 == "A" && b2 == "U") ||
    (b1 == "U" && b2 == "A") ||
    (b1 == "C" && b2 == "G") ||
    (b1 == "G" && b2 == "C")
  }
  
  # Initialize DP table
  dp <- matrix(0, n, n)
  
  for (length_sub in 0:(n-1)) {
    for (i in 1:(n - length_sub)) {
      j <- i + length_sub
      if (i > j) {
        dp[i,j] <- 1
      } else if ((j - i + 1) %% 2 != 0) {
        dp[i,j] <- 0
      } else if (i == j) {
        dp[i,j] <- 1
      } else {
        total <- 0
        for (k in seq(i+1, j, by=2)) {
          if (can_pair(bases[i], bases[k])) {
            left <- if (i+1 <= k-1) dp[i+1, k-1] else 1
            right <- if (k+1 <= j) dp[k+1, j] else 1
            total <- (total + left * right) %% MOD
          }
        }
        dp[i,j] <- total
      }
    }
  }
  
  return(dp[1,n])
}

# Main
main <- function() {
  lines <- readLines("Dataset.txt", warn=FALSE)
  rna <- paste0(lines[!startsWith(lines, ">")], collapse="")
  result <- count_rna_matchings_dp(rna)
  cat(result, "\n")
}

main()
