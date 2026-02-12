count_rna_matchings_fast <- function(rna) {
  MOD <- 1000000
  n <- nchar(rna)
  bases <- strsplit(rna, "")[[1]]

  can_pair <- function(b1, b2) {
    (b1 == "A" && b2 == "U") ||
    (b1 == "U" && b2 == "A") ||
    (b1 == "C" && b2 == "G") ||
    (b1 == "G" && b2 == "C")
  }

  pair_ok <- matrix(FALSE, n, n)
  for (i in 1:n)
    for (j in 1:n)
      pair_ok[i, j] <- can_pair(bases[i], bases[j])

  dp <- matrix(0, n, n)

  for (i in 1:n)
    dp[i, i] <- 1

  for (len in 2:n) {
    for (i in 1:(n - len + 1)) {
      j <- i + len - 1

      dp[i, j] <- dp[i + 1, j]

      if (i + 1 <= j) {
        for (k in (i + 1):j) {
          if (pair_ok[i, k]) {
            left  <- if (i + 1 <= k - 1) dp[i + 1, k - 1] else 1
            right <- if (k + 1 <= j) dp[k + 1, j] else 1
            dp[i, j] <- (dp[i, j] + left * right) %% MOD
          }
        }
      }
    }
  }

  return(dp[1, n])
}

main <- function() {
  lines <- readLines("Dataset.txt",
                     warn = FALSE)

  rna <- paste(lines[!startsWith(lines, ">")], collapse = "")
  result <- count_rna_matchings_fast(rna)
  cat(result, "\n")
}

main()
