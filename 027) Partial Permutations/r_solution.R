partial_permutations_r <- function(n, k, mod = 1000000) {
  if (k > n) return(0)
  
  result <- 1
  # Calculate n × (n-1) × ... × (n-k+1) mod 1,000,000
  for (i in (n - k + 1):n) {
    result <- (result * i) %% mod
  }
  
  return(result)
}

input_data <- readLines("Dataset.txt", warn = FALSE)
values <- strsplit(input_data, " ")[[1]]
n <- as.integer(values[1])
k <- as.integer(values[2])
print(paste("R Basic:", partial_permutations_r(n, k)))