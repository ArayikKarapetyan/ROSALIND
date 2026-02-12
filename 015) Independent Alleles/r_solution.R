probability_at_least_n_r <- function(k, N) {
  # Total organisms in generation k
  n <- 2^k
  
  # Probability a single offspring is AaBb
  # Aa x Aa → 0.5 Aa, Bb x Bb → 0.5 Bb
  # Independent: 0.5 * 0.5 = 0.25
  p <- 0.25
  
  # Calculate probability of less than N using binomial distribution
  # P(X < N) = sum_{i=0}^{N-1} choose(n, i) * p^i * (1-p)^{n-i}
  prob_less_than_n <- 0
  
  for (i in 0:(N-1)) {
    prob_i <- choose(n, i) * (p^i) * ((1-p)^(n-i))
    prob_less_than_n <- prob_less_than_n + prob_i
  }
  
  # Probability of at least N
  prob_at_least_n <- 1 - prob_less_than_n
  return(prob_at_least_n)
}

# Using R's built-in binomial functions
probability_at_least_n_r_builtin <- function(k, N) {
  n <- 2^k
  p <- 0.25
  
  # pbinom gives cumulative probability P(X <= q)
  # So P(X < N) = P(X <= N-1)
  prob_less_than_n <- pbinom(N - 1, size = n, prob = p)
  
  prob_at_least_n <- 1 - prob_less_than_n
  return(prob_at_least_n)
}


input_data <- readLines('Dataset.txt', warn = FALSE)
values <- strsplit(input_data, " ")[[1]]
k <- as.integer(values[1])
N <- as.integer(values[2])

print(paste("R Basic:", round(probability_at_least_n_r(k, N), 6)))
print(paste("R Built-in:", round(probability_at_least_n_r_builtin(k, N), 6)))