rabbit_pairs_r <- function(n, k) {
  # Base cases
  if (n <= 2) {
    return(1)
  }
  
  # Initialize DP vector
  dp <- numeric(n)
  dp[1] <- 1
  dp[2] <- 1
  
  # Fill DP table
  for (i in 3:n) {
    dp[i] <- dp[i-1] + k * dp[i-2]
  }
  
  return(dp[n])
}

# Alternative functional approach
rabbit_pairs_r_functional <- function(n, k) {
  # Base cases
  if (n <= 2) return(1)
  
  # Recursive with memoization
  memo <- new.env()
  memo$fib <- function(n) {
    if (n <= 2) return(1)
    if (exists(as.character(n), envir = memo)) {
      return(get(as.character(n), envir = memo))
    }
    result <- memo$fib(n-1) + k * memo$fib(n-2)
    assign(as.character(n), result, envir = memo)
    return(result)
  }
  
  return(memo$fib(n))
}

# R usage
input_data <- readLines("Dataset.txt", warn = FALSE)
values <- strsplit(input_data, " ")[[1]]
n <- as.integer(values[1])
k <- as.integer(values[2])
print(paste("R DP result:", rabbit_pairs_r(n, k)))
print(paste("R functional result:", rabbit_pairs_r_functional(n, k)))