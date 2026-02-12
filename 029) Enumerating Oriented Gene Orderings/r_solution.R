library(combinat)

generate_signed_permutations_r <- function(n) {
  # Calculate total count directly: n! * 2^n
  total <- factorial(n) * (2^n)
  
  # Generate numbers
  numbers <- 1:n
  
  # Generate all permutations
  perms <- combinat::permn(numbers)
  
  # Generate all sign combinations
  sign_grid <- as.matrix(expand.grid(replicate(n, c(-1, 1), simplify = FALSE)))
  
  # Create all combinations
  all_perms <- list()
  
  for (perm in perms) {
    perm_matrix <- matrix(rep(perm, each = nrow(sign_grid)), nrow = nrow(sign_grid))
    signed_matrix <- perm_matrix * sign_grid
    for (i in 1:nrow(signed_matrix)) {
      all_perms <- c(all_perms, list(signed_matrix[i, ]))
    }
  }
  
  return(list(count = total, permutations = all_perms))
}



n <- as.integer(readLines("Dataset.txt", warn = FALSE))
result <- generate_signed_permutations_r(n)

cat(result$count, "\n")
for (perm in result$permutations) {
  cat(paste(perm, collapse = " "), "\n")
}