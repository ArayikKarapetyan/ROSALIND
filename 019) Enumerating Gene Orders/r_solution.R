generate_permutations_r <- function(n) {
  # Generate all permutations using expand.grid or permutations from gtools
  if (n <= 7) {
    # Simple approach for small n
    nums <- 1:n
    
    # Generate permutations using expand.grid (not efficient but works for n≤7)
    if (n <= 4) {
      # Simple expand.grid approach for very small n
      perms <- expand.grid(rep(list(nums), n))
      # Filter to keep only rows with all unique values
      perms <- perms[apply(perms, 1, function(x) length(unique(x)) == n), ]
    } else {
      # Use gtools package for larger n (within limit)
      if (!requireNamespace("gtools", quietly = TRUE)) {
        # Manual implementation if gtools not available
        perms <- manual_permutations_r(n)
      } else {
        perms <- gtools::permutations(n, n)
        # Convert from indices to actual numbers
        perms <- matrix(nums[perms], ncol = n)
      }
    }
    
    # Format output
    total <- nrow(perms)
    output <- character(total + 1)
    output[1] <- as.character(total)
    
    for (i in 1:total) {
      output[i + 1] <- paste(perms[i, ], collapse = " ")
    }
    
    return(paste(output, collapse = "\n"))
  } else {
    stop("n must be ≤ 7")
  }
}

# Manual permutations in R
manual_permutations_r <- function(n) {
  # Recursive function to generate permutations
  permute_recursive <- function(arr) {
    if (length(arr) == 1) {
      return(matrix(arr, nrow = 1))
    }
    
    results <- matrix(nrow = 0, ncol = length(arr))
    
    for (i in 1:length(arr)) {
      first <- arr[i]
      rest <- arr[-i]
      
      sub_perms <- permute_recursive(rest)
      for (j in 1:nrow(sub_perms)) {
        results <- rbind(results, c(first, sub_perms[j, ]))
      }
    }
    
    return(results)
  }
  
  return(permute_recursive(1:n))
}


n <- as.integer(readLines("Dataset.txt", warn = FALSE))


cat("R Solution:\n")
cat(generate_permutations_r(n), "\n")
