spectral_convolution_r <- function(S1, S2) {
  # Compute spectral convolution S1 ⊖ S2
  # Returns table of differences and their counts
  
  # Generate all differences
  differences <- numeric(0)
  
  for (s1 in S1) {
    for (s2 in S2) {
      diff <- round(s1 - s2, 5)  # Round to handle precision
      differences <- c(differences, diff)
    }
  }
  
  # Count multiplicities
  conv_table <- table(differences)
  
  return(conv_table)
}

find_max_multiplicity_r <- function(conv_table) {
  # Find maximum multiplicity and corresponding x
  
  if (length(conv_table) == 0) {
    return(list(multiplicity = 0, x = 0.0))
  }
  
  # Find maximum multiplicity
  max_mult <- max(conv_table)
  
  # Find x values with this multiplicity
  max_x <- as.numeric(names(conv_table)[conv_table == max_mult])
  
  # Take absolute value and return first one
  abs_x <- abs(max_x)
  
  return(list(multiplicity = max_mult, x = abs_x[1]))
}

main <- function() {
  # Read input from Dataset.txt
  lines <- readLines("54) Comparing Spectra with the Spectral Convolution\\Dataset.txt", warn = FALSE)
  
  # Parse the two multisets
  S1 <- as.numeric(strsplit(lines[1], " ")[[1]])
  S2 <- as.numeric(strsplit(lines[2], " ")[[1]])
  
  # Compute spectral convolution
  conv_table <- spectral_convolution_r(S1, S2)
  
  # Find maximum multiplicity and x
  result <- find_max_multiplicity_r(conv_table)
  
  # Output results
  cat(result$multiplicity, "\n")
  cat(sprintf("%.5f", result$x), "\n")
}

# Run the program
main()