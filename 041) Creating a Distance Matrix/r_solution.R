compute_p_distance <- function(s1, s2) {
  # Calculate p-distance between two strings
  chars1 <- strsplit(s1, "")[[1]]
  chars2 <- strsplit(s2, "")[[1]]
  
  if (length(chars1) != length(chars2)) {
    stop("Strings must have equal length")
  }
  
  mismatches <- sum(chars1 != chars2)
  return(mismatches / length(chars1))
}

compute_distance_matrix_r <- function(sequences) {
  # Compute p-distance matrix for list of sequences
  n <- length(sequences)
  matrix <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in i:n) {
      if (i == j) {
        matrix[i, j] <- 0
      } else {
        dist <- compute_p_distance(sequences[i], sequences[j])
        matrix[i, j] <- dist
        matrix[j, i] <- dist
      }
    }
  }
  
  return(matrix)
}

main <- function() {
  # Read FASTA file
  lines <- readLines("41) Creating a Distance Matrix\\Dataset.txt", warn = FALSE)
  
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
  
  # Compute distance matrix
  dist_matrix <- compute_distance_matrix_r(sequences)
  
  # Output with 5 decimal places
  for (i in 1:nrow(dist_matrix)) {
    row_values <- sprintf("%.5f", dist_matrix[i, ])
    cat(paste(row_values, collapse = " "), "\n")
  }
}

# Run the program
main()