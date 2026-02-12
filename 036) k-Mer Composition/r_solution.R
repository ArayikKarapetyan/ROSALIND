kmer_composition_r <- function(dna, k = 4) {
  # Generate all possible k-mers in lexicographic order
  generate_kmers <- function(k) {
    bases <- c("A", "C", "G", "T")
    expand.grid(rep(list(bases), k)) %>% 
      apply(1, paste, collapse = "") %>%
      sort()
  }
  
  all_kmers <- generate_kmers(k)
  
  # Initialize count vector
  counts <- integer(length(all_kmers))
  names(counts) <- all_kmers
  
  # Count occurrences of each k-mer
  n <- nchar(dna)
  for (i in 1:(n - k + 1)) {
    kmer <- substr(dna, i, i + k - 1)
    if (kmer %in% all_kmers) {
      counts[kmer] <- counts[kmer] + 1
    }
  }
  
  return(as.vector(counts))
}

main <- function() {
  # Load required library
  if (!require("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
    library(dplyr)
  }
  
  # Read FASTA file
  lines <- readLines("36) k-Mer Composition\\Dataset.txt", warn = FALSE)
  
  # Extract DNA sequence
  dna <- ""
  in_sequence <- FALSE
  for (line in lines) {
    if (startsWith(line, ">")) {
      in_sequence <- TRUE
      next
    }
    if (in_sequence) {
      dna <- paste0(dna, line)
    }
  }
  
  # Calculate 4-mer composition
  composition <- kmer_composition_r(dna, 4)
  
  # Output result
  cat(paste(composition, collapse = " "), "\n")
}

# Run the program
main()