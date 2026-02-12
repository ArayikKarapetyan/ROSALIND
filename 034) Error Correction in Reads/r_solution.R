library(stringr)

# Function to compute reverse complement
reverse_complement <- function(dna) {
  bases <- strsplit(dna, "")[[1]]
  rev_bases <- rev(bases)
  comp <- chartr("ACGT", "TGCA", rev_bases)
  paste(comp, collapse = "")
}

# Function to compute Hamming distance
hamming_distance <- function(s1, s2) {
  chars1 <- strsplit(s1, "")[[1]]
  chars2 <- strsplit(s2, "")[[1]]
  sum(chars1 != chars2)
}

# Main error correction function
correct_errors_r <- function(reads) {
  # Count frequencies of reads and their reverse complements
  freq <- table(reads)
  
  # Add reverse complements to frequency count
  all_seqs <- c(reads, sapply(reads, reverse_complement))
  total_freq <- table(all_seqs)
  
  # Identify correct reads (total count >= 2)
  correct <- names(total_freq)[total_freq >= 2]
  
  # Find corrections
  corrections <- character(0)
  
  for (read in reads) {
    if (read %in% correct) {
      next
    }
    
    found_correction <- NULL
    ambiguous <- FALSE
    
    # Generate all single mutations
    chars <- strsplit(read, "")[[1]]
    n <- length(chars)
    
    for (i in 1:n) {
      original_base <- chars[i]
      
      for (base in c("A", "C", "G", "T")) {
        if (base == original_base) next
        
        # Create mutated sequence
        mutated <- chars
        mutated[i] <- base
        candidate <- paste(mutated, collapse = "")
        
        # Check if candidate is correct
        if (candidate %in% correct) {
          # Verify Hamming distance
          if (hamming_distance(read, candidate) == 1) {
            if (!is.null(found_correction)) {
              ambiguous <- TRUE
              break
            }
            found_correction <- candidate
          }
        }
      }
      
      if (ambiguous) break
    }
    
    if (!ambiguous && !is.null(found_correction)) {
      corrections <- c(corrections, paste0(read, "->", found_correction))
    }
  }
  
  return(corrections)
}

# Main function
main <- function() {
  # Read FASTA file
  lines <- readLines("Dataset.txt", warn = FALSE)
  
  reads <- character(0)
  current_seq <- ""
  
  for (line in lines) {
    if (str_starts(line, ">")) {
      if (nchar(current_seq) > 0) {
        reads <- c(reads, current_seq)
        current_seq <- ""
      }
    } else {
      current_seq <- line
    }
  }
  
  if (nchar(current_seq) > 0) {
    reads <- c(reads, current_seq)
  }
  
  # Correct errors
  corrections <- correct_errors_r(reads)
  
  # Output results
  for (corr in corrections) {
    cat(corr, "\n")
  }
}

# Run the program
main()