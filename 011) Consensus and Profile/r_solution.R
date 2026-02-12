parse_fasta_r <- function(fasta_string) {
  lines <- strsplit(fasta_string, "\n")[[1]]
  sequences <- list()
  current_id <- ""
  
  for (line in lines) {
    if (startsWith(line, ">")) {
      current_id <- substring(line, 2)
      sequences[[current_id]] <- ""
    } else {
      sequences[[current_id]] <- paste0(sequences[[current_id]], line)
    }
  }
  
  return(sequences)
}

compute_profile_and_consensus_r <- function(sequences) {
  # Convert sequences to matrix of characters
  seq_list <- as.character(sequences)
  n <- nchar(seq_list[1])  # sequence length
  m <- length(seq_list)    # number of sequences
  
  # Initialize profile matrix
  profile <- list(
    A = integer(n),
    C = integer(n),
    G = integer(n),
    T = integer(n)
  )
  
  names(profile) <- c("A", "C", "G", "T")
  
  # Fill profile matrix
  for (seq in seq_list) {
    chars <- strsplit(seq, "")[[1]]
    for (j in 1:n) {
      nucleotide <- chars[j]
      profile[[nucleotide]][j] <- profile[[nucleotide]][j] + 1
    }
  }
  
  # Compute consensus string
  consensus_chars <- character(n)
  for (j in 1:n) {
    counts <- c(profile$A[j], profile$C[j], profile$G[j], profile$T[j])
    max_idx <- which.max(counts)
    consensus_chars[j] <- c("A", "C", "G", "T")[max_idx]
  }
  
  consensus <- paste(consensus_chars, collapse = "")
  
  return(list(consensus = consensus, profile = profile))
}

fasta_data <- readLines("C:\\Users\\arayi\\Desktop\\ROSALIND\\11) Consensus and Profile\\Dataset.txt", warn = FALSE)
fasta_string <- paste(fasta_data, collapse = "\n")
sequences <- parse_fasta_r(fasta_string)
result <- compute_profile_and_consensus_r(sequences)

cat("Consensus:", result$consensus, "\n")
cat("\nProfile Matrix:\n")
for (nuc in c("A", "C", "G", "T")) {
  cat(nuc, ": ", paste(result$profile[[nuc]], collapse = " "), "\n")
  }