calculate_tt_ratio_r <- function(s1, s2) {
  # Convert to character vectors
  s1_chars <- strsplit(s1, "")[[1]]
  s2_chars <- strsplit(s2, "")[[1]]
  
  transitions <- 0
  transversions <- 0
  
  # Define transition pairs
  transition_pairs <- list(
    c("A", "G"), c("G", "A"),
    c("C", "T"), c("T", "C")
  )
  
  for (i in seq_along(s1_chars)) {
    if (s1_chars[i] == s2_chars[i]) next
    
    pair <- c(s1_chars[i], s2_chars[i])
    
    # Check if it's a transition
    is_transition <- any(sapply(transition_pairs, function(x) all(x == pair)))
    
    if (is_transition) {
      transitions <- transitions + 1
    } else {
      transversions <- transversions + 1
    }
  }
  
  if (transversions == 0) {
    return(Inf)
  }
  
  return(transitions / transversions)
}


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


fasta_data <- readLines("Dataset.txt", warn = FALSE)
fasta_string <- paste(fasta_data, collapse = "\n")
sequences <- parse_fasta_r(fasta_string)
s1 <- sequences[[1]]
s2 <- sequences[[2]]

ratio <- calculate_tt_ratio_r(s1, s2)
cat(sprintf("R result: %.11f\n", ratio))
