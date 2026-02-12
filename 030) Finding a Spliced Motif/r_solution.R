find_subsequence_indices_r <- function(s, t) {
  # Convert strings to character vectors
  s_chars <- strsplit(s, "")[[1]]
  t_chars <- strsplit(t, "")[[1]]
  
  indices <- integer(length(t_chars))
  s_pos <- 1
  t_pos <- 1
  
  while (t_pos <= length(t_chars) && s_pos <= length(s_chars)) {
    if (s_chars[s_pos] == t_chars[t_pos]) {
      indices[t_pos] <- s_pos
      t_pos <- t_pos + 1
    }
    s_pos <- s_pos + 1
  }
  
  if (t_pos > length(t_chars)) {
    return(indices)
  } else {
    return(NA)
  }
}

# FASTA parsing function for R
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
s <- sequences[[1]]
t <- sequences[[2]]
indices <- find_subsequence_indices_r(s, t)
cat("R result:", paste(indices, collapse = " "), "\n")