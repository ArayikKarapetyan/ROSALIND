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
  
  return(as.character(sequences))
}

longest_common_substring_r <- function(sequences) {
  if (length(sequences) == 0) return("")
  
  first_seq <- sequences[1]
  n <- nchar(first_seq)
  
  longest <- ""
  
  # Try all substrings of first sequence
  for (i in 1:n) {
    for (j in i:n) {
      substr <- substr(first_seq, i, j)
      
      # Check if in all other sequences
      in_all <- TRUE
      for (k in 2:length(sequences)) {
        if (!grepl(substr, sequences[k], fixed = TRUE)) {
          in_all <- FALSE
          break
        }
      }
      
      if (in_all && nchar(substr) > nchar(longest)) {
        longest <- substr
      }
    }
  }
  
  return(longest)
}

fasta_data <- readLines("C:\\Users\\arayi\\Desktop\\ROSALIND\\14) Finding a Shared Motif\\Dataset.txt", warn = FALSE)
fasta_string <- paste(fasta_data, collapse = "\n")
sequences <- parse_fasta_r(fasta_string)

print(paste("R Basic:", longest_common_substring_r(sequences)))
