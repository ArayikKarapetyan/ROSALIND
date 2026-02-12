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

calculate_gc_content_r <- function(dna_string) {
  dna_chars <- strsplit(dna_string, "")[[1]]
  gc_count <- sum(dna_chars == "G" | dna_chars == "C")
  return((gc_count / length(dna_chars)) * 100)
}

find_highest_gc_r <- function(fasta_string) {
  sequences <- parse_fasta_r(fasta_string)
  
  highest_id <- ""
  highest_gc <- 0
  
  for (id in names(sequences)) {
    dna <- sequences[[id]]
    gc_content <- calculate_gc_content_r(dna)
    
    if (gc_content > highest_gc) {
      highest_gc <- gc_content
      highest_id <- id
    }
  }
  
  return(list(id = highest_id, gc_content = highest_gc))
}


input_file <- "Dataset.txt"  # Make sure the path is correct

file_content <- readLines(input_file, warn = FALSE)
data <- paste(file_content, collapse = "\n")


result <- find_highest_gc_r(data)
cat(sprintf("R result: %s\n%.6f\n", result$id, result$gc_content))

