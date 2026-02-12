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

overlap_length_r <- function(a, b, min_overlap = 1) {
  max_len <- min(nchar(a), nchar(b))
  
  # Try from largest to smallest overlap
  for (overlap in max_len:min_overlap) {
    a_end <- substring(a, nchar(a) - overlap + 1, nchar(a))
    b_start <- substring(b, 1, overlap)
    
    if (a_end == b_start) {
      return(overlap)
    }
  }
  
  return(0)
}

merge_strings_r <- function(a, b, overlap) {
  if (overlap == 0) {
    return(paste0(a, b))
  } else {
    return(paste0(a, substring(b, overlap + 1)))
  }
}

find_max_overlap_pair_r <- function(strings) {
  n <- length(strings)
  max_overlap <- 0
  best_i <- -1
  best_j <- -1
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) next
      
      overlap <- overlap_length_r(strings[i], strings[j])
      if (overlap > max_overlap) {
        max_overlap <- overlap
        best_i <- i
        best_j <- j
      }
    }
  }
  
  return(list(pair = c(best_i, best_j), overlap = max_overlap))
}

shortest_superstring_r <- function(strings) {
  current_strings <- strings
  
  while (length(current_strings) > 1) {
    result <- find_max_overlap_pair_r(current_strings)
    i <- result$pair[1]
    j <- result$pair[2]
    overlap <- result$overlap
    
    if (i == -1 || j == -1) {
      # No overlaps found, concatenate all
      return(paste(current_strings, collapse = ""))
    }
    
    merged <- merge_strings_r(current_strings[i], current_strings[j], overlap)
    
    # Remove the two strings and add merged one
    # Remove larger index first
    idx <- c(i, j)
    idx <- sort(idx, decreasing = TRUE)
    current_strings <- current_strings[-idx]
    current_strings <- c(current_strings, merged)
  }
  
  return(current_strings[1])
}


strings <- parse_fasta_r(paste(readLines("Dataset.txt"), collapse = "\n"))

result <- shortest_superstring_r(strings)
cat(result, "\n")