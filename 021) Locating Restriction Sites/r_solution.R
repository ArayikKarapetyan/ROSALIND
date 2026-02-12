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

reverse_complement_r <- function(dna) {
  complement <- list(A = "T", T = "A", C = "G", G = "C")
  dna_chars <- strsplit(dna, "")[[1]]
  rev_chars <- rev(dna_chars)
  comp_chars <- sapply(rev_chars, function(x) complement[[x]])
  return(paste(comp_chars, collapse = ""))
}

find_reverse_palindromes_r <- function(dna, min_len = 4, max_len = 12) {
  results <- data.frame(position = integer(), length = integer())
  n <- nchar(dna)
  
  complement <- list(A = "T", T = "A", C = "G", G = "C")
  
  for (i in 1:n) {
    for (length in min_len:max_len) {
      if (i + length - 1 > n) break
      
      substring <- substr(dna, i, i + length - 1)
      
      # Check if palindrome
      is_palindrome <- TRUE
      half_len <- floor(length / 2)
      
      for (k in 1:half_len) {
        left <- substr(substring, k, k)
        right <- substr(substring, length - k + 1, length - k + 1)
        
        if (complement[[left]] != right) {
          is_palindrome <- FALSE
          break
        }
      }
      
      if (is_palindrome) {
        results <- rbind(results, data.frame(position = i, length = length))
      }
    }
  }
  
  return(results)
}

# More efficient R version
find_reverse_palindromes_efficient_r <- function(dna, min_len = 4, max_len = 12) {
  n <- nchar(dna)
  results <- list()
  
  complement <- c(A = "T", T = "A", C = "G", G = "C")
  
  # Split DNA into character vector for faster access
  dna_chars <- strsplit(dna, "")[[1]]
  
  for (i in 1:n) {
    max_possible <- min(max_len, n - i + 1)
    
    for (length in min_len:max_possible) {
      # Quick check: first and last characters must be complements
      if (complement[dna_chars[i]] != dna_chars[i + length - 1]) {
        next
      }
      
      # Check the rest
      is_palindrome <- TRUE
      for (k in 1:(floor(length / 2) - 1)) {
        if (complement[dna_chars[i + k]] != dna_chars[i + length - k]) {
          is_palindrome <- FALSE
          break
        }
      }
      
      if (is_palindrome) {
        results[[length(results) + 1]] <- c(i, length)
      }
    }
  }
  
  if (length(results) > 0) {
    results_matrix <- do.call(rbind, results)
    colnames(results_matrix) <- c("position", "length")
    return(results_matrix)
  } else {
    return(matrix(nrow = 0, ncol = 2, dimnames = list(NULL, c("position", "length"))))
  }
}




fasta_data <- paste(
  readLines("21) Locating Restriction Sites\\Dataset.txt", warn = FALSE),
  collapse = "\n"
)

sequences <- parse_fasta_r(fasta_data)
dna <- sequences[[1]]

results <- find_reverse_palindromes_r(dna)

for (i in 1:nrow(results)) {
    cat(results$position[i], results$length[i], "\n")
    }