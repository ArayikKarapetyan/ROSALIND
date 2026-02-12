parse_set_r <- function(set_str) {
  # Parse set string like "{1, 2, 3, 4, 5}"
  if (nchar(set_str) == 0) {
    return(integer(0))
  }
  
  # Remove braces
  set_str <- gsub("[{}]", "", set_str)
  
  if (nchar(trimws(set_str)) == 0) {
    return(integer(0))
  }
  
  # Split by comma and convert to integers
  elements <- as.integer(strsplit(set_str, ",")[[1]])
  return(elements)
}

format_set_r <- function(vec) {
  # Format integer vector as set string
  if (length(vec) == 0) {
    return("{}")
  }
  
  sorted_vec <- sort(vec)
  elements <- paste(sorted_vec, collapse = ", ")
  return(paste0("{", elements, "}"))
}

set_operations_r <- function(n, A_vec, B_vec) {
  # Convert to sets (use logical vectors for efficiency)
  A_set <- logical(n)
  if(length(A_vec) > 0) A_set[A_vec] <- TRUE
  
  B_set <- logical(n)
  if(length(B_vec) > 0) B_set[B_vec] <- TRUE
  
  # Universal set
  U_set <- rep(TRUE, n)
  
  # Perform operations
  union_set <- which(A_set | B_set)
  intersect_set <- which(A_set & B_set)
  A_minus_B_set <- which(A_set & !B_set)
  B_minus_A_set <- which(B_set & !A_set)
  A_complement_set <- which(!A_set & U_set)
  B_complement_set <- which(!B_set & U_set)
  
  return(list(
    union = union_set,
    intersect = intersect_set,
    A_minus_B = A_minus_B_set,
    B_minus_A = B_minus_A_set,
    A_complement = A_complement_set,
    B_complement = B_complement_set
  ))
}

main <- function() {
  # Read input from Dataset.txt
  lines <- readLines("51) Introduction to Set Operations\\Dataset.txt", warn = FALSE)
  lines <- lines[lines != ""]  # Remove empty lines
  
  n <- as.integer(lines[1])
  A_vec <- parse_set_r(lines[2])
  B_vec <- parse_set_r(lines[3])
  
  # Perform set operations
  results <- set_operations_r(n, A_vec, B_vec)
  
  # Prepare output lines
  output_lines <- c(
    format_set_r(results$union),
    format_set_r(results$intersect),
    format_set_r(results$A_minus_B),
    format_set_r(results$B_minus_A),
    format_set_r(results$A_complement),
    format_set_r(results$B_complement)
  )
  
  # Write to Output.txt
  writeLines(output_lines, "output.txt")
}

# Run the program
main()
