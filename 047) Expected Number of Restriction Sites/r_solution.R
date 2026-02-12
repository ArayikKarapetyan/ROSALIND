expected_occurrences_r <- function(n, s, gc_contents) {
  results <- numeric(length(gc_contents))
  length_s <- nchar(s)
  
  # Count GC and AT bases
  bases <- strsplit(s, "")[[1]]
  gc_count <- sum(bases %in% c("G", "C"))
  at_count <- length_s - gc_count
  
  # Calculate for each GC-content
  for (i in seq_along(gc_contents)) {
    gc_content <- gc_contents[i]
    
    # Probability of generating s at a specific position
    prob_s <- (gc_content / 2)^gc_count * ((1 - gc_content) / 2)^at_count
    
    # Expected number of occurrences
    expected <- prob_s * (n - length_s + 1)
    results[i] <- expected
  }
  
  return(results)
}

main <- function() {
  # Read input from Dataset.txt
  lines <- readLines("47) Expected Number of Restriction Sites\\Dataset.txt", warn = FALSE)
  lines <- lines[lines != ""]  # Remove empty lines
  
  n <- as.integer(lines[1])
  s <- lines[2]
  gc_contents <- as.numeric(strsplit(lines[3], " ")[[1]])
  
  # Calculate expected occurrences
  results <- expected_occurrences_r(n, s, gc_contents)
  
  # Format output with 3 decimal places
  formatted <- sprintf("%.3f", results)
  cat(paste(formatted, collapse = " "), "\n")
}

# Run the program
main()