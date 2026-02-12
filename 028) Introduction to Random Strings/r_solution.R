calculate_log_probabilities_r <- function(dna_string, gc_contents) {
  # Split DNA string into characters
  dna_chars <- strsplit(dna_string, "")[[1]]
  
  # Count A/T and C/G
  at_count <- sum(dna_chars == "A" | dna_chars == "T")
  cg_count <- sum(dna_chars == "C" | dna_chars == "G")
  
  results <- numeric(length(gc_contents))
  
  for (i in seq_along(gc_contents)) {
    gc <- gc_contents[i]
    prob_at <- (1 - gc) / 2
    prob_cg <- gc / 2
    
    # log10 probability
    results[i] <- at_count * log10(prob_at) + cg_count * log10(prob_cg)
  }
  
  return(results)
}

input_data <- readLines("28) Introduction to Random Strings/Dataset.txt", warn = FALSE)

dna_string <- input_data[1]
gc_contents <- as.numeric(strsplit(input_data[2], " ")[[1]])


result <- round(calculate_log_probabilities_r(dna_string, gc_contents), 3)

cat(result, sep = " ")
