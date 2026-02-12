probability_at_least_one_match_r <- function(N, gc_content, motif) {
  # Calculate probability of single string matching motif
  prob_single <- 1.0
  
  for (base in strsplit(motif, "")[[1]]) {
    if (base %in% c("G", "C")) {
      prob_single <- prob_single * (gc_content / 2)
    } else {
      prob_single <- prob_single * ((1 - gc_content) / 2)
    }
  }
  
  # Probability that at least one of N strings matches
  if (prob_single == 0) {
    return(0.0)
  } else if (prob_single == 1) {
    return(1.0)
  } else {
    prob_none <- (1 - prob_single) ^ N
    return(1 - prob_none)
  }
}

main <- function() {
  # Read input
  lines <- readLines("Dataset.txt", warn = FALSE)
  
  # Parse first line: N and GC content
  first_line <- strsplit(lines[1], " ")[[1]]
  N <- as.integer(first_line[1])
  gc_content <- as.numeric(first_line[2])
  
  # Parse motif
  motif <- lines[2]
  
  # Calculate probability
  prob <- probability_at_least_one_match_r(N, gc_content, motif)
  
  # Output with 3 decimal places
  cat(sprintf("%.3f\n", prob))
}

# Run the program
main()