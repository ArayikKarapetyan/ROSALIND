library(gmp)

count_maximum_matchings_r <- function(rna) {
  bases <- strsplit(rna, "")[[1]]
  counts <- table(bases)
  
  A <- ifelse(is.na(counts["A"]), 0, counts["A"])
  U <- ifelse(is.na(counts["U"]), 0, counts["U"])
  C <- ifelse(is.na(counts["C"]), 0, counts["C"])
  G <- ifelse(is.na(counts["G"]), 0, counts["G"])
  
  # Convert to big integers
  A <- as.bigz(A); U <- as.bigz(U)
  C <- as.bigz(C); G <- as.bigz(G)
  
  # Permutations
  au_ways <- factorial(max(A, U)) / factorial(abs(A - U))
  cg_ways <- factorial(max(C, G)) / factorial(abs(C - G))
  
  return(au_ways * cg_ways)
}

main <- function() {
  lines <- readLines("Dataset.txt", warn = FALSE)
  
  rna <- ""
  for (line in lines) {
    if (!startsWith(line, ">")) {
      rna <- paste0(rna, line)
    }
  }
  
  result <- count_maximum_matchings_r(rna)
  
  # Correct output
  cat(as.character(result), "\n")
}

main()

