expected_dominant_offspring_r <- function(counts) {
  # Unpack the six counts
  AA_AA <- counts[1]
  AA_Aa <- counts[2]
  AA_aa <- counts[3]
  Aa_Aa <- counts[4]
  Aa_aa <- counts[5]
  aa_aa <- counts[6]
  
  # Probabilities of dominant phenotype for each mating type
  p_AA_AA <- 1.0
  p_AA_Aa <- 1.0
  p_AA_aa <- 1.0
  p_Aa_Aa <- 0.75
  p_Aa_aa <- 0.5
  p_aa_aa <- 0.0
  
  # Each couple produces exactly 2 offspring
  offspring_per_couple <- 2
  
  # Calculate expected dominant offspring
  expected <- (
    AA_AA * p_AA_AA * offspring_per_couple +
    AA_Aa * p_AA_Aa * offspring_per_couple +
    AA_aa * p_AA_aa * offspring_per_couple +
    Aa_Aa * p_Aa_Aa * offspring_per_couple +
    Aa_aa * p_Aa_aa * offspring_per_couple +
    aa_aa * p_aa_aa * offspring_per_couple
  )
  
  return(expected)
}

input_data <- readLines("Dataset.txt", warn = FALSE)
counts <- as.integer(strsplit(input_data, " ")[[1]])

print(paste("R Basic:", expected_dominant_offspring_r(counts)))