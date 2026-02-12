dominant_probability_r <- function(k, m, n) {
  total <- k + m + n
  
  # Calculate probabilities for all pair types
  total_pairs <- total * (total - 1)
  
  # Calculate number of dominant offspring from each pair type
  dominant_count <- 
    k * (k - 1) * 1 +           # KK: always dominant
    2 * k * m * 1 +             # KM or MK: always dominant
    2 * k * n * 1 +             # KN or NK: always dominant
    m * (m - 1) * 0.75 +        # MM: 75% dominant
    2 * m * n * 0.5 +           # MN or NM: 50% dominant
    n * (n - 1) * 0             # NN: never dominant
  
  probability <- dominant_count / total_pairs
  return(probability)
}

# R usage
input_data <- readLines("Dataset.txt", warn = FALSE)
values <- strsplit(input_data, " ")[[1]]
k <- as.integer(values[1])
m <- as.integer(values[2])
n <- as.integer(values[3])
print(paste("R:", round(dominant_probability_r(k, m, n), 5)))
