hamming_distance_r <- function(s, t) {
  # Check lengths
  if (nchar(s) != nchar(t)) {
    stop("Strings must be of equal length")
  }
  
  # Split strings into characters
  s_chars <- strsplit(s, "")[[1]]
  t_chars <- strsplit(t, "")[[1]]
  
  # Calculate mismatches
  mismatches <- sum(s_chars != t_chars)
  return(mismatches)
}

# Alternative using mapply
hamming_distance_r_mapply <- function(s, t) {
  if (nchar(s) != nchar(t)) {
    stop("Strings must be of equal length")
  }
  
  s_chars <- strsplit(s, "")[[1]]
  t_chars <- strsplit(t, "")[[1]]
  
  # Use mapply for pairwise comparison
  mismatches <- sum(mapply(function(x, y) x != y, s_chars, t_chars))
  return(mismatches)
}

input_file <- "Dataset.txt"  # Make sure the path is correct

# Read the entire file
file_content <- readLines(input_file, warn = FALSE)

# Combine all lines into a single string
data <- trimws(file_content)
# Clean up: remove any whitespace or newlines
s <- data[1]
t <- data[2]

print(paste("R Basic:", hamming_distance_r(s, t)))
print(paste("R mapply:", hamming_distance_r_mapply(s, t)))