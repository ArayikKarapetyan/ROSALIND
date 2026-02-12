reverse_complement_r <- function(dna_string) {
  # Split string into characters
  dna_chars <- strsplit(dna_string, "")[[1]]
  
  # Reverse the characters
  reversed_chars <- rev(dna_chars)
  
  # Create complement mapping
  complement <- function(base) {
    switch(base,
           "A" = "T",
           "T" = "A",
           "C" = "G",
           "G" = "C",
           base)
  }
  
  # Apply complement to each character
  comp_chars <- sapply(reversed_chars, complement)
  
  # Combine back to string
  return(paste(comp_chars, collapse = ""))
}

# Alternative concise version using chartr
reverse_complement_concise <- function(dna_string) {
  # Reverse the string
  reversed_string <- paste(rev(strsplit(dna_string, "")[[1]]), collapse = "")
  
  # Translate A<->T and C<->G
  complement_string <- chartr("ATCG", "TAGC", reversed_string)
  
  return(complement_string)
}


input_file <- "Dataset.txt"  # Make sure the path is correct

# Read the entire file
file_content <- readLines(input_file, warn = FALSE)

# Combine all lines into a single string
data <- paste(file_content, collapse = "")

# Clean up: remove any whitespace or newlines
data <- gsub("\\s", "", data)

print(paste("R result:", reverse_complement_r(data)))
# print(paste("R concise:", reverse_complement_concise(data)))