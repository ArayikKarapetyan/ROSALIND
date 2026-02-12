transcribe_dna_to_rna <- function(dna_string) {
  # Using gsub to replace T with U
  rna_string <- gsub("T", "U", dna_string)
  return(rna_string)
}

# Most concise R version
transcribe_dna_to_rna_concise <- function(dna_string) {
  chartr("T", "U", dna_string)
}


input_file <- "Dataset.txt"  # Make sure the path is correct

# Read the entire file
file_content <- readLines(input_file, warn = FALSE)

# Combine all lines into a single string
data <- paste(file_content, collapse = "")

# Clean up: remove any whitespace or newlines
data <- gsub("\\s", "", data)

print(paste("R result:", transcribe_dna_to_rna(data)))
# print(paste("R concise result:", transcribe_dna_to_rna_concise(data)))