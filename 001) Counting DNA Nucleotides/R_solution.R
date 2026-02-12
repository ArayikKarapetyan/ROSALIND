# FUNCTION 1: Count nucleotides in DNA string
count_nucleotides_r <- function(dna_string) {

  nucleotides <- strsplit(dna_string, "")[[1]]
  
  # Count each nucleotide
  a_count <- sum(nucleotides == "A")
  c_count <- sum(nucleotides == "G")
  g_count <- sum(nucleotides == "C")
  t_count <- sum(nucleotides == "T")
  
  return(c(A = a_count, C = c_count, G = g_count, T = t_count))
}

# FUNCTION 2: Alternative concise R version
count_nucleotides_r_concise <- function(dna_string) {

  nucleotides <- strsplit(dna_string, "")[[1]]
  counts <- table(nucleotides)
  
  # Ensure all nucleotides are present (some might be missing)
  result <- c(
    A = ifelse("A" %in% names(counts), counts["A"], 0),
    C = ifelse("C" %in% names(counts), counts["C"], 0),
    G = ifelse("G" %in% names(counts), counts["G"], 0),
    T = ifelse("T" %in% names(counts), counts["T"], 0)
  )
  
  return(as.vector(result))
}

# FUNCTION 3: Another alternative using factor levels
count_nucleotides_r_factors <- function(dna_string) {
  
  nucleotides <- strsplit(dna_string, "")[[1]]
  counts <- table(factor(nucleotides, levels = c("A", "C", "G", "T")))
  return(as.vector(counts))
}

input_file <- "Dataset.txt"  # Make sure the path is correct

# Read the entire file
file_content <- readLines(input_file, warn = FALSE)

# Combine all lines into a single string
data <- paste(file_content, collapse = "")

# Clean up: remove any whitespace or newlines
data <- gsub("\\s", "", data)

# Method 1: Using basic R
result1 <- count_nucleotides_r(data)
cat("Method 1 (Basic R):\n")
cat(paste(result1, collapse = " "), "\n\n")

# # Method 2: Using table()
# result2 <- count_nucleotides_r_concise(data)
# cat("Method 2 (Using table()):\n")
# cat(paste(result2, collapse = " "), "\n\n")

# # Method 3: Using factors
# result3 <- count_nucleotides_r_factors(data)
# cat("Method 3 (Using factors):\n")
# cat(paste(result3, collapse = " "), "\n\n")

