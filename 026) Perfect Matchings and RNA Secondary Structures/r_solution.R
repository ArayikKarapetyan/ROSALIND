parse_fasta_r <- function(fasta_string) {
  lines <- strsplit(fasta_string, "\n")[[1]]
  sequences <- list()
  current_id <- ""
  
  for (line in lines) {
    if (startsWith(line, ">")) {
      current_id <- substring(line, 2)
      sequences[[current_id]] <- ""
    } else {
      sequences[[current_id]] <- paste0(sequences[[current_id]], line)
    }
  }
  
  return(as.character(sequences)[1])
}

count_perfect_matchings_r <- function(rna) {
  # Count bases
  bases <- strsplit(rna, "")[[1]]
  a_count <- sum(bases == "A")
  u_count <- sum(bases == "U")
  c_count <- sum(bases == "C")
  g_count <- sum(bases == "G")
  
  # Check condition
  if (a_count != u_count || c_count != g_count) {
    return(0)
  }
  
  # Calculate factorials
  au_matchings <- factorial(a_count)
  cg_matchings <- factorial(c_count)
  
  return(au_matchings * cg_matchings)
}

# Using table for counting
count_perfect_matchings_table_r <- function(rna) {
  # Create frequency table
  bases <- strsplit(rna, "")[[1]]
  freq <- table(bases)
  
  # Get counts (handle missing bases by returning 0)
  a_count <- ifelse("A" %in% names(freq), freq["A"], 0)
  u_count <- ifelse("U" %in% names(freq), freq["U"], 0)
  c_count <- ifelse("C" %in% names(freq), freq["C"], 0)
  g_count <- ifelse("G" %in% names(freq), freq["G"], 0)
  
  # Check condition
  if (a_count != u_count || c_count != g_count) {
    return(0)
  }
  
  # Calculate
  return(factorial(a_count) * factorial(c_count))
}

rna <- parse_fasta_r(paste(readLines("26) Perfect Matchings and RNA Secondary Structures\\Dataset.txt"), collapse = "\n"))

cat("Perfect matchings:", count_perfect_matchings_r(rna), "\n")