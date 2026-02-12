library(httr)
library(stringr)

fetch_protein_sequence_r <- function(protein_id) {
  # FIX: Strip the ID for the URL (e.g., P01866_GCB_MOUSE -> P01866)
  accession <- strsplit(protein_id, "_")[[1]][1]
  url <- paste0("https://rest.uniprot.org/uniprotkb/", accession, ".fasta")
  
  response <- tryCatch({
    GET(url, timeout(15))
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(response) && status_code(response) == 200) {
    fasta_text <- content(response, "text", encoding = "UTF-8")
    lines <- strsplit(fasta_text, "\n")[[1]]
    
    if (length(lines) > 0 && startsWith(lines[1], ">")) {
      sequence <- paste(lines[-1], collapse = "")
      return(sequence)
    }
  }
  return(NULL)
}

find_n_glycosylation_motif_stringr <- function(sequence) {
  # N-glycosylation motif: N{P}[ST]{P}
  # Pattern: N followed by (not P), then (S or T), then (not P)
  # (?=...) ensures we find overlapping matches like 'NNST'
  pattern <- "(?=(N[^P][ST][^P]))"
  
  # str_locate_all returns a matrix of start/end positions
  matches <- str_locate_all(sequence, regex(pattern))[[1]]
  
  if (nrow(matches) > 0) {
    return(matches[, 1]) # Return 1-based start positions
  } else {
    return(integer(0))
  }
}

# --- Main Execution ---

# Load IDs from the specific Rosalind filename provided
protein_ids <- readLines("Dataset.txt", warn = FALSE) [cite: 2]
protein_ids <- protein_ids[protein_ids != ""]

for (protein_id in protein_ids) {
  sequence <- fetch_protein_sequence_r(protein_id)
  
  if (!is.null(sequence)) {
    positions <- find_n_glycosylation_motif_stringr(sequence)
    
    # Rosalind output: Print ID and then the space-separated positions
    if (length(positions) > 0) {
      cat(protein_id, "\n")
      cat(paste(positions, collapse = " "), "\n")
    }
  }
}