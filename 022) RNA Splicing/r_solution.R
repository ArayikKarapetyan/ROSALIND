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
  
  return(sequences)
}

remove_introns_r <- function(dna, introns) {
  # Remove all introns from DNA
  result <- dna
  
  for (intron in introns) {
    result <- gsub(intron, "", result, fixed = TRUE)
  }
  
  return(result)
}

dna_to_rna_r <- function(dna) {
  # Transcribe DNA to RNA (T -> U)
  return(gsub("T", "U", dna))
}

translate_rna_to_protein_r <- function(rna) {
  # Genetic code translation table
  genetic_code <- list(
    UUU = "F", UUC = "F", UUA = "L", UUG = "L",
    CUU = "L", CUC = "L", CUA = "L", CUG = "L",
    AUU = "I", AUC = "I", AUA = "I", AUG = "M",
    GUU = "V", GUC = "V", GUA = "V", GUG = "V",
    UCU = "S", UCC = "S", UCA = "S", UCG = "S",
    CCU = "P", CCC = "P", CCA = "P", CCG = "P",
    ACU = "T", ACC = "T", ACA = "T", ACG = "T",
    GCU = "A", GCC = "A", GCA = "A", GCG = "A",
    UAU = "Y", UAC = "Y", UAA = "*", UAG = "*",
    CAU = "H", CAC = "H", CAA = "Q", CAG = "Q",
    AAU = "N", AAC = "N", AAA = "K", AAG = "K",
    GAU = "D", GAC = "D", GAA = "E", GAG = "E",
    UGU = "C", UGC = "C", UGA = "*", UGG = "W",
    CGU = "R", CGC = "R", CGA = "R", CGG = "R",
    AGU = "S", AGC = "S", AGA = "R", AGG = "R",
    GGU = "G", GGC = "G", GGA = "G", GGG = "G"
  )
  
  protein <- ""
  rna_len <- nchar(rna)
  
  for (i in seq(1, rna_len - 2, by = 3)) {
    codon <- substr(rna, i, i + 2)
    amino_acid <- genetic_code[[codon]]
    
    if (is.null(amino_acid)) {
      break  # Invalid codon
    }
    
    if (amino_acid == "*") {
      break  # Stop codon
    }
    
    protein <- paste0(protein, amino_acid)
  }
  
  return(protein)
}

rna_splicing_r <- function(fasta_string) {
  sequences <- parse_fasta_r(fasta_string)
  
  # First sequence is DNA, rest are introns
  seq_ids <- names(sequences)
  dna <- sequences[[seq_ids[1]]]
  introns <- sequences[seq_ids[-1]]
  
  # Step 1: Remove introns
  exons_only <- remove_introns_r(dna, introns)
  
  # Step 2: Transcribe to RNA
  rna <- dna_to_rna_r(exons_only)
  
  # Step 3: Translate to protein
  protein <- translate_rna_to_protein_r(rna)
  
  return(protein)
}


fasta_data <- readLines("Dataset.txt", warn = FALSE)
protein <- rna_splicing_r(paste(fasta_data, collapse = "\n"))

cat("R Solution:\n")
cat(protein)