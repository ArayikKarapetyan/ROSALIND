rna_to_protein_r <- function(rna_string) {
  # RNA codon table
  codon_table <- list(
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
  
  # Split RNA into codons
  rna_len <- nchar(rna_string)
  protein <- character()
  
  for (i in seq(1, rna_len - 2, by = 3)) {
    codon <- substr(rna_string, i, i + 2)
    amino_acid <- codon_table[[codon]]
    
    if (is.null(amino_acid)) {
      break  # Invalid codon
    }
    
    if (amino_acid == "*") {
      break  # Stop codon
    }
    
    protein <- c(protein, amino_acid)
  }
  
  return(paste(protein, collapse = ""))
}

rna_string <- readLines("Dataset.txt", warn = FALSE)
rna_string <- gsub("\\s", "", rna_string)  # Remove whitespace
print(paste("R Basic:", rna_to_protein_r(rna_string)))