count_rna_strings_r <- function(protein) {
  # RNA codon counts for each amino acid
  codon_counts <- list(
    A = 4,  # Alanine: GCU, GCC, GCA, GCG
    C = 2,  # Cysteine: UGU, UGC
    D = 2,  # Aspartic acid: GAU, GAC
    E = 2,  # Glutamic acid: GAA, GAG
    F = 2,  # Phenylalanine: UUU, UUC
    G = 4,  # Glycine: GGU, GGC, GGA, GGG
    H = 2,  # Histidine: CAU, CAC
    I = 3,  # Isoleucine: AUU, AUC, AUA
    K = 2,  # Lysine: AAA, AAG
    L = 6,  # Leucine: UUA, UUG, CUU, CUC, CUA, CUG
    M = 1,  # Methionine: AUG
    N = 2,  # Asparagine: AAU, AAC
    P = 4,  # Proline: CCU, CCC, CCA, CCG
    Q = 2,  # Glutamine: CAA, CAG
    R = 6,  # Arginine: CGU, CGC, CGA, CGG, AGA, AGG
    S = 6,  # Serine: UCU, UCC, UCA, UCG, AGU, AGC
    T = 4,  # Threonine: ACU, ACC, ACA, ACG
    V = 4,  # Valine: GUU, GUC, GUA, GUG
    W = 1,  # Tryptophan: UGG
    Y = 2   # Tyrosine: UAU, UAC
  )
  
  mod <- 1000000
  total <- 1
  
  # Split protein into individual amino acids
  aa_chars <- strsplit(protein, "")[[1]]
  
  for (aa in aa_chars) {
    count <- codon_counts[[aa]]
    if (is.null(count)) {
      return(0)  # Invalid amino acid
    }
    total <- (total * count) %% mod
  }
  
  # Multiply by stop codon possibilities (3 stop codons)
  total <- (total * 3) %% mod
  
  return(total)
}


protein <- readLines("Dataset.txt", warn = FALSE)
protein <- gsub("\\s", "", protein)  # Remove whitespace

print(paste("R Basic:", count_rna_strings_r(protein)))
