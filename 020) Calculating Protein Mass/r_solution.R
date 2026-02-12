calculate_protein_mass_r <- function(protein) {
  # Monoisotopic mass table
  mass_table <- list(
    A = 71.03711,   # Alanine
    C = 103.00919,  # Cysteine
    D = 115.02694,  # Aspartic acid
    E = 129.04259,  # Glutamic acid
    F = 147.06841,  # Phenylalanine
    G = 57.02146,   # Glycine
    H = 137.05891,  # Histidine
    I = 113.08406,  # Isoleucine
    K = 128.09496,  # Lysine
    L = 113.08406,  # Leucine
    M = 131.04049,  # Methionine
    N = 114.04293,  # Asparagine
    P = 97.05276,   # Proline
    Q = 128.05858,  # Glutamine
    R = 156.10111,  # Arginine
    S = 87.03203,   # Serine
    T = 101.04768,  # Threonine
    V = 99.06841,   # Valine
    W = 186.07931,  # Tryptophan
    Y = 163.06333   # Tyrosine
  )
  
  # Split protein into individual amino acids
  aa_chars <- strsplit(protein, "")[[1]]
  
  total_mass <- 0.0
  
  for (aa in aa_chars) {
    mass <- mass_table[[aa]]
    if (is.null(mass)) {
      stop(paste("Invalid amino acid:", aa))
    }
    total_mass <- total_mass + mass
  }
  
  return(total_mass)
}


protein <- readLines("Dataset.txt", warn = FALSE)
protein <- gsub("\\s", "", protein)  # Remove whitespace

print(paste("R Basic:", round(calculate_protein_mass_r(protein), 3)))
