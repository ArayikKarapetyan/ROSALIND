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

reverse_complement_r <- function(dna) {
  complement <- list(A = "T", T = "A", C = "G", G = "C")
  dna_chars <- strsplit(dna, "")[[1]]
  rev_chars <- rev(dna_chars)
  comp_chars <- sapply(rev_chars, function(x) complement[[x]])
  return(paste(comp_chars, collapse = ""))
}

translate_codon_r <- function(codon) {
  genetic_code <- list(
    ATA = "I", ATC = "I", ATT = "I", ATG = "M",
    ACA = "T", ACC = "T", ACG = "T", ACT = "T",
    AAC = "N", AAT = "N", AAA = "K", AAG = "K",
    AGC = "S", AGT = "S", AGA = "R", AGG = "R",
    CTA = "L", CTC = "L", CTG = "L", CTT = "L",
    CCA = "P", CCC = "P", CCG = "P", CCT = "P",
    CAC = "H", CAT = "H", CAA = "Q", CAG = "Q",
    CGA = "R", CGC = "R", CGG = "R", CGT = "R",
    GTA = "V", GTC = "V", GTG = "V", GTT = "V",
    GCA = "A", GCC = "A", GCG = "A", GCT = "A",
    GAC = "D", GAT = "D", GAA = "E", GAG = "E",
    GGA = "G", GGC = "G", GGG = "G", GGT = "G",
    TCA = "S", TCC = "S", TCG = "S", TCT = "S",
    TTC = "F", TTT = "F", TTA = "L", TTG = "L",
    TAC = "Y", TAT = "Y", TAA = "*", TAG = "*",
    TGC = "C", TGT = "C", TGA = "*", TGG = "W"
  )
  
  return(genetic_code[[codon]])
}

translate_orf_r <- function(dna_sequence, start_pos) {
  protein <- ""
  i <- start_pos
  
  while (i + 3 <= nchar(dna_sequence)) {
    codon <- substr(dna_sequence, i, i + 2)
    amino_acid <- translate_codon_r(codon)
    
    if (is.null(amino_acid)) {
      return("")  # Invalid codon
    }
    
    if (amino_acid == "*") {
      return(protein)  # Stop codon
    }
    
    protein <- paste0(protein, amino_acid)
    i <- i + 3
  }
  
  return("")  # No stop codon found
}

find_all_orfs_r <- function(dna) {
  proteins <- character(0)
  
  # Check forward and reverse complement
  sequences <- c(dna, reverse_complement_r(dna))
  
  for (seq in sequences) {
    # Check 3 reading frames
    for (frame in 0:2) {
      i <- frame + 1  # R is 1-indexed
      seq_len <- nchar(seq)
      
      while (i + 2 <= seq_len) {
        codon <- substr(seq, i, i + 2)
        
        if (codon == "ATG") {
          protein <- translate_orf_r(seq, i)
          
          if (nchar(protein) > 0) {
            proteins <- c(proteins, protein)
          }
          
          i <- i + 3  # Move to next codon
        } else {
          i <- i + 1  # Check next position
        }
      }
    }
  }
  
  return(unique(proteins))
}


parse_fasta_r <- function(fasta_string) {
  lines <- strsplit(fasta_string, "\n")[[1]]
  sequence <- paste(lines[!startsWith(lines, ">")], collapse = "")
  return(sequence)
}

fasta_data <- readLines("Dataset.txt", warn = FALSE)
dna <- parse_fasta_r(paste(fasta_data, collapse = "\n"))

proteins <- find_all_orfs_r(dna[[1]])

for (protein in proteins) {
  cat(protein, "\n")
}