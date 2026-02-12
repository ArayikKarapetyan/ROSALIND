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

build_overlap_graph_r <- function(sequences, k = 3) {
  edges <- matrix(character(), ncol = 2, nrow = 0)
  ids <- names(sequences)
  
  for (i in seq_along(ids)) {
    for (j in seq_along(ids)) {
      if (i == j) next  # Skip self-loops
      
      s_id <- ids[i]
      t_id <- ids[j]
      s_seq <- sequences[[s_id]]
      t_seq <- sequences[[t_id]]
      
      # Check if suffix of s matches prefix of t
      if (nchar(s_seq) >= k && nchar(t_seq) >= k) {
        s_suffix <- substr(s_seq, nchar(s_seq) - k + 1, nchar(s_seq))
        t_prefix <- substr(t_seq, 1, k)
        
        if (s_suffix == t_prefix) {
          edges <- rbind(edges, c(s_id, t_id))
        }
      }
    }
  }
  
  return(edges)
}

# More efficient R version using lists
build_overlap_graph_efficient_r <- function(sequences, k = 3) {
  ids <- names(sequences)
  edges <- list()
  
  # Precompute prefixes and suffixes
  prefixes <- list()
  suffixes <- list()
  
  for (id in ids) {
    seq <- sequences[[id]]
    if (nchar(seq) >= k) {
      prefix <- substr(seq, 1, k)
      suffix <- substr(seq, nchar(seq) - k + 1, nchar(seq))
      
      if (!(prefix %in% names(prefixes))) {
        prefixes[[prefix]] <- character()
      }
      prefixes[[prefix]] <- c(prefixes[[prefix]], id)
      
      if (!(suffix %in% names(suffixes))) {
        suffixes[[suffix]] <- character()
      }
      suffixes[[suffix]] <- c(suffixes[[suffix]], id)
    }
  }
  
  # Find matches
  for (suffix in names(suffixes)) {
    if (suffix %in% names(prefixes)) {
      s_ids <- suffixes[[suffix]]
      t_ids <- prefixes[[suffix]]
      
      for (s_id in s_ids) {
        for (t_id in t_ids) {
          if (s_id != t_id) {
            edges <- c(edges, list(c(s_id, t_id)))
          }
        }
      }
    }
  }
  
  return(edges)
}

# Using data.frame for edges
build_overlap_graph_df_r <- function(sequences, k = 3) {
  ids <- names(sequences)
  edges <- data.frame(from = character(), to = character())
  
  for (i in seq_along(ids)) {
    for (j in seq_along(ids)) {
      if (i == j) next
      
      s_id <- ids[i]
      t_id <- ids[j]
      s_seq <- sequences[[s_id]]
      t_seq <- sequences[[t_id]]
      
      if (nchar(s_seq) >= k && nchar(t_seq) >= k) {
        s_suffix <- substr(s_seq, nchar(s_seq) - k + 1, nchar(s_seq))
        t_prefix <- substr(t_seq, 1, k)
        
        if (s_suffix == t_prefix) {
          edges <- rbind(edges, data.frame(from = s_id, to = t_id))
        }
      }
    }
  }
  
  return(edges)
}


fasta_data <- readLines("Dataset.txt", warn = FALSE)
fasta_string <- paste(fasta_data, collapse = "\n")
sequences <- parse_fasta_r(fasta_string)

edges <- build_overlap_graph_r(sequences)

cat("R Basic:\n")
for (i in 1:nrow(edges)) {
  cat(edges[i, 1], edges[i, 2], "\n")
}

edges_eff <- build_overlap_graph_efficient_r(sequences)
cat("\nR Efficient:\n")
for (edge in edges_eff) {
  cat(edge[1], edge[2], "\n")
}

edges_df <- build_overlap_graph_df_r(sequences)
cat("\nR Data Frame:\n")
for (i in 1:nrow(edges_df)) {
  cat(edges_df$from[i], edges_df$to[i], "\n")
}