build_trie_r <- function(strings) {
  nodes <- list()
  nodes[[1]] <- list()  # root
  
  edges <- data.frame(
    parent = integer(0),
    child = integer(0),
    char = character(0),
    stringsAsFactors = FALSE
  )
  
  next_id <- 2
  
  for (s in strings) {
    current_id <- 1
    chars <- strsplit(s, "")[[1]]
    
    for (char in chars) {
      if (!is.null(nodes[[current_id]][[char]])) {
        current_id <- nodes[[current_id]][[char]]
      } else {
        nodes[[current_id]][[char]] <- next_id
        nodes[[next_id]] <- list()
        
        edges <- rbind(
          edges,
          data.frame(
            parent = current_id,
            child = next_id,
            char = char,
            stringsAsFactors = FALSE
          )
        )
        
        current_id <- next_id
        next_id <- next_id + 1
      }
    }
  }
  
  edges
}

main <- function() {
  # Read input
  strings <- readLines("53) Introduction to Pattern Matching\\Dataset.txt", warn = FALSE)
  strings <- strings[strings != ""]
  
  edges <- build_trie_r(strings)
  
  # Write to output.txt
  writeLines(
    apply(edges, 1, function(row)
      paste(row["parent"], row["child"], row["char"])
    ),
    "output.txt"
  )
}

# Run
main()
