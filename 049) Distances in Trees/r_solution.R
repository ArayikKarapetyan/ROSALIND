library(ape)

main <- function() {
  lines <- readLines("Dataset.txt", warn = FALSE)
  lines <- lines[lines != ""]

  res <- character()
  i <- 1

  while (i <= length(lines)) {
    tree_str <- lines[i]
    i <- i + 1

    nodes <- strsplit(lines[i], " ")[[1]]
    a <- nodes[1]
    b <- nodes[2]
    i <- i + 1

    tree <- read.tree(text = tree_str)

    # ape automatically sets branch length = 1
    D <- dist.nodes(tree)

    # tips are 1:N
    ta <- which(tree$tip.label == a)
    tb <- which(tree$tip.label == b)

    res <- c(res, as.character(D[ta, tb]))
  }

  cat(paste(res, collapse = " "))
  cat(collapse = "\n")
}

main()
