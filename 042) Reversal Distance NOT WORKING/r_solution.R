library(collections)

N <- 10
TARGET <- 1:N

# Encode permutation as integer
encode <- function(p) {
  sum(p * 11^(seq_len(N)-1))
}

# Generate neighbors
neighbors <- function(p) {
  res <- list()
  idx <- 1
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      q <- p
      q[i:j] <- rev(q[i:j])
      res[[idx]] <- q
      idx <- idx + 1
    }
  }
  res
}

reversal_distance_bidirectional <- function(start) {
  if (all(start == TARGET)) return(0)

  q1 <- deque()
  q2 <- deque()

  d1 <- new.env(hash = TRUE)
  d2 <- new.env(hash = TRUE)

  e_start <- encode(start)
  e_target <- encode(TARGET)

  q1$push(list(p = start, d = 0))
  q2$push(list(p = TARGET, d = 0))

  d1[[as.character(e_start)]] <- 0
  d2[[as.character(e_target)]] <- 0

  while (q1$size() > 0 && q2$size() > 0) {

    # Expand smaller frontier
    if (q1$size() <= q2$size()) {
      cur <- q1$popleft()
      for (nxt in neighbors(cur$p)) {
        e <- encode(nxt)
        key <- as.character(e)
        if (!is.null(d1[[key]])) next
        if (!is.null(d2[[key]])) {
          return(cur$d + 1 + d2[[key]])
        }
        d1[[key]] <- cur$d + 1
        q1$push(list(p = nxt, d = cur$d + 1))
      }
    } else {
      cur <- q2$popleft()
      for (nxt in neighbors(cur$p)) {
        e <- encode(nxt)
        key <- as.character(e)
        if (!is.null(d2[[key]])) next
        if (!is.null(d1[[key]])) {
          return(cur$d + 1 + d1[[key]])
        }
        d2[[key]] <- cur$d + 1
        q2$push(list(p = nxt, d = cur$d + 1))
      }
    }
  }

  -1
}

# Normalize permutation
normalize <- function(p1, p2) {
  pos <- match(p1, p2)
  pos
}

main <- function() {
  lines <- readLines("Dataset.txt")
  lines <- lines[lines != ""]

  results <- integer()

  for (i in seq(1, length(lines), by = 2)) {
    p1 <- as.integer(strsplit(lines[i], " ")[[1]])
    p2 <- as.integer(strsplit(lines[i+1], " ")[[1]])

    perm <- normalize(p1, p2)
    results <- c(results, reversal_distance_bidirectional(perm))
  }

  cat(paste(results, collapse = " "), "\n")
}

main()
