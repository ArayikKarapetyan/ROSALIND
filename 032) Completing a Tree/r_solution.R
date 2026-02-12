count_edges_to_complete_tree_r <- function(n, edges_list) {
  # Build adjacency list
  adj_list <- vector("list", n)
  for (i in 1:n) {
    adj_list[[i]] <- integer(0)
  }
  
  for (edge in edges_list) {
    u <- edge[1]
    v <- edge[2]
    adj_list[[u]] <- c(adj_list[[u]], v)
    adj_list[[v]] <- c(adj_list[[v]], u)
  }
  
  # Count connected components using DFS
  visited <- rep(FALSE, n)
  components <- 0
  
  for (node in 1:n) {
    if (!visited[node]) {
      components <- components + 1
      stack <- c(node)
      visited[node] <- TRUE
      
      while (length(stack) > 0) {
        current <- stack[1]
        stack <- stack[-1]
        
        for (neighbor in adj_list[[current]]) {
          if (!visited[neighbor]) {
            visited[neighbor] <- TRUE
            stack <- c(neighbor, stack)
          }
        }
      }
    }
  }
  
  # Minimum edges needed to connect components
  return(components - 1)
}

read_input <- function(filename) {
  lines <- readLines(filename, warn = FALSE)
  n <- as.integer(lines[1])
  
  edges_list <- list()
  for (i in 2:length(lines)) {
    if (nchar(lines[i]) > 0) {
      values <- as.integer(strsplit(lines[i], " ")[[1]])
      edges_list[[length(edges_list) + 1]] <- c(values[1], values[2])
    }
  }
  
  return(list(n = n, edges = edges_list))
}

data <- read_input("Dataset.txt")

n <- data$n
edges_list <- data$edges

result <- count_edges_to_complete_tree_r(n, edges_list)
cat(result, "\n")
