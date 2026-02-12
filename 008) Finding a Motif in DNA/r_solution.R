find_substring_positions_r <- function(s, t) {
  s_chars <- strsplit(s, "")[[1]]
  t_chars <- strsplit(t, "")[[1]]
  t_len <- length(t_chars)
  s_len <- length(s_chars)
  
  positions <- numeric()
  
  for (i in 1:(s_len - t_len + 1)) {
    if (all(s_chars[i:(i + t_len - 1)] == t_chars)) {
      positions <- c(positions, i)
    }
  }
  
  return(positions)
}

# Alternative using gregexpr (doesn't find overlapping matches)
find_substring_positions_gregexpr <- function(s, t) {
  # gregexpr finds non-overlapping matches by default
  matches <- gregexpr(t, s, fixed = TRUE)[[1]]
  
  # Convert to positions (filter out -1 for no match)
  positions <- as.numeric(matches)
  positions <- positions[positions > 0]
  
  return(positions)
}

input_data <- readLines("Dataset.txt", warn = FALSE)
s <- input_data[1]
t <- input_data[2]

print(paste("R Manual:", paste(find_substring_positions_r(s, t), collapse = " ")))