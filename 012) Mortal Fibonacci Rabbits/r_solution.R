mortal_rabbits_r <- function(n, m) {
  if (n <= 0) return(0)
  if (n == 1) return(1)
  
  # Initialize age vector
  ages <- integer(m)
  ages[1] <- 1  # R is 1-indexed: ages[1] = newborns
  
  for (month in 2:n) {
    # Newborns = sum of rabbits age 2 to m (reproducing ages)
    newborns <- sum(ages[2:m])
    
    # Age rabbits: shift down
    # Rabbits at age m die
    for (age in m:2) {
      ages[age] <- ages[age-1]
    }
    
    # Set newborns
    ages[1] <- newborns
  }
  
  return(sum(ages))
}

input_data <- readLines("Dataset.txt", warn = FALSE)
values <- strsplit(input_data, " ")[[1]]
n <- as.integer(values[1])
m <- as.integer(values[2])

print(paste("R Basic:", mortal_rabbits_r(n, m)))