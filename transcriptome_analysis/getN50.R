# N50 calculation
nFifty <- function(df) {
  L <- c(2, 2, 3, 3, 5, 7) # Store length of each conitg in a vector.
  Lr <- unlist(tapply(L, L, function(x) rep(x[1], sum(x))))
  median(Lr) 
}

# Get N50
nFifty(comboFrame)
