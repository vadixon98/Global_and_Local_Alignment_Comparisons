find_gsa <- function(v, w, sigma, mu)  
{
  v <- strsplit(v, "")    
  w <- strsplit(w, "")  
  v <- as.vector(v[[1]]) 
  w <- as.vector(w[[1]])
  n <- length(v)  # length of first sequence
  m <- length(w)   # length of second sequence
  print (n)
  print(m)
  s <- matrix(0, nrow=n+1, ncol=m+1)  # s is the dynamic programming table
  b <- matrix(0, nrow=n+1, ncol=m+1)  # b is the matrix of backtracing pointers
  # BACKTRACING POINTERS: 1 to represent up arrow, 2 for left arrow, and 3 for diagonal arrow
  for (i in 2:(n+1)){
    for (j in 2:(m+1)) {
      # dynamic programming table, s
      s[i,j] <- max(s[i-1,j]-sigma, s[i,j-1]-sigma)
      if (v[i-1] == w[j-1]) {
        s[i,j] <- max(s[i,j], (s[i-1, j-1] + 1))
      } else {
        s[i,j] <- max(s[i,j], s[i-1,j-1] - mu)
      }
      # backtracing pointer matrix, b
      if (v[i - 1] == w[j - 1]) {
        # Match: diagonal is s[i-1,j-1], up is s[i-1,j] - sigma, left is s[i,j-1] - sigma
        if (s[i,j] == s[i - 1, j - 1]) {
          b[i,j] <- 3  # diagonal
        } else if (s[i,j] == s[i - 1, j] - sigma) {
          b[i,j] <- 1  # up
        } else {
          b[i,j] <- 2  # left
        }
      } else {
        # Mismatch: diagonal is s[i-1,j-1] - mu, up is s[i-1,j] - sigma, left is s[i,j-1] - sigma
        if (s[i,j] == s[i - 1, j - 1] - mu) {
          b[i,j] <- 3  # diagonal
        } else if (s[i,j] == s[i - 1, j] - sigma) {
          b[i,j] <- 1  # up
        } else {
          b[i,j] <- 2  # left
        }
      }
    }
  }
  return(list("s"=s,"b"=b, "vectorV"=v, "vectorW"=w))   
}

print_alignment <- function(b,v,w,i,j) # Prints aligned sequence v in left column and aligned sequence w in right column
{
  if (i==1 | j==1)  # stop if current cell is in first row or column
    return()
  if (b[i,j]==3) { #diagonal
    print_alignment(b,v,w,i-1,j-1)
    print(paste(v[i-1], w[j-1], sep = " "))
  } else if (b[i,j]==1) { #up
    print_alignment(b,v,w,i-1,j)
    print(paste(v[i-1], "-", sep = " "))
  } else if (b[i,j]==2) { #left 
    print_alignment(b,v,w,i,j-1)
    print(paste("-", w[j-1], sep=" "))
  }
}


