s1 <- "TCCCAGTTATGTCAGGGGACACGAGCATGCAGAGAC"   # sequence v in book (vertical sequence in Figure 6.16)
s2 <- "AATTGCCGCCGTCGTTTTCAGCAGTTATGTCAGATC"    # sequence w in book (horizontal sequence in Figure 6.16)

s5 <- "TCTCATCCACCACTCCTAATACTGTGTGGACCGACCTAACCCTGGGCGTGTGCGTCTTTTCTGGGGTCACCGAACGCGGTGGACCGTACCTTTCCCCATAAATCAGCACATATTCACGTGACTTCTCTCTCCCAAAGCGATCCATTATACGACGGTTATAAGTGCCGTTTCTGAGTACCTTAAAACCGTCGCAATTTATACAGGGCGCATCTACCTTCAGTTGGAAACTCGCTCTCCAACGTACCAGCCGACCGAGGCGCCAACACTGCACGCAACCGGTAGCCGAGTGGTAGTCTGTGTCTTTCATGGAATTCCCACTTGAAGAGCAGGGACCGCCACATTATCCCAACAAGTAAGTTCTGACGGCTAACCCTAGTTCCCGTACGGGCCACTCCGGTTCCCCTAGTACTTTTGAATTAAAGTGTGCATCGTTACGTCTGAAACCAGGGTTGGCTCTTTTTTCACGGTGTCTTTCGCAGCGCACTTGAGGCAGGCCCGCGCGAGGGAATACCGGGGCACGCCCTTGTAACCCCTGGCTGTTGCTAATTTCTACTCTCCTGTCTGATTAAGCAGGTTTCGGTTCGATTCAATCCAATAGTTTAGTTCTTGATGCCATCCGTCCTTATGTGTTTCTAGGTAGATCCGTATCACTTTGTCTGTTGGTTTTGACTGGCTCCCGTCGCCTTTGAACGAACCTTGTAGCCCCCCGATCTTGGGGAAGCTATTAAAAGTACTCAACCAGACTAAGCACGGGCAGCTGGGAATATGATTGCCGTCTAGGTGGCGGTTCGTGATCGTTCTCGGCCGCGCCGTTAGCTTACGTAG"    #na20766a haplotype (ACE2)
s6 <- "TCTCATCCACCACTCCTAATACTGTGTGGACCGACCTAACCCTGGGCGTGTGCGCCTTTTTTGGGATCACCGAACGCGGTGGACCGTACCTTTCCCCATAAATCAGCACATATTCACGTGACTTCAAAGCGATCCATTGTATGACGGTTATAAGTGCCGTTTCTGAGTAGCTTAAAACCGATCGCAATTTATACAGGGCGCATCTACCTTCAGTTGGAAACTCGGTCTCCAACGTACCAGCCGACCGAGGCGCCAACACTGCACGCAACCGGTAGCCGAGTGGTAGTCTGTGTCTTTCATGGAATTCCCACTTGAAGAGCAGGGACCGCCACATAAGATTATCCCAACAAGTATGTTCTGACGGCTAACCTTGGTTCCCGTACGGGCCACTCCGGTTCCCCTAGTACTTTTAATTGACTTAAAGTGTGTAACGTTACGTCTTAAACCAGGGCTGGCTCTTTTTTCACGGTGTCTTTCGCAGCGCACTTGAGGCAGGCCCGCGCGAGAGAATACCGGGGCACGCCCTTGTAACCCCTGGCTGTTGCCAATTTCTACTCTTCTGGCTGATTAAGTAGGTTTCGGTTCGATTCAATCCAATAGTTTAGTTCTTGATGCCAAAGGAAGGTCCGTCCTTATGTGTTTCTAGGTAGTCCGTATCACTTTGTCTGTTAGTTTTGACTGGTCTCCCGTCGCCTTTGAACGGACACTTGTAGCCCCACCGATCTTGGGGAAGCTATTAAAAGCACTCAACCAGACTAAGCCGAGCAGCTGGGAACATGATTGCCGTCTAGGTGGCGGTTCGTGATCGTTCTCGGCCGCGCCGTTAGCTTACGTAG" #hg03367a

find_lsa <- function(v, w, sigma, mu)  
{
  v <- strsplit(v, "")    
  w <- strsplit(w, "")  
  v <- as.vector(v[[1]]) 
  w <- as.vector(w[[1]])
  n <- length(v)  # length of first sequence
  m <- length(w)   # length of second sequence
  print(paste("length of sequence v: ", n, sep = ""), quote = F)
  print(paste("length of sequence w: ", m, sep = ""),quote = F)
  s <- matrix(0, nrow=n+1, ncol=m+1)  # s is the dynamic programming table
  b <- matrix(0, nrow=n+1, ncol=m+1)  # b is the matrix of backtracing pointers
  # BACKTRACING POINTERS: 1 to represent up arrow, 2 for left arrow, and 3 for diagonal arrow
  for (i in 2:(n+1)){
    for (j in 2:(m+1)) {
      # dynamic programming table, s
      s[i,j] <- max(0, s[i-1,j]-sigma, s[i,j-1]-sigma) 
      if (v[i-1] == w[j-1]) {
        s[i,j] <- max(s[i,j], (s[i-1, j-1] + 1))
      } else {
        s[i,j] <- max(s[i,j], s[i-1,j-1] - mu)
      }
      # backtracing pointer matrix, b
      if(s[i,j] == s[i-1,j] - sigma) 
        b[i,j] <- 1
      else if(s[i,j] == s[i,j-1] - sigma)
        b[i,j] <- 2
      else { 
        if(v[i-1] == w[j-1]) {
          if (s[i,j] == s[i-1,j-1] + 1)
            b[i,j] <- 3
        } else {
          if (s[i,j] == s[i-1,j-1] - mu)
            b[i,j] <- 3
        }  
      }  
    }
  }
  return(list("s"=s,"b"=b, "vectorV"=v, "vectorW"=w))   
}

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
      if(s[i,j] == s[i-1,j] - sigma) 
        b[i,j] <- 1
      else if(s[i,j] == s[i,j-1] - sigma)
        b[i,j] <- 2
      else { 
        if(v[i-1] == w[j-1]) {
          if (s[i,j] == s[i-1,j-1] + 1)
            b[i,j] <- 3
        } else {
          if (s[i,j] == s[i-1,j-1] - mu)
            b[i,j] <- 3
        }  
      }  
    }
  }
  return(list("s"=s,"b"=b, "vectorV"=v, "vectorW"=w))   
}

plot_edit_graph <- function(dst, backtrace)  
  # PARAMETERS
  # pass <object>$s -- returned by find_lsa() -- as the argument to the *dst* parameter
  # pass <object>$b -- returned by find_lsa() -- as the argument to the *backtrace* parameter, 
  #  <object> is the name you gave to the list returned by find_lsa() -- e.g., q$s and q$b if you ran q <- find_lsa(...)
{
  len <- dim(backtrace)
  v <- vector(mode="numeric")
  w <- vector(mode="numeric")
  for(i in 2:len[1]) { # rows
    for (j in 2:len[2]) { #cols
      if (backtrace[i,j] == 3) {
        v <- c(v,i)
        w <- c(w,j)
      }
    }
  }
  plot(1, type="n",xlim=range(1:len[2]), ylim=rev(range(1:len[1])), xlab="sequence w", ylab="sequence v", font.lab=2) 
  text(w, v, "\\#H2298" , vfont=c("serif", "plain"), xpd=TRUE )
  r <- which(dst==max(dst), arr.ind=T)
  x <- length(r[,1]) # number of equally maximal cells
  for (i in 1:x)
    points(r[i,2], r[i,1], cex = 4, col = "magenta")
}

print_local_alignment <- function(q)  #Pass list returned by find_lsa() as argument to parameter q
{
  dst <- q$s
  backtrace <- q$b
  v <- q$vectorV
  w <- q$vectorW
  terminate <- which(dst==max(dst), arr.ind=T) # coordinates of max entry in dynamic sequencing table
#  terminate <- terminate[1,]
  begin <- terminate
  while (1) {
    if(backtrace[begin] == 3) {
      begin <- begin - 1
    } else {
      begin <- begin + 1
      break
    } 
  }
  a <- seq(begin[1],terminate[1],1)
  b <- seq(begin[2],terminate[2],1)
  for(i in 1:length(a)) {
    print(paste(v[a[i]-1], w[b[i]-1], sep = " "))
  }
  print(paste("Best local alignment on diagonal between ", begin[1], ", " , begin[2], " and ", terminate[1], ", ", terminate[2], " of edit graph.", sep = "" ), quote=  F)
}

gather_gsa <- function(b,v,w,i,j,fn)  #filename to output temporary alignment; note that  
{
  if (i==1 | j==1)  # stop if current cell is in first row or column
    return()
  if (b[i,j]==3) { #diagonal
    gather_gsa(b,v,w,i-1,j-1,fn)
    line = paste(v[i-1], w[j-1], sep = "\t") 
    write(line, fn, append=T)
  } else if (b[i,j]==1) { #up
    gather_gsa(b,v,w,i-1,j, fn)
    line = paste(v[i-1], "-", sep = "\t")
    write(line, fn, append=T)
  } else if (b[i,j]==2) { #left 
    gather_gsa(b,v,w,i,j-1,fn)
    line =  paste("-", w[j-1], sep="\t")
    write(line, fn, append=T)
  }
}

print_gsa <- function(b,v,w,i,j,fn)
{
  gather_gsa(b,v,w,i,j,fn)
  d <- read.table(file = fn, header = F)
  seq1 <- paste(d[,1], collapse=" ")
  seq2 <- paste(d[,2], collapse=" ")
  line <- ">v"
  write(line, fn, append=F)
  write(seq1, fn, append=T)
  line <- ">w"
  write(line, fn, append=T)
  write(seq2, fn, append = T)
}

