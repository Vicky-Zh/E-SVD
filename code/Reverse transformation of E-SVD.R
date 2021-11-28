Decompression <- function(theta){
  n = dim(theta)[1]
  l = dim(theta)[2]
  c = cos(theta)
  d = sin(theta)
  #Create the iteration matrix
  A_itr = matrix(0,nrow = n, ncol = l)
  A_itr[1:l, 1:l] = diag(1, l, l)
  B= A_itr
  B_new = B

  for (i in as.numeric(n-l==0):(l-1)){
    print(paste('colomn = ',i))
    h = l-i
    for (j in 0:(i-1+n-l)){
      v = n-j
      for (k in 1:l){
        B_new[h,k] = c[v,h] * B[h,k] - d[v,h] * B[v,k]
        B_new[v,k] = d[v,h] * B[h,k] + c[v,h] * B[v,k]
      }
      B = B_new
    }
  }
  return(B)
}



