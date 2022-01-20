Givens_transformation <- function(A){
  #A_itr is used to store matrix after iteration of each coloum.
  n = dim(A)[1]
  l = dim(A)[2]
  A_itr = A
  #Create matrices to store parameters.
  s = matrix(nrow = n, ncol = l)
  c = matrix(nrow = n, ncol = l)
  d = matrix(nrow = n, ncol = l)
  theta = matrix(nrow = n, ncol = l)
  #Givens transformation
  #For each column
  for (i in 1:l){
    A_new = A_itr
    if ((i+1)<=n){
      s_temp = A_itr[i,i]
      for (j in (i+1):n){
        s_new = sqrt(s_temp^2 + A_itr[j,i]^2)
        c_temp = s_temp/s_new
        d_temp = A_itr[j,i]/s_new
        s[j,i] = s_new
        c[j,i] = c_temp
        d[j,i] = d_temp
        
        tan_theta = d_temp/c_temp
        if (tan_theta<0 && d_temp>0){
          theta[j,i] = atan(tan_theta)+pi
        }else{
          if (tan_theta>0 && d_temp<0){
            theta[j,i] = atan(tan_theta)-pi
          }else{
            theta[j,i] = atan(tan_theta)
          }
        }
        if((i+1)<=l){
          for (q in (i+1):l){
            A_new[j,q] = -d_temp / s_temp * sum(A_itr[,i][1:(j-1)] * A_itr[,q][1:(j-1)]) + c_temp * A_itr[j,q]
          }
        }
        s_temp = s_new
      }
    }
    # Assign the i^th row and column.
    A_itr = A_new
    A_itr[i,] = c(rep(0, i-1), 1, rep(0, l-i))
    A_itr[,i] = c(rep(0, i-1), 1, rep(0, n-i))
  }
  return(theta = theta)
}










  
  