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
    #For each row below the diagonal element
    if((i+1)<=n){
      j = i+1
      s_temp = A_itr[i,i]
      s_new = sqrt(s_temp^2 + A_itr[i+1,i]^2)
      #Compute cos and sin of the rotation angle
      c_temp = s_temp/s_new
      d_temp = A_itr[j,i]/s_new
      s[j,i] = s_new
      c[j,i] = c_temp
      d[j,i] = d_temp
      s_temp = s_new
      #Figure out the rotation angle theta
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
     #Compute A^{(k+1)} = A_new
      if((i+1)<=l){
        for (k in (i+1):l){
          A_new[i+1,k] = -d_temp * A_itr[i,k] + c_temp * A_itr[i+1,k]
        }
      }
    }
    
    if ((i+2)<=n){
      for (j in (i+2):n){
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










  
  