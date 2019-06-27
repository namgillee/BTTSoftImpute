get_bttframe <- function(v, dir_rem = 0)
{
  #Some neighboring core tensors are removed from TT tensor v. 
  #If dir_rem = 0,    remove G[[bk]]
  #If dir_rem = 1,    remove G[[bk]], G[[bk+1]]
  #If dir_rem >= 2,   remove G[[bk]], G[[bk+1]], ..., G[[bk+dir_rem]]
  #If dir_rem = -1,   remove G[[bk]], G[[bk-1]]
  #If dir_rem <= -2,  remove G[[bk]], G[[bk-1]], ..., G[[bk+dir_rem]]
  #
  #The removed TT cores are replaced by reshaping of identity matrices.
  #Still, G[[bk]] has the block mode K.
  
  
  num_removed_core = abs(dir_rem) + 1
  
  
  if (v$N < v$block + dir_rem) {
    warning('get_btt_frame(): Cannot remove the required number of TT-cores (to the right)')
    
  } else if (1 > v$block + dir_rem) {
    warning('get_btt_frame(): Cannot remove the required number of TT-cores (to the left)')
  }
    
  
  bk = v$block
  J = v$J
  
  if (dir_rem >= 0) {
    
      R_start = v$R[bk]
      R_end = v$R[bk+dir_rem+1]
      J_all = J[bk:(bk+dir_rem)]
      
      #Remove the first TT-core: which is a 4th-order tensor
      R1 = R_start
      R2 = prod(J_all[-1]) * R_end
      J1 = J_all[1]
      
      crnew = diag(1, R1*J1*R2)
      dim(crnew) = c(R1, J1, R2, R1*J1*R2)
      crnew = aperm(crnew, c(1,2,4,3))
      
      v$G[[bk]] = crnew
      
      v$K = R1*J1*R2
    
      #Remove the rest TT-cores
      if (dir_rem >= 1) {
      
          #Increasing order, until N
          for (cur_core in 2:(dir_rem+1)) {
            v$R[bk+cur_core-1] = R2    #Update R[+1]
            J_all = J_all[-1]
            
            R1 = R2
            R2 = prod(J_all[-1]) * R_end
            J1 = J_all[1]
            
            crnew = diag(1, J1*R2) #R1==J1*R2
            dim(crnew) = c(R1, J1, R2)
        
            v$G[[bk+cur_core-1]] = crnew
            
          }
      }
      
  } else {
      # dir_rem < 0
      
      R_start = v$R[bk+dir_rem]
      R_end = v$R[bk+1]
      J_all = J[(bk+dir_rem):bk]
      
      #Remove the first TT-core: which is a 4th-order tensor
      R1 = prod(J_all[-(-dir_rem+1)]) * R_start
      R2 = R_end
      J1 = J_all[-dir_rem+1]
      
      crnew = diag(1, R1*J1*R2)
      dim(crnew) = c(R1, J1, R2, R1*J1*R2)
      crnew = aperm(crnew, c(1,2,4,3))
      
      v$G[[bk]] = crnew
      
      v$K = R1*J1*R2
      
      #Remove the rest TT-cores
      if (-dir_rem >= 1) {
        
        #Increasing order, until N
        for (cur_core in 1:(-dir_rem)) {
          v$R[bk-cur_core+1] = R1    #Update R[+0]
          J_all = J_all[-length(J_all)]
          
          R2 = R1
          R1 = prod(J_all[-length(J_all)]) * R_start
          J1 = J_all[length(J_all)]
          
          crnew = diag(1, R1*J1) #R2==R1*J1
          dim(crnew) = c(R1, J1, R2)
          
          v$G[[bk-cur_core]] = crnew
          
        }
      }
    
    
  }

  return(v)
}
  