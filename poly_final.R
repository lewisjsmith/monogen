library(dplyr)
pind<-function(x, p){
  
  if(p == 0 || p == 1){
    return(NULL)
  }
  
  # Test Matrix

  # The polynomial power to which you wish to generate
  # p = 6
  
  # Column count
  n = dim(x)[2]
  
  # Currently overestimates column size
  tcol <- 0
  
  # for(pc in c(1:p)){
  #   tcol <- tcol + n^pc
  # }
  
  # For each polynomial set
  tcol <- n^p
  
  # Indices are saved here
  indices <- matrix(0, nrow = 1, ncol = tcol)
  
  # To save indices
  it <- 1
  
  # list of 1, 4, 5, which is found from g - 1 subgroups starting index
  # These are the starting points of the known repeating pattern
  st_ind <- matrix(0, nrow=1, ncol=n)
  st_ind_c <- 1
  gc <- 0
  for(g in c(2:n)){
    for(h in c(1:(g-1))){
      st_ind[st_ind_c] <- h + gc
      st_ind_c <- st_ind_c + 1
    }
    gc <- gc + n
  }
  
  # Part 1
  # Original group size, g
  # Stacking group size, G
  # polynomial power, p
  # Following the first input, groups are of increasing size:
  # They are groups of a size of the original group squared, G
  # The amount of these groups then increases by the product of the original group...
  # for each successive polynomial power
  # e.g. the positions within these G groups are constant, 
  # for g = 3, positions 4, 7, 8 within G = 9 are non-unique, polynomial power = 2
  for(si in c(1:((n^p)/(n^2)))){
    for(sta in st_ind){
      indices[it] <- (n^p)/(n^(p-1)) + (sta + ((n^2) * (si-1)))
      it <- it + 1
    }
  }

  # Part 2 - TO CHANGE THE LIST OF
  # The positions within G become groups of size g each upon the next polynomial power
  # e.g. a previous highlighted y*x, yx, would become a xy group (xy*x, xy*y, xy*z)
  # there currently exists overlap within these two "Parts"
  for(sj in c(0:(n-1))){
    for(sta in pind(x, p-1)){
      if(sta > 0){
        indices[it] <- (((sta - 1))*n) + 1 + sj
        it <- it + 1
      }
    }
  }
  
  #print(which(indices==0))
  
  # Order
  indices <- sort(indices)
  
  # Need to delete duplicates
  
  indices <- indices[!duplicated(indices)]
  
  if(length(which(indices==0))==0){
    return(indices)
  }
  else{
    return(indices[-1])
  }
}
pgen<-function(x, p){

  n = dim(x)[2]
  poly = p - 1
  
  # fin = matrix(x, ncol=dim(x)[2])
  fin <- x
  
  last = 0
  r = 1
  
  rmv_fin <- matrix(0)
  itmat <- matrix(0)
  itmatstart <- 0
  rmv_list <- matrix()
  prev <- 0
  
  nrtotal = 0
  counter = 0
  
  # loop ends when the desired power is reached
  while(r != poly + 1){
    
    prev <- prev+(n^r)
    
    # List to remove
    rmv <- pind(x, r+1)
    
    #print(rmv)
    rmv_list <- cbind(rmv_list, matrix(rmv + prev, ncol=length(rmv)))
    
    it <- 1
    
    # iterate through the columns
    for(i in c(1:dim(x)[2])){
      
      # iterate through the latest products
      for(j in c((1+last):((n^r)+last))){
        
        # if the iteration number is not in hte list, proceed, 
        # the iteration must also be less than the number of permutations
        if(it %in% rmv == FALSE && it <= n^(r+1)){
          fin <- cbind(fin, x[,i]*fin[,j]) # the real code
          it <- it + 1
        }
        
        else{
          fin <- cbind(fin, matrix(NA, nrow = dim(x)[1], ncol=1))
          # itmat <- cbind(itmat, (it+last))
          it <- it + 1
        }
      }
    }
    last = last + n^r
    # itmatstart <- itmatstart + (it - 1)
    # rmv_fin <- cbind(rmv_fin, matrix((rmv + last), nrow = 1))
    r = r + 1
  }
  
  rmv_list <- rmv_list[-1]
  
  # itmat <- itmat[-1]
  # itmat <- itmat[!duplicated(itmat)]
  
  # rmv_fin <- rmv_fin[,-1]
  # rmv_fin <- rmv_fin[!duplicated(rmv_fin)]
  # print(rmv_fin)
  # print(rmv_fin)
  # print(itmat)
  # print(length(itmat))
  #print(counter)
  
  # 
  fin <- subset.matrix(fin, drop=TRUE, select=-rmv_list)
  return(fin)
}