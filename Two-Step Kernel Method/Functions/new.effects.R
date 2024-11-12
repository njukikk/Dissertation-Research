
# a function that will generate new effects for any given prop vector
new.effects <- function(props2change, maxeffect, matching_indices, prop.new){
  
      n <- length(props2change) # number of props to use
      addeff <- rep(NA, n) # a vector to store new effects to be added to props
      
      U_i <- min(maxeffect, 1-props2change) # minimum upper bound for the new effects to be generated
      L_i <- max(-maxeffect, -props2change) # maximum lower bound for the new effects to be generated
      
      # generating random uniform numbers between -1 and 1
      r_i <- runif(n-1, -1, 1)
      
      r.k_2 <- -sum(r_i) # this ensures sum(r_i)=0
      
      r_i <- c(r_i, r.k_2) # the complete vector of r_i's
      
      r.min <- min(r_i) # minimum r_i
      r.max <- max(r_i) # maximum r_i
      
      # find a large enough c such that L_i <= e_i=r_i/c <= U_i
      c1 <- r.min/L_i
      c2 <- r.max/U_i
      c <- max(c1, c2) # c is at least max(c1, c2)
      
      # computing new effects
      neweff <- r_i/c # divide r_i's by a constant c
      
      # finding index whose absolute value corresponds to sub.maxeffect
      ind.with.maxeff <- which(round(abs(neweff),9)==round(maxeffect,9))
      
      # just in case I end up with two |maxeffect| values, this code will pick one, it doesn't matter which one. 
      # Here I picked the minimum index
      if (length(ind.with.maxeff) > 1) {ind.with.maxeff <- min(ind.with.maxeff)}
      
      # assigning this sub.maxeffect to the maximum prop among the small props we started with
      addeff[which.max(props2change)] <- neweff[ind.with.maxeff]
      
      # new effects generated exclusive of the one that corresponded to sub.maxeffect
      neweff <- neweff[-ind.with.maxeff]
      
      # tracing the addeff that are not NAs
      adj <- is.na(addeff)
      
      # all new effects generated for small props
      addeff[adj] <- neweff
      
      # adding the generated new effects to the small proportions we started with
      new.props <- props2change + addeff
      
      # replacing the changed props in the initial prop vector with new props
      prop.new[matching_indices] <- new.props
  
  return(prop.new)
}

