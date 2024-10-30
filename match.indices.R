
# a function to trace indices of the props that will change
match.indices <- function(prop, ind.nochange, props2change){
  
  # for loop to generate indices of the props that will change
  prop.ind <- prop # renaming the initial props for a purpose
  prop.ind[ind.nochange] <- 5 # replacing the props that I don't intend to change, 
  #this will help when matching indices if the props I intend to change (see the immediate below for loop)
  
  
  # a container to store props vector indices that will change and corresponds to small props
  matching_indices <- vector()
  
  for (val in props2change) {
    
    # find index of closest match in prop.ind vector
    index <- which.min(abs(prop.ind - val)) # big absolute differences will never win here, that's why I replaced props not 
    # changing with 5 to give big differences
    
    # check if this index has been used already to avoid having two or more same indices
    if (!index %in% matching_indices) {
      
      matching_indices <- c(matching_indices, index) # if index is not used, assign it
      
    } else {
      
      # if index already used, find next closest match
      remaining_indices <- setdiff(seq_along(prop.ind), matching_indices) # eliminating index already used
      
      next_index <- remaining_indices[which.min(abs(prop.ind[remaining_indices] - val))] # using the next closest index
      
      matching_indices <- c(matching_indices, next_index) # adding the index to my matching_indices vector
    }
  }
  
  return(matching_indices)
}