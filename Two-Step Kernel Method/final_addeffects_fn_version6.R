
# A function to generate new effects

add.new.effects <- function(maxeffect, prop, n.nochange = 0, ind.nochange = NULL){
  
  
  # maxeffect: maximum effect size in the prop scale
  # prop: proportions of components
  # n.nochange: how many celltypes abundance will not change
  # ind.nochange: indices of the nochange cell types, a vector of length n.nochange
  
  
  # calling the function that will generate new effects for any given proportions
  source("D:/A.School/OSU/Class/1. Dissertation/6. Fall 2024/Simulations/Two-Step Kernel Method/new.effects.R")
  source("D:/A.School/OSU/Class/1. Dissertation/6. Fall 2024/Simulations/Two-Step Kernel Method/match.indices.R")
  
  # a vector that will be updated with the new props generated. Props that won't change will remain the same.
  prop.new <- prop
  
      
  # extracting the props I intend to change from the original props
  if (n.nochange == 0) {prop2use <- prop} else {prop2use <- prop[-(ind.nochange)]}
  
    
  # extracting the indices of the props that are less than or equal to the predefined maxeffect
  small.props_ind <- which(prop2use <= maxeffect)
  
      
  # extracting the indices of the props that are greater than the predefined maxeffect
  big.props_ind <- which(prop2use > maxeffect)
  
  ##############################################################################################################
  ##############################################################################################################
  
      # if else algorithm that will generate effects sizes for the small and big props  
  
      if (length(small.props_ind)>=1){ # checking if I have at least one small prop, if yes, I will generate effects sizes
        # for the two groups (small and large props), otherwise I will generate effect sizes for only the big props
        
        
        # extracting the props less than or equal to predefined maxeffect
        small.props <- prop2use[small.props_ind]
        
        
        # extracting the props greater than predefined maxeffect
        big.props <- prop2use[-(small.props_ind)] 
        
        
        # just in case I end up with only one small prop, this line of code will add one more prop (minimum prop) from the big props.
        if (length(small.props_ind)==1){small.props <- c(small.props, min(big.props))} else {small.props <- small.props}
        
        
        # after adding the minimum big props to the small props, I will remove it from the big props. this ONLY happens if I have one small prop
        if (length(small.props_ind)==1){big.props <- big.props[-(which.min(big.props))]} else {big.props <- big.props}
        
        #  indices of the props that will change and corresponds to small props    
        matching_indices1 <- match.indices(prop, ind.nochange, small.props)
        
        
        # indices of the props that will change and corresponds to big props
        prop.ind <- prop 
        prop.ind[ind.nochange] <- 5 # replacing props not changing with 5
        prop.ind[matching_indices1] <- 5 # replacing props corresponding to small props with
        matching_indices2 <- match.indices(prop.ind, ind.nochange, big.props)
        
        
        # defining a new maxeffect for props<=maxeffect
        sub.maxeffect <- min(small.props)*0.99
        
        # generating new effects for the small props    
        prop.new <- new.effects(small.props, sub.maxeffect, matching_indices1, prop.new)
        
        # generating new effects for the big props
        new.props <- new.effects(big.props, maxeffect, matching_indices2, prop.new)
        
        
    } else {
    
        
        # extracting the props greater than predefined maxeffect
        big.props <- prop2use
                                
        #  indices of the props that will change and corresponds to big props 
        matching_indices2 <- match.indices(prop, ind.nochange, big.props)  
        
        
        # generating new effects for big props of there were no small props  
        new.props <- new.effects(big.props, maxeffect, matching_indices2, prop.new)
         
        } # end of the else statement 
      ##############################################################################################################
  
  return(new.props) # return the new proportions
  
} # end of the function


