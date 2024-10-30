# Load the necessary packages
library("MASS")
library("coda.base")
library("lme4")
library("vegan")
library("MiRKAT")


# Function to compute pairwise Aitchison/Bray-Curtis/Euclidean etc distances where both compositions are nonzero
pairwise_nonzero_dist <- function(propmat, compdist2="aitchison", closure="No") {
  
  n <- nrow(propmat)
  dist.mat <- matrix(NA, n, n) # container to store the distance results
  
  # looping over different pairs of samples
  for (i in 1:n) {
    
    for (j in 1:n) {
      
      # Finding indices where both samples (i & j) have nonzero compositions
      nonzero.indices <- which(propmat[i, ] > 0 & propmat[j, ] > 0)
        
        # Sub-compositions for rows i and j
        sub.comp.i <- propmat[i, nonzero.indices]
        sub.comp.j <- propmat[j, nonzero.indices]
        
        # Combine the two samples with nonzero compositions
        sub.comps <- rbind(sub.comp.i, sub.comp.j)
        
        # Deciding if to normalize the subcompositions or not
        if (closure=="No"){
          
          sub.comps <- sub.comps #maintain the two subcompositions as they are
          
        } else {
          
          sub.comps <- sub.comps/rowSums(sub.comps) # take closure
          
        }
        
        # Compute Aitchison/Bray-Curtis/Euclidean etc distance for the sub-compositions
        if (compdist2 == "BCD"){
          
          dist.mat[i, j] <- dist.mat[j, i] <- vegdist(sub.comps, method="bray")
          
        } else if (compdist2=="canberra"){ 
          
          dist.mat[i, j] <- dist.mat[j, i] <- vegdist(sub.comps, method="canberra")
          
        } else {
          
          dist.mat[i, j] <- dist.mat[j, i] <- vegdist(sub.comps, method=compdist2)
          
        }
    }
  }
  
  return(as.dist(dist.mat))
}


#### Function for the two-step KDC test
dcor.comp.2step <- function(propmat, xvec, num.perm = 1e4, Xdist = "categorical", closure = "No",
                      compdist1 = "hamming", compdist2 = "aitchison", combinethru="sum", 
                      gaussian = F, result = "test"){
  
  # propmat: n x num.cat matrix of proportions (n = sample size, num.cat = number of cell types)
  # xvec: an n x 1 vector containing the values of the predictor
  # num.perm: number of permutations
  # Xdist: specifying the kernel used for the predictor. 
  #       "categorical" - Hamming distance kernel, "LK" - Linear kernel, "gaussian" - Gaussian kernel
  # compdist1: specifying the kernel used for the presence absence compositions in step one. 
  #       "hamming" - Hamming distance, "jaccard" - Jaccard index distance,
  # compdist2: specifying the kernel used for the compositions that are both nonzero in a pair of samples. 
  #       "aitchison" - Aitchison distance,"euclidean" - Euclidean distance, "BCD" - Bray Curtis distance
  # gaussian: T/F, whether to put the distance into a Gaussian Kernel (if FALSE, use D2K like mirkat)
  #         Don't use gaussian = T for compdist2 = "BCD"
  # result: a vector, "test" - provides test results, "cor" - provides value of dcor

  
  
# STEP ONE: Use kernel based on Hamming distance for presence/absence

  # Values > 0 assign 1 else assign 0.
  yvec <- matrix(NA, nrow = nrow(propmat), ncol = ncol(propmat))
  
  for (k in 1:length(propmat)){
    if (propmat[k] > 0) {yvec[k] <- 1} else {yvec[k] <- 0} 
  }
  
  #sum(yvec==0) 
  
  # Creating a Kernel Matrix for step 1
  
  if (compdist1=="hamming"){
    
    # Hamming distance function
    hamming_dist <- function(x, y) {
      sum(x != y)
    }
    
    n <- nrow(yvec) # yvec contains 0's and 1's
    
    ham.dist <- matrix(0, n, n) # distance matrix for step one
    
    # Computing kernel matrix using hamming distance
    for (i in 1:n) {
      
      for (j in 1:n) {
        
        # Compute Hamming distance
        #ham.dist <- hamming_dist(yvec[i, ], yvec[j, ])
        ham.dist[i, j] <- hamming_dist(yvec[i, ], yvec[j, ])
        
        # Step One Kernel: Based on Hamming distance since the data is binary.
        #K11[i, j] <- exp(-ham.dist) # exp(-|X^i - X^j|) where | | is hamming distance
        }
    }
    
    # Step One Kernel: Based on Hamming distance since the data is binary.
    K1 <- D2K(as.matrix(ham.dist))
    
  } else {
    
    # Computing kernel matrix using jaccard distance
    jaccard_dist <- vegdist(yvec, method = "jaccard")
    
    # Step One Kernel: Based on Jaccard distance since the data is binary.
    K1 <- D2K(as.matrix(jaccard_dist))
    
  }

    
  
  # STEP TWO: Defining a Kernel for compositions that are both nonzero in any pair of samples
  
  # Create the Kernel matrix (K2)
  d <- pairwise_nonzero_dist(propmat, compdist2=compdist2, closure=closure) 
  
  if (gaussian == F){ # by default, gaussian=F
    
    K2 <- D2K(as.matrix(d))
    
  }else{
    
    if (compdist2 == "BCD"){
      print("Gaussian kernel not recommended for BCD, using D2K instead")
      K2 <- D2K(as.matrix(d))
      
    }else{
      alldist <- as.numeric(d)
      K2 <- exp(-as.matrix(d^2)/median(alldist)^2) # Gaussian Kernel using Aitchson distance
    }
  }
  

  # Combine K1 and K2 Kernels
  if (combinethru=="sum"){
    K <- K1 + K2 # Sum
  } else {
    K <- K1 %*% K2 # Product
  }
  
  
  # Defining Kernel for the predictor variable
  if (Xdist == "categorical"){
    d2 <- outer(xvec, xvec, FUN = function(x,y){1 - as.numeric(x == y)})
    L <- D2K(as.matrix(d2))
  }
  
  
  m <- nrow(K)
  H <- diag(m) - (1/m) * matrix(1, nrow=m, ncol=m) #centering matrix
  
  
  if ("test" %in% result){
    
    kdcobs <- (1/(m^2))*sum(K * H %*% L %*% H)
    #kdcobs <- (1/(m^2))*sum(K %*% H %*% L %*% H)
    
    #Do permutation
    kdcstore <- numeric(num.perm)
    
    for (i in 1:num.perm){
      ind <- sample(m)
      K.curr <- K[ind, ind]
      
      kdcstore[i] <- (1/(m^2))*sum(K.curr * H %*% K2 %*% H)
    }
    pval <- (sum(kdcstore[-1] >= kdcobs)+1)/(num.perm)
  }else{
    kdcobs <- pval <- NA
  }
  
  return(list(K1=K1, K2=K2, K=K, kdcstat = kdcobs, pval = pval))
}
