

###### Functions for simulating data for cell type abundance comparison #######

library("MASS")
library("coda.base")
library("lme4")
library("SummarizedExperiment")
library("diffcyt")
library("energy")
library("vegan")
library("MiRKAT")



############### Function to obtain effect sizes from one given proportion vector using multilogit ###################

#Getting the effect size betas in the mlogit scale starting from the prop scale

getES.mlogit <- function(prop, ref = "last"){
  
  # ref: reference cell type
  
  n <-length(prop)
  if (ref == "last"){ref <- n}
  
  beta <- numeric(n)
  beta[ref] <- 0
  
  beta[-ref] <- log(prop[-ref]/prop[ref])
  
  return(beta)
}



############### Function to obtain proportions from effect sizes using multilogit ###################

getProps.mlogit <- function(beta, betaX = NULL, X = NULL){
  
  # beta: vector of coefficients for categories
  # betaX: matrix of coefficients for the predictors:
  #       rows are categories, columns are different predictors, 
  # X is a column vector of the predictor values, num.pred x 1
  # Note that although we use one predictor of interest, 
  #       a covariate is essentially another predictor
  
  n <- length(beta) #number of categories
  
  if(is.null(betaX) == 1){
    betaX <- matrix(0, nrow = n,ncol = 1)
    X <- 0
  }
  
  ref <- which(beta == 0)
  
  prop <- numeric(n)
  
  prop[ref] <- 1
  prop[-ref] <- exp(as.numeric(beta[-ref] + as.matrix(betaX[-ref,])%*%X))*prop[ref]
  
  prop <- prop/sum(prop)
  return(prop)
}



######## Function to simulate proportions based on effect sizes #########

simProps <- function(beta, betaX, Xmat, nvec, sigb = NULL){
  
  # rows of Xmat denote sample units, columns denote different X variables
  # nvec: a vector containing total number of cells per file
  # sigb: covariance matrix of the random effects (num.cat x num.cat)
  
  num.cat <- length(beta)
  y <- matrix(NA, nrow = nrow(Xmat), ncol = num.cat)
  
  if(is.null(sigb)){
    for (i in 1:nrow(Xmat)){
      probs <- getProps.mlogit(beta, betaX, as.matrix(Xmat[i, ]))
      y[i, ] <- t(rmultinom(n = 1, size = nvec[i], prob = probs))
    }
  }else{
    bmat <- matrix(0, nrow = nrow(Xmat), ncol = num.cat)
    ref <- which(beta == 0)
    bmat[, ref] <- rep(0, nrow(Xmat))
    bmat[, -ref] <- mvrnorm(nrow(Xmat), mu = rep(0, num.cat-1), sigb)
    for (i in 1:nrow(Xmat)){
      probs <- getProps.mlogit(beta + bmat[i, ], betaX, as.matrix(Xmat[i, ]))
      y[i, ] <- t(rmultinom(n = 1, size = nvec[i], prob = probs))
    }
  }
  
  return(y)
}

##########################################################################

# Two-Step Kernel Method and Zero-Imputed Simulations


simulation_study_pseudo <- function(beta, betaX, Xmat, nvec, sigb, num.sim=10000, num.perm = 10000, 
                                    printafter = 100, Xdist = "categorical", method = "Two.Step", 
                                    closure = "No", compdist1 = "hamming", compdist2 = "aitchison", 
                                    combinethru="sum", standardize = TRUE){
  
  #Contains the random effect, sigb is the sd of the random effect
  #rows of Xmat denote sample units, columns denote different X variables
  #sigb is the covariance matrix of the random effects (num.cat x num.cat)
  
  pval.kernel <- numeric(num.sim)
  
  #### Simulation begins ####
  
  for (i in 1:num.sim){
    
    if (floor(i/printafter) == i/printafter){print(i)}
    
    #x is the predictor of interest. Works for only one predictor.
    x <- as.vector(Xmat)
    xvec <- as.vector(Xmat)
    
    if (method == "Two.Step"){
      
      y <- simProps(beta, betaX, Xmat, nvec, sigb)
      ng <- nrow(y)/2
      num.cat <- ncol(y)
      propmat <- y/rep(nvec, num.cat)
      
      out <- dcor.comp.2step(propmat, xvec, num.perm = num.perm, Xdist = Xdist, 
                             closure = closure, compdist1 = compdist1, compdist2 = compdist2, 
                             combinethru = combinethru, result = "test")
      
      pval.kernel[i] <- out$pval
      
    } else {
      
      y <- simProps(beta, betaX, Xmat, nvec, sigb)
      
      # Replacing the zeros with 1
      propmat.0.impute <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
      
      for (j in 1:length(y)){
        if (y[j] == 0) {propmat.0.impute[j] <- 1} else {propmat.0.impute[j] <- y[j]} 
      }
      
      # Taking closure after replacing zeros
      propmat <- propmat.0.impute/rowSums(propmat.0.impute)
      
      out <- dcor.comp(propmat, xvec, num.perm = num.perm, Xdist = Xdist, compdist = "aitchison",
                       gaussian = F, result = "test") 
      
      pval.kernel[i] <- out$pval
    }
  }  
  
  return(list(pval = pval.kernel))
}

##############################################################################




