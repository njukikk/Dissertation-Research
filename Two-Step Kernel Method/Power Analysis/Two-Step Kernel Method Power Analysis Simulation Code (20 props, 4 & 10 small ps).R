
# Power Analysis for the Two-Step Kernel Method (20 props, 4 & 10 small ps Compositional Data)

### Source files
source("D:/A.School/OSU/Class/1. Dissertation/6. Fall 2024/Simulations/Two-Step Kernel Method/CODAK_functions.R")
source("D:/A.School/OSU/Class/1. Dissertation/6. Fall 2024/Simulations/Two-Step Kernel Method/two_step_KDC_v2.R")
source("D:/A.School/OSU/Class/1. Dissertation/6. Fall 2024/Simulations/Two-Step Kernel Method/simulation_functions.R")
source("D:/A.School/OSU/Class/1. Dissertation/6. Fall 2024/Simulations/Two-Step Kernel Method/final_addeffects_fn_version6.R")

# Set working directory
setwd("D:/A.School/OSU/Class/1. Dissertation/6. Fall 2024/Simulations/Two-Step Kernel Method/Dist1 = Hamming, Dist2 = Aitchison/20 ps 4 small ps")

#Proportions for generating data without zeros
props.kernel <- c(0.15, 0.13, 0.07, 0.05, 0.1, 0.05, 0.05, 0.05, 0.02, 0.1, 0.02, 
                  0.003, 0.002, 0.004, 0.12, 0.04, 0.004, 0.007, 0.02, 0.01)

#Comments this part when generating data with 10 small ps
#Proportions for generating data with zeros
props.pseudo <- props.kernel

props.pseudo[4] <- 9e-5
props.pseudo[6] <- props.pseudo[6] + 0.05 - 9e-5
props.pseudo[9] <- 7e-5
props.pseudo[11] <- props.pseudo[11] + 0.02 - 7e-5
props.pseudo[12] <- 5e-5
props.pseudo[15] <- props.pseudo[15] + 0.003 - 5e-5
props.pseudo[17] <- 6e-5
props.pseudo[18] <- props.pseudo[18] + 0.004 - 6e-5


# # Comments this part when generating data with 3 small ps
# #Proportions for generating data with zeros
# props.pseudo <- props.kernel
# 
# i <- 1
# while (i <= 20){
#   props.pseudo[i] <- props.kernel[i]*1e-3
#   props.pseudo[i+1] <- props.kernel[i+1] + props.kernel[i] - props.pseudo[i]
#   i <- i + 2
# }
############################################################################################################

num.cat <- length(props.pseudo)

reps <- 20
num.sim1 <- 10000
num.sim2 <- 500

#12 observations in each group 
ng <- 12
Xmat <- matrix(rep(c(0, 1), each = ng), ncol = 1)

#number of events for the individual samples
set.seed(92486)
nvec <- sample(c(30000, 40000, 50000), 2*ng, replace = T) # cells per individual samples

#covariance matrix of random effects
set.seed(92486)
sigbvec <- rep(0.2, num.cat-1) 
sigb <- diag(num.cat-1)*sigbvec^2

############################################################################################################

## Case 1: No difference in all
set.seed(123)
props.pseudo.new <- props.pseudo

betarev1 <- betarev2 <- getES.mlogit(props.pseudo)

betaX <- as.matrix(betarev2 - betarev1)

two.step.case1 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb, num.sim=num.sim1, num.perm = 10000,
                                          printafter = 1, Xdist = "categorical", method = "Two.Step", 
                                          closure = "No", compdist1 = "hamming", compdist2 = "aitchison", 
                                          combinethru="sum", standardize = TRUE)

zero.imputed.case1 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb, num.sim=num.sim1, num.perm = 10000,
                                              printafter = 1, Xdist = "categorical", method = "Zero.Imputed", 
                                              closure = "No", compdist1 = "hamming", compdist2 = "aitchison", 
                                              combinethru="sum", standardize = TRUE)


# Case 1 data: Type I error computation  
logalphavec <- seq(log10(0.00005), log10(0.05), 0.05)
alphavec <- 10^logalphavec

powervec.two.step.case1 <- powervec.zero.imputed.case1 <- numeric(length(alphavec))
for(k in 1:length(alphavec)){
  powervec.two.step.case1[k] <- mean(two.step.case1$pval <= alphavec[k], na.rm = T)
  powervec.zero.imputed.case1[k] <- mean(zero.imputed.case1$pval <= alphavec[k], na.rm = T)
}

# saving case 1 data for future use
case1_20props_4small_2step_vs_0imputation_initeffsize <- cbind.data.frame(alphavec, powervec.two.step.case1, 
                                                              powervec.zero.imputed.case1)

write.csv(case1_20props_4small_2step_vs_0imputation_initeffsize, 
          paste0(getwd(), "/case1_20props_4small_2step_vs_0imputation_initeffsize.csv"), row.names = TRUE)

############################################################################################################

## Case 2: Small difference in all

set.seed(123)

ds.pseudo.case2 <- numeric(reps)
pval.two.step.case2 <- matrix(NA, nrow = reps, ncol = num.sim2)
pval.zero.imputed.case2 <- matrix(NA, nrow = reps, ncol = num.sim2)

for (i in 1:reps){
  print(paste("reps =", i))
  
  betarev1 <- getES.mlogit(props.pseudo)
  props.pseudo.new <- add.new.effects(maxeffect = 0.002, prop = props.pseudo, n.nochange = 0, ind.nochange = NULL)
  
  betarev2 <- getES.mlogit(props.pseudo.new)
  betaX <- as.matrix(betarev2 - betarev1)
  
  set.seed(1)
  two.step.case2 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb, num.sim=num.sim2, num.perm = 10000,
                                            Xdist = "categorical", method = "Two.Step", closure = "No",
                                            compdist1 = "hamming", compdist2 = "aitchison", combinethru="sum",
                                            printafter = 1, standardize = TRUE)
  
  zero.imputed.case2 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb, num.sim=num.sim2, num.perm = 10000,
                                                Xdist = "categorical", method = "Zero.Imputed", closure = "No",
                                                compdist1 = "hamming", compdist2 = "aitchison", combinethru="sum",
                                                printafter = 1, standardize = TRUE)
  
  pval.two.step.case2[i,] <- two.step.case2$pval
  pval.zero.imputed.case2[i, ] <- zero.imputed.case2$pval
  
  ds.pseudo.case2[i] <- dist(rbind(props.pseudo, props.pseudo.new), method = 'aitchison')
}


### Case 2: Computing power values
cat.num <- nrow(pval.two.step.case2)
powervec.two.step.case2 <- powervec.zero.imputed.case2 <- numeric(cat.num)

for (i in 1:cat.num){
  powervec.two.step.case2[i] <- mean(pval.two.step.case2[i,] < 0.05)
  powervec.zero.imputed.case2[i] <- mean(pval.zero.imputed.case2[i,] < 0.05)
}

# saving case 2 data for future use
case2_20props_4small_2step_vs_0imputation_initeffsize <- cbind.data.frame(ds.pseudo.case2, powervec.two.step.case2, 
                                                              powervec.zero.imputed.case2)

write.csv(case2_20props_4small_2step_vs_0imputation_initeffsize, 
          paste0(getwd(), "/case2_20props_4small_2step_vs_0imputation_initeffsize.csv"), row.names = TRUE)

############################################################################################################

## Case 3: Small difference in some

set.seed(123)
ds.pseudo.case3 <- numeric(reps)
pval.two.step.case3 <- matrix(NA, nrow = reps, ncol = num.sim2)
pval.zero.imputed.case3 <- matrix(NA, nrow = reps, ncol = num.sim2)


for (i in 1:reps){
  
  print(paste("reps =", i))
  
  betarev1 <- getES.mlogit(props.pseudo)
  props.pseudo.new <- add.new.effects(maxeffect = 0.0015, prop = props.pseudo, n.nochange = 10, ind.nochange = 3:12)
  
  betarev2 <- getES.mlogit(props.pseudo.new)
  
  betaX <- as.matrix(betarev2 - betarev1)
  
  set.seed(1)
  two.step.case3 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb, num.sim=num.sim2, num.perm = 10000,
                                            Xdist = "categorical", method = "Two.Step", closure = "No",
                                            compdist1 = "hamming", compdist2 = "aitchison", combinethru="sum",
                                            printafter = 1, standardize = TRUE)
  
  zero.imputed.case3 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb, num.sim=num.sim2, num.perm = 10000,
                                                Xdist = "categorical", method = "Zero.Imputed", closure = "No",
                                                compdist1 = "hamming", compdist2 = "aitchison", combinethru="sum",
                                                printafter = 1, standardize = TRUE)
  
  pval.two.step.case3[i,] <- two.step.case3$pval
  pval.zero.imputed.case3[i, ] <- zero.imputed.case3$pval
  
  ds.pseudo.case3[i] <- dist(rbind(props.pseudo, props.pseudo.new), method = 'aitchison')
}

### Case 3: Computing power values
cat.num <- nrow(pval.two.step.case3)
powervec.two.step.case3 <- powervec.zero.imputed.case3 <- numeric(cat.num)

for (i in 1:cat.num){
  powervec.two.step.case3[i] <- mean(pval.two.step.case3[i,] < 0.05)
  powervec.zero.imputed.case3[i] <- mean(pval.zero.imputed.case3[i,] < 0.05)
}

# saving case 2 data for future use
case3_20props_4small_2step_vs_0imputation_initeffsize <- cbind.data.frame(ds.pseudo.case3, powervec.two.step.case3, 
                                                              powervec.zero.imputed.case3)

write.csv(case3_20props_4small_2step_vs_0imputation_initeffsize, 
          paste0(getwd(), "/case3_20props_4small_2step_vs_0imputation_initeffsize.csv"), row.names = TRUE)

############################################################################################################

## Case 4: Large difference in few

set.seed(123)
ds.pseudo.case4 <- numeric(reps)
pval.two.step.case4 <- matrix(NA, nrow = reps, ncol = num.sim2)
pval.zero.imputed.case4 <- matrix(NA, nrow = reps, ncol = num.sim2)

for (i in 1:reps){
  print(paste("reps =", i))
  
  betarev1 <- getES.mlogit(props.pseudo)
  props.pseudo.new <- add.new.effects(maxeffect = 0.003, prop = props.pseudo, n.nochange = 15, ind.nochange = 2:16)
  
  betarev2 <- getES.mlogit(props.pseudo.new)
  
  betaX <- as.matrix(betarev2 - betarev1)
  
  set.seed(1)
  two.step.case4 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb, num.sim=num.sim2, num.perm = 10000,
                                            Xdist = "categorical", method = "Two.Step", closure = "No",
                                            compdist1 = "hamming", compdist2 = "aitchison", combinethru="sum",
                                            printafter = 1, standardize = TRUE)
  
  zero.imputed.case4 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb, num.sim=num.sim2, num.perm = 10000,
                                                Xdist = "categorical", method = "Zero.Imputed", closure = "No",
                                                compdist1 = "hamming", compdist2 = "aitchison", combinethru="sum",
                                                printafter = 1, standardize = TRUE)
  
  pval.two.step.case4[i,] <- two.step.case4$pval
  pval.zero.imputed.case4[i, ] <- zero.imputed.case4$pval
  
  ds.pseudo.case4[i] <- dist(rbind(props.pseudo, props.pseudo.new), method = 'aitchison')
}


### Case 4: Computing power values
cat.num <- nrow(pval.two.step.case4)
powervec.two.step.case4 <- powervec.zero.imputed.case4 <- numeric(cat.num)

for (i in 1:cat.num){
  powervec.two.step.case4[i] <- mean(pval.two.step.case4[i,] < 0.05)
  powervec.zero.imputed.case4[i] <- mean(pval.zero.imputed.case4[i,] < 0.05)
}

# saving case 2 data for future use
case4_20props_4small_2step_vs_0imputation_initeffsize <- cbind.data.frame(ds.pseudo.case4, powervec.two.step.case4, 
                                                              powervec.zero.imputed.case4)

write.csv(case4_20props_4small_2step_vs_0imputation_initeffsize, 
          paste0(getwd(), "/case4_20props_4small_2step_vs_0imputation_initeffsize.csv"), row.names = TRUE)

###########################################################################################################

# saving all cases data for future use
ds.pseudo <- c(ds.pseudo.case2, ds.pseudo.case3, ds.pseudo.case4)
powervec.two.step.3cases <- c(powervec.two.step.case2, powervec.two.step.case3, powervec.two.step.case4)
powervec.zero.imputed.3cases <- c(powervec.zero.imputed.case2, powervec.zero.imputed.case3, 
                                  powervec.zero.imputed.case4)


all3cases_20props_4small_2step_vs_0imputation_initeffsize <- cbind.data.frame(ds.pseudo, powervec.two.step.3cases, 
                                                                  powervec.zero.imputed.3cases)
write.csv(all3cases_20props_4small_2step_vs_0imputation_initeffsize, 
          paste0(getwd(), "/all3cases_20props_4small_2step_vs_0imputation_initeffsize.csv"), row.names = TRUE)
