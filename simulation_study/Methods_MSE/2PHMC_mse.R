##===============================================================================
## PHMC (CoxPH) Method in Table 1 and Table 2
##===============================================================================

###############################################################################
#  Running this script you will save the results for 
#  PHMC Method: ACE_est, SPCE_cure_est

#  The script by parallel for all settings may take ~30 mins to execute. 
#  To save time, the results for both scenarios have been saved 
#  in the folder "MSE_intermediate_results/..":
#  1.-Table 1 and Table 2 for weak association in the manuscript
#     PHMC_mse_weak.RData
#     ......
#
#  2.-Table 1 and Table 2 for strong association in the manuscript
#     PHMC_mse_strong.RData
#     ......
#
#  The saved results can be loaded with load(). 
#  Example: load("MSE_intermediate_results/weak/PHMC_mse_weak.RData")
#  
#  Using these intermediate results, Table 1 and Table 2 in the paper
#  can be represented using the master script in "MSE_table1&2.R"
#  without manual intervention
#################################################################################

###############################################################################
#  Remark: (a) The following code gives the all results for  
#              gamma0_true = 1.5 (PE = 20%), 1 (PE = 35%) and 0 (PE = 50%)
#              alpha0_true = 0, 0.5, 1,  
#              and alpha11_true = beta11_true = log(1.5)  
#              to obtain the other results for weak association:
#              (PHMC_mse_weak.RData)
#
#          (b) To obtain the all results for strong association for  
#              gamma0_true = 1.5 (PE = 20%), 1 (PE = 35%) and 0 (PE = 50%)
#              alpha0_true = 0, 0.5, 1,  
#              and alpha11_true = beta11_true = log(2)  
#              to obtain the other results for strong association:
#              (PHMC_mse_strong.RData)
###############################################################################

#setwd()
source("functions/function_definitions.R")
load("true/True_value.RData")

# Necessary packages ------------------------------------------------------
packagelist = c('smcure', 'foreign','parallel','doParallel')
lapply(packagelist, function(x) do.call("require", list(x)))
timestart <- Sys.time()

# replication ------------------------------------------------------------------
Getace <- function (ii){
  set.seed(30 * ii)
  
  data.ran <- getdata() 
  
  Cox_est <- function(DATA){
 
    DATA <- DATA[order(DATA$obse_time), ]
    phmc <- smcure( Surv(obse_time,delta) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+z,
                    cureform = ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+z,
                    data = DATA, model = "ph", Var = FALSE)
   
    #---Estimate the survival probability under treat Z=1 ------------------
    Z1 <- cbind(DATA$x1, DATA$x2, DATA$x3, DATA$x4, DATA$x5, DATA$x6, 
                DATA$x7, DATA$x8, DATA$x9, DATA$x10, 1)
    pred1 <- predictsmcure(phmc, newX = Z1, newZ = Z1, model = "ph")
    
    id <- dim(pred1$prediction)[2]
    t1 <- as.numeric(pred1$prediction[, id])
    surv1_t <- pred1$prediction[,-id]
    Cox_S1_t <- as.numeric(apply(surv1_t, 1, mean))
    
    #---Estimate the survival probability under treat Z=0 ------------------
    Z0 <- cbind(DATA$x1, DATA$x2, DATA$x3, DATA$x4, DATA$x5, DATA$x6, 
                DATA$x7, DATA$x8, DATA$x9, DATA$x10, 0)
    pred0 <- predictsmcure(phmc, newX = Z0, newZ = Z0, model = "ph")
    
    id <- dim(pred0$prediction)[2]
    t0 <- as.numeric(pred0$prediction[, id])
    surv0_t <- pred0$prediction[,-id]
    Cox_S0_t <- as.numeric(apply(surv0_t, 1, mean))
    
    #---Caculate the area under the survival curves ------------------
    surv1_area <- trapezoid(t1, Cox_S1_t)
    surv0_area <- trapezoid(t0, Cox_S0_t)
    
    #---caculate the RACE and SPCE for cure rate---------------------------------
    ace <- surv1_area - surv0_area
    spce <- min(Cox_S1_t[Cox_S1_t!=0]) - min(Cox_S0_t[Cox_S0_t!=0])
    
    return(list(ACE_est = ace, SPCE_est = spce))
  }
  summ_est <- Cox_est(data.ran)
  return(list(summ_est = summ_est))
}

######----------------------SET--------------------------------------------------
u1 <- 5   # Upper bound of censored time, Tau_C
tau <- u1-0.5
nn <- 500         #sample size
repli <- 500       # the number of replications

mu1 <- 2  #used in the baseline hazard function under treat z=1
mu0 <- 6  #used in the baseline hazard function under treat z=0
lamb <- 0.5  #used in the baseline hazard function 

#======= alpha0 = 0, PE = 20%, alpha11=beta11=log(2) strong association ==================================
summ_true <- summ_weak_alpha0[c(1,4)]

#set the coefficients of propensity score 
gamma0_true <- 1.5
gamma1_true <- log(1.5)
gamma2_true <- -log(1.5)
gamma3_true <- log(1.5)
gamma4_true <- -log(1.5)
gamma5_true <- log(1.5)
gamma6_true <- -log(1.5)
gamma7_true <- log(1.5)
gamma8_true <- -log(1.5)
gamma9_true <-  0
gamma10_true <- 0
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)

# set the coefficients of cure rate 
# alpha11_true is the exposure effect
alpha0_true <- 0
alpha1_true <- log(1.5)
alpha2_true <- -log(1.5)
alpha3_true <- log(1.5)
alpha4_true <- -log(1.5)
alpha5_true <- log(2)
alpha6_true <- log(2)
alpha7_true <- log(1.5)
alpha8_true <- log(1.5)
alpha9_true <-  log(1.25)
alpha10_true <- log(1.25)
alpha11_true <- log(1.5)
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 

# set the coefficients of survival model for the uncured subjects
# beta11_true is the exposure effect
beta1_true <- log(1.5)
beta2_true <- -log(1.5)
beta3_true <- log(1.5)
beta4_true <- -log(1.5)
beta5_true <- log(2)
beta6_true <- log(2)
beta7_true <- log(1.5)
beta8_true <- log(1.5)
beta9_true <- log(1.25)
beta10_true <- log(1.25)
beta11_true <- log(1.5) 
bet <- c(beta1_true, beta2_true, beta3_true, beta4_true, beta5_true,
         beta6_true, beta7_true, beta8_true, beta9_true, beta10_true, beta11_true)

#------- Parallel------------------------------------------------------------
no_cores <- detectCores(logical = FALSE)  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE20 <- Get_MSE()

stopImplicitCluster()  


##====== alpha0=0, PE=35%, strong association ===============
gamma0_true <- 1
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 0
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 
bet <- c(beta1_true, beta2_true, beta3_true, beta4_true, beta5_true,
         beta6_true, beta7_true, beta8_true, beta9_true, beta10_true, beta11_true)


no_cores <- detectCores(logical = FALSE)  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE35 <- Get_MSE()

stopImplicitCluster()  


##====== alpha0=0, PE=50%, strong association ==================
gamma0_true <- 0
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 0
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 
bet <- c(beta1_true, beta2_true, beta3_true, beta4_true, beta5_true,
         beta6_true, beta7_true, beta8_true, beta9_true, beta10_true, beta11_true)


no_cores <- detectCores(logical = FALSE)  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE50 <- Get_MSE()

stopImplicitCluster()  

#-------------summary-------------------------------------
Results_mse_alp0 <- list(PE20 = Results_mse_PE20$Summary_summ,
                         PE35 = Results_mse_PE35$Summary_summ,
                         PE50 = Results_mse_PE50$Summary_summ)
#----------------------------------------------------------

##====== alpha0 = 0.5, PE = 20%, strong association===============
summ_true <- summ_weak_alpha05[c(1,4)]

gamma0_true <- 1.5
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 0.5
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 
bet <- c(beta1_true, beta2_true, beta3_true, beta4_true, beta5_true,
         beta6_true, beta7_true, beta8_true, beta9_true, beta10_true, beta11_true)

no_cores <- detectCores(logical = FALSE)  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE20 <- Get_MSE()

stopImplicitCluster()  

##====== alpha0 = 0.5, PE = 35%, strong association===============
gamma0_true <- 1
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 0.5
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 
bet <- c(beta1_true, beta2_true, beta3_true, beta4_true, beta5_true,
         beta6_true, beta7_true, beta8_true, beta9_true, beta10_true, beta11_true)


no_cores <- detectCores(logical = FALSE)  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE35 <- Get_MSE()

stopImplicitCluster()  


##====== alpha0 = 0.5, PE = 50%, strong association==================
gamma0_true <- 0
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 0.5
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 
bet <- c(beta1_true, beta2_true, beta3_true, beta4_true, beta5_true,
         beta6_true, beta7_true, beta8_true, beta9_true, beta10_true, beta11_true)


no_cores <- detectCores(logical = FALSE)  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE50 <- Get_MSE()

stopImplicitCluster()  

#-------------summary-------------------------------------
Results_mse_alp05 <- list(PE20 = Results_mse_PE20$Summary_summ,
                          PE35 = Results_mse_PE35$Summary_summ,
                          PE50 = Results_mse_PE50$Summary_summ)
#----------------------------------------------------------

##====== alpha0 = 1, PE = 20%, strong association ==================
summ_true <- summ_weak_alpha1[c(1,4)]

gamma0_true <- 1.5
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 1
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 
bet <- c(beta1_true, beta2_true, beta3_true, beta4_true, beta5_true,
         beta6_true, beta7_true, beta8_true, beta9_true, beta10_true, beta11_true)


no_cores <- detectCores(logical = FALSE)  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE20 <- Get_MSE()

stopImplicitCluster()  

##====== alpha0 = 1, PE = 35%, strong association===============
gamma0_true <- 1
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 1
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 
bet <- c(beta1_true, beta2_true, beta3_true, beta4_true, beta5_true,
         beta6_true, beta7_true, beta8_true, beta9_true, beta10_true, beta11_true)

no_cores <- detectCores(logical = FALSE)  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE35 <- Get_MSE()

stopImplicitCluster()  

##====== alpha0 = 1, PE = 50%, strong association===============
gamma0_true <- 0
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 1
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 
bet <- c(beta1_true, beta2_true, beta3_true, beta4_true, beta5_true,
         beta6_true, beta7_true, beta8_true, beta9_true, beta10_true, beta11_true)


no_cores <- detectCores(logical = FALSE)
cl <- makeCluster(no_cores)
registerDoParallel(cl)  

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE50 <- Get_MSE()

stopImplicitCluster()  

#------------- summary -------------------------------------
Results_mse_alp1 <- list(
  PE20 = Results_mse_PE20$Summary_summ,
  PE35 = Results_mse_PE35$Summary_summ,
  PE50 = Results_mse_PE50$Summary_summ
)
#----------------------------------------------------------

Results_mse <- list(
  alp0 = as.data.frame(Results_mse_alp0),
  alp05 = as.data.frame(Results_mse_alp05),
  alp1 = as.data.frame(Results_mse_alp1)
)

##----------save results --- weak association --------------------------------------------------------------------------
#save(Results_mse_weak, file = "MSE_intermediate_results/weak/PHMC_mse_weak.RData")
