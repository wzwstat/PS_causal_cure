##===============================================================================
## Native Method in Table 1 and Table 2
##===============================================================================

###############################################################################
#  Running this script you will save the results for 
#  Native Method: ACE_est, SPCE_cure_est

#  The script by parallel for one setting may take ~2 hours to execute. 
#  To save time, the results for both scenarios have been saved 
#  in the folder "MSE_intermediate_results/..":
#  1.-Table 1 and Table 2 for weak association in the manuscript
#     Native_mse_weak.RData
#     ......
#
#  2.-Table 1 and Table 2 for strong association in the manuscript
#     Native_mse_strong.RData
#     ......
#
#  The saved results can be loaded with load(). 
#  Example:load("MSE_intermediate_results/weak/Native_mse_weak.RData")
#          load("MSE_intermediate_results/strong/Native_mse_strong.RData")
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
#              (Native_mse_weak.RData)
#
#          (b) To obtain the all results for strong association for  
#              gamma0_true = 1.5 (PE = 20%), 1 (PE = 35%) and 0 (PE = 50%)
#              alpha0_true = 0, 0.5, 1,  
#              and alpha11_true = beta11_true = log(2)  
#              to obtain the other results for strong association:
#              (Native_mse_strong.RData)
###############################################################################

#setwd()
source("functions/function_definitions.R")
load("true/True_value.RData")

# Necessary packages ------------------------------------------------------
packagelist = c('mixcure','foreign','parallel','doParallel')
lapply(packagelist, function(x) do.call("require", list(x)))
timestart <- Sys.time()

Getace <- function (ii){
  set.seed(30*ii)
  
  data_ran <- getdata()
  
  Native <- function(data){
    spce <- ace <- rep(0,J)
    tgrid <- seq(min(data$obse_time), max(data$obse_time), length = 101)[-1]
    
    for(j in 1:J){
      DATA <- datastratify(data)[[j]]
      DATA1 <- DATA[DATA$z == 1, ]
      DATA0 <- DATA[DATA$z == 0, ]
      
      pd1 <- mixcure(Surv(obse_time, delta) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,
                     ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = DATA1 )
      pd0 <- mixcure(Surv(obse_time, delta)~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,
                     ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = DATA0)
      
      SS1 <- predict(pd1, newdata = DATA, times = tgrid)
      SS0 <- predict(pd0, newdata = DATA, times = tgrid)
      surv1_t <- surv0_t <- matrix(0, nrow = nrow(DATA), ncol = length(tgrid)) 
      
      for(i in 1:nrow(DATA)){
        surv1_t[i, ] <- as.numeric( SS1$surv[[i]]) 
        surv0_t[i, ] <- as.numeric( SS0$surv[[i]]) 
      }
      a1 <- apply(surv1_t, 2, mean)
      a0 <- apply(surv0_t, 2, mean)
      ace[j] <- trapezoid(tgrid, a1) - trapezoid(tgrid, a0)
      spce[j] <-  min(a1[a1 != 0]) - min(a0[a0 != 0])
    }
    
    ACE <- mean(ace)
    SPCE <- mean(spce)
    return(list(ACE_est = ACE, SPCE_est = SPCE))
  }
  summ_est <- Native(data_ran)
  return(list(summ_est = summ_est))
}

######----------------------SET--------------------------------------------------
J <- 1
mu1 <- 2
mu0 <- 6
u1 <- 5
lamb <- 0.5
tau <- u1-0.5
nn <- 500          #sample size
repli <- 500

#======= alpha0 = 0, PE = 20%, alpha11=beta11=log(2) weak association ==================================
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

##====== alpha0=0, PE=35%, weak association =======================
gamma0_true <- 1
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE35 <- Get_MSE()

##====== alpha0=0, PE=50%, weak association =======================
gamma0_true <- 0
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE50 <- Get_MSE()

#------------- summary ---------------------------------------------
Results_mse_alp0 <- list(
  PE20 = Results_mse_PE20$Summary_summ,
  PE35 = Results_mse_PE35$Summary_summ,
  PE50 = Results_mse_PE50$Summary_summ
)
#-------------------------------------------------------------------

##====== alpha0 = 0.5, PE = 20%, weak association===================
summ_true <- summ_weak_alpha05[c(1,4)]

gamma0_true <- 1.5
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 0.5
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE20 <- Get_MSE()

##====== alpha0 = 0.5, PE = 35%, weak association===================
gamma0_true <- 1
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 0.5
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 
bet <- c(beta1_true, beta2_true, beta3_true, beta4_true, beta5_true,
         beta6_true, beta7_true, beta8_true, beta9_true, beta10_true, beta11_true)

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE35 <- Get_MSE()

##====== alpha0 = 0.5, PE = 50%, weak association=====================
gamma0_true <- 0
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 0.5
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE50 <- Get_MSE()

#------------- summary -----------------------------------------------
Results_mse_alp05 <- list(
  PE20 = Results_mse_PE20$Summary_summ,
  PE35 = Results_mse_PE35$Summary_summ,
  PE50 = Results_mse_PE50$Summary_summ
)
#-------------------------------------------------------------------

##====== alpha0 = 1, PE = 25%, weak association ====================
summ_true <- summ_weak_alpha1[c(1,4)]

gamma0_true <- 1.5
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 1
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE20 <- Get_MSE()

##====== alpha0 = 1, PE = 35%, weak association======================
gamma0_true <- 1
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 1
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE35 <- Get_MSE()

##====== alpha0 = 1, PE = 50%, weak association========================
gamma0_true <- 0
gam <- c(gamma0_true, gamma1_true, gamma2_true, gamma3_true, gamma4_true,
         gamma5_true, gamma6_true, gamma7_true, gamma8_true, gamma9_true, gamma10_true)
alpha0_true <- 1
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

Results_mse_PE50 <- Get_MSE()

#-----------summary-----------------------------------------------------
Results_mse_alp1 <- list(
  PE20 = Results_mse_PE20$Summary_summ,
  PE35 = Results_mse_PE35$Summary_summ,
  PE50 = Results_mse_PE50$Summary_summ
)
#-----------summary_weak------------------------------------------------

Results_mse_weak <- list(
  alp0 = as.data.frame(Results_mse_alp0),
  alp05 = as.data.frame(Results_mse_alp05),
  alp1 = as.data.frame(Results_mse_alp1)
)


#save(Results_mse_weak, file = "MSE_intermediate_results/weak/Native_mse_weak.RData")

