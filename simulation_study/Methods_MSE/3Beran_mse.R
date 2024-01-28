##===============================================================================
## Beran Method in Table 1 and Table 2
##===============================================================================

###############################################################################
#  Running this script you will save the results for
#  Beran Method: ACE_est, SPCE_cure_est
 
#  The script by parallel for one setting may take ~13 hours to execute. 
#  To save time, the results for both scenarios have been saved 
#  in the folder "MSE_intermediate_results/..":
#  1.-Table 1 and Table 2 for weak association in the manuscript
#      (1) MSE_ACE_SPCEcure_weak_alp0_PE20_Beran.RData
#      (2) MSE_ACE_SPCEcure_weak_alp0_PE35_Beran.RData
#      (3) MSE_ACE_SPCEcure_weak_alp0_PE50_Beran.RData
#      (4) MSE_ACE_SPCEcure_weak_alp05_PE20_Beran.RData
#      (5) MSE_ACE_SPCEcure_weak_alp05_PE35_Beran.RData
#      (6) MSE_ACE_SPCEcure_weak_alp05_PE50_Beran.RData
#      (7) MSE_ACE_SPCEcure_weak_alp1_PE20_Beran.RData
#      (8) MSE_ACE_SPCEcure_weak_alp1_PE35_Beran.RData
#      (9) MSE_ACE_SPCEcure_weak_alp1_PE50_Beran.RData
#      ......
#  2.-Table 1 and Table 2 for strong association in the manuscript
#      (1) MSE_ACE_SPCEcure_strong_alp0_PE20_Beran.RData
#      (2) MSE_ACE_SPCEcure_strong_alp0_PE35_Beran.RData
#      (3) MSE_ACE_SPCEcure_strong_alp0_PE50_Beran.RData
#      (4) MSE_ACE_SPCEcure_strong_alp05_PE20_Beran.RData
#      (5) MSE_ACE_SPCEcure_strong_alp05_PE35_Beran.RData
#      (6) MSE_ACE_SPCEcure_strong_alp05_PE50_Beran.RData
#      (7) MSE_ACE_SPCEcure_strong_alp1_PE20_Beran.RData
#      (8) MSE_ACE_SPCEcure_strong_alp1_PE35_Beran.RData
#      (9) MSE_ACE_SPCEcure_strong_alp1_PE50_Beran.RData
#      ......
#
#  The saved results can be loaded with load(). 
#  Example: load("MSE_intermediate_results/weak/MSE_ACE_SPCEcure_weak_alp0_PE20_Beran.RData")
#  
#  Using these intermediate results, Table 1 and Table 2 in the paper
#  can be represented using the master script in "MSE_table1&2.R" 
#  without manual intervention
#################################################################################

###############################################################################
#  Remark 1: The following code gives the results for gamma0_true = 1.5 (PE = 20%), 
#            alpha0_true = 0, and alpha11_true = beta11_true = log(1.5)    
#            (MSE_ACE_SPCEcure_weak_alp0_PE20_Beran.RData)
#
#            Change gamma0_true = 1, 0 alpha0_true = 0.5, 1 in lines 144 - 190 to obtain the other results for weak association:
#            MSE_ACE_SPCEcure_weak_alp0_PE35_Beran.RData
#            MSE_ACE_SPCEcure_weak_alp0_PE50_Beran.RData
#            MSE_ACE_SPCEcure_weak_alp05_PE35_Beran.RData
#            MSE_ACE_SPCEcure_weak_alp05_PE50_Beran.RData
#            MSE_ACE_SPCEcure_weak_alp1_PE35_Beran.RData
#            MSE_ACE_SPCEcure_weak_alp1_PE50_Beran.RData
###############################################################################
#  Remark 2: To obtain the results for strong association 
#           (MSE_ACE_SPCEcure_strong_alp0_PE35_Beran.RData)
#            In this script, change the values of alpha11_true and beta11_true in lines 44-190
# 
#            Change gamma0_true = 1, 0 and alpha0_true = 0.5, 1 to obtain the results for 
#            MSE_ACE_SPCEcure_strong_alp05_PE35_Beran.RData
#            MSE_ACE_SPCEcure_strong_alp05_PE50_Beran.RData
#            MSE_ACE_SPCEcure_strong_alp1_PE35_Beran.RData
#            MSE_ACE_SPCEcure_strong_alp1_PE50_Beran.RData
###############################################################################

#setwd()
source("functions/function_definitions.R")
load("true/True_value.RData")

# Necessary packages ------------------------------------------------------
packagelist = c('condSURV', 'npcure', 'foreign', 'parallel', 'doParallel')
lapply(packagelist, function(x) do.call("require", list(x)))
timestart <- Sys.time()

Getace <- function (ii){
  set.seed(30*ii)
  time0 <- Sys.time()
  
  # generate survival data
  data_ran <- getdata()
  
  PS_Beran <- function(data){
    ace <- spce_cure <- rep(0, J)
    tgrid <- seq(min(data$obse_time), max(data$obse_time), length = 101)[-1]
    
    for(j in 1:J){

      DATA <- datastratify(data)[[j]]
      DATA1 <- DATA[DATA$z == 1,]
      DATA0 <- DATA[DATA$z == 0,]
      
      # calculate the bandwidth
      hcv <- berancv(pscore, obse_time, delta, DATA, x0 = DATA$pscore)
      hncure <- probcurehboot(pscore, obse_time, delta, DATA, x0 = DATA$pscore)
      
      # Beran estimator hat_S(t|Xi)
      surv1_t <- surv0_t <- matrix(0, nrow(DATA), length(tgrid)) #hat S(t|Xi)
      for (ix in 1:nrow(DATA)) {
        for(it in 1:length(tgrid)){
          surv1_t[ix, it] <- Beran(time = DATA1$obse_time, status = DATA1$delta,
                                   covariate = DATA1$pscore, x = hcv$x0[ix], y = tgrid[it], bw = hcv$h[ix] )
          surv0_t[ix, it] <- Beran(time = DATA0$obse_time, status = DATA0$delta,
                                   covariate = DATA0$pscore, x = hcv$x0[ix], y = tgrid[it], bw = hcv$h[ix] )
        }
      }
      
      # Estimate the area under the survival curve using trapezoidal integration
      surv1_area <- surv0_area <- rep(0, length(tgrid))
      for(it in 1:length(tgrid) ){
        surv1_area[it] <- trapezoid(tgrid, surv1_t[it,])
        surv0_area[it] <- trapezoid(tgrid, surv0_t[it,])
      }
      ace[j] <- mean(surv1_area - surv0_area)
      
      # Estimate the cure probability
      cure1_x <- probcure(pscore, obse_time, delta, DATA1, x0 = hncure$x0, h = hncure$h)$q
      cure0_x <- probcure(pscore, obse_time, delta, DATA0, x0 = hncure$x0, h = hncure$h)$q
      spce_cure[j] <- mean(cure1_x) - mean(cure0_x)
    }
    
    ACE <- mean(ace)
    SPCE_cure <- mean(spce_cure)
    return(list(ACE_est = ACE, SPCE_cure_est = SPCE_cure))
  }
  
  summ_est <- PS_Beran(data_ran)
  return(list(summ_est = summ_est))
}

######----------------------SET--------------------------------------------------
J <- 1
u1 <- 5   # upper bound of censored time, Tau_C
tau <- u1 - 0.5
nn <- 500         # sample size
repli <- 500       # the number of replications

mu1 <- 2  #used in the baseline hazard function under treat z=1
mu0 <- 6  #used in the baseline hazard function under treat z=0
lamb <- 0.5  #used in the baseline hazard function

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

summ_true <- summ_weak_alpha0[c(1,4)]

##-------Parallel--------------------------------------
no_cores <- detectCores(logical = FALSE) 
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

MySummary <- Get_MSE()
stopImplicitCluster()

MySummary

time1 <- Sys.time()
runningave <- time1 - timestart


##----------save results --- weak association --------------------------------------------------------------------------
#save(MyRsult, file = "MSE_intermediate_results/weak/MSE_ACE_SPCEcure_weak_alp0_PE20_Beran.RData")
