##===============================================================================
#  Script:      Table1_Supplementary.R
#
#  Authors:     Ziwen Wang | Chenguang Wang | Xiaoguang Wang
#
#  Purpose:     Create all results of Table 1 in Supplementary Material 
#
#  Note:        This script may take ~40 mins to execute
##===============================================================================

## Source functions and load true value -----------------------------------------
source("functions/function_definitions.R")
load("true/True_value.RData")

## Necessary package ------------------------------------------------------------
packagelist = c('condSURV','npcure','foreign','parallel','doParallel')
lapply(packagelist, function(x) do.call("require", list(x)))
timestart <- Sys.time()

Getace <- function (ii){
  set.seed(30 * ii)
  time0 <- Sys.time()
  
  # generate survival data
  data_ran <- getdata()
  
  PS_Beran <- function(data){
    ace <- spce_cure <- rep(0, J)
    tgrid <- seq(min(data$obse_time), max(data$obse_time), length = 101)[-1]
    
    for(j in 1:J){
      # stratification using propensity score
      DATA <- datastratify(data)[[j]]
      DATA1 <- DATA[DATA$z == 1, ]
      DATA0 <- DATA[DATA$z == 0, ]
      
      # calculate the bandwidth
      hn <- cc * (nrow(DATA)) ^ (-1 / 5)
       
      # Beran estimator hat_S(t|Xi)
      surv1_t <- surv0_t <- matrix(0, nrow(DATA), length(tgrid)) #hat S(t|Xi)
      for (ix in 1:nrow(DATA)) {
        for(it in 1:length(tgrid)){
          surv1_t[ix, it] <- Beran(time = DATA1$obse_time, status = DATA1$delta,
                                   covariate = DATA1$pscore, x = DATA$pscore[ix], y = tgrid[it], bw = hn )
          surv0_t[ix, it] <- Beran(time = DATA0$obse_time, status = DATA0$delta,
                                   covariate = DATA0$pscore, x = DATA$pscore[ix], y = tgrid[it], bw = hn )
        }
      }
      
      # Estimate the area under the survival curve using trapezoidal integration
      surv1_area <- surv0_area <- rep(0, length(tgrid))
      for (ix in 1:nrow(DATA)) {
        surv1_area[ix] <- trapezoid(tgrid, surv1_t[ix, ])
        surv0_area[ix] <- trapezoid(tgrid, surv0_t[ix, ])
      }
      ace[j] <- mean(surv1_area - surv0_area)
      
      # Estimate the cure probability
      cure1_x <- probcure(pscore, obse_time, delta, DATA1, x0 = DATA$pscore, h = rep(hn, 100))$q
      cure0_x <- probcure(pscore, obse_time, delta, DATA0, x0 = DATA$pscore, h = rep(hn, 100))$q
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
cc <- 0.05
J <- 5
u1 <- 5   # upper bound of censored time, Tau_C
tau <- u1 - 0.5
nn <- 500         # sample size
repli <- 500       # the number of replications

mu1 <- 2  #used in the baseline hazard function under treat z=1
mu0 <- 6  #used in the baseline hazard function under treat z=0
lamb <- 0.5  #used in the baseline hazard function

#set the coefficients of propensity score 
gamma0_true <- 0
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

##------- Parallel --------------------------------------
no_cores <- detectCores(logical = FALSE) 
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

##====== weak association ===============================
summ_true <- summ_weak_alpha0[c(1,4)]

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

MySummary_C005 <- Get_MSE()

C <- 0.1
MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)
MySummary_C01 <- Get_MSE()

C <- 0.3
MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)
MySummary_C03 <- Get_MSE()

C <- 0.5
MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)
MySummary_C05 <- Get_MSE()

C <- 1
MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)
MySummary_C1 <- Get_MSE()

load("simulation_study/MSE_intermediate_results/weak/MSE_ACE_SPCEcure_weak_alp0_PE50_PSBeran.RData")
Results_hx_weak <- Get_MSE()

# -------------summary------------------------------------
Results_mse_weak<-list(
  h_x = Results_hx_weak$Summary_summ,
  C005 = MySummary_C005$Summary_summ,
  C01 = MySummary_C01$Summary_summ,
  C03 = MySummary_C03$Summary_summ,
  C05 = MySummary_C05$Summary_summ,
  C1 = MySummary_C1$Summary_summ)

## ====== strong association ===============================
summ_true <- summ_strong_alpha0[c(1,4)]
alpha11_true <- log(2)
alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, alpha5_true,
         alpha6_true, alpha7_true, alpha8_true, alpha9_true, alpha10_true, alpha11_true) 
beta11_true <- log(2)
bet <- c(beta1_true, beta2_true, beta3_true, beta4_true, beta5_true,
         beta6_true, beta7_true, beta8_true, beta9_true, beta10_true, beta11_true)

C <- 0.05
MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

MySummary_C005 <- Get_MSE()

C <- 0.1
MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)
MySummary_C01 <- Get_MSE()

C <- 0.3
MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)
MySummary_C03 <- Get_MSE()

C <- 0.5
MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)
MySummary_C05 <- Get_MSE()

C <- 1
MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)
MySummary_C1 <- Get_MSE()

load("simulation_study/MSE_intermediate_results/strong/MSE_ACE_SPCEcure_strong_alp0_PE50_PSBeran.RData")
Results_hx_strong <- Get_MSE()

stopImplicitCluster() 

#-------------summary-------------------------------------
Results_mse_strong <- list(
  h_x = Results_hx_strong$Summary_summ,
  C005 = MySummary_C005$Summary_summ,
  C01 = MySummary_C01$Summary_summ,
  C03 = MySummary_C03$Summary_summ,
  C05 = MySummary_C05$Summary_summ,
  C1 = MySummary_C1$Summary_summ
)

## ===============================================================================
##   The results about table 1 in Supplementary
## ===============================================================================
print(Results_mse_weak)
print(Results_mse_strong)

timeend <- Sys.time()
runningave <- timeend - timestart

#save(Results_mse_weak, Results_mse_strong, file = "simulation_study/results/Table1_supplementary.Rdata")

