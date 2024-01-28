##===============================================================================
## PS_Beran Method in Table 3 and Table 4
##===============================================================================

##################################################################################
#  Running this script you will save the results for
#    (a) Proposed estimator: ACE_est, SPCE1_est, SPCE2_est, SPCE_cure_est,
#    (b) the corresponding standard error: Boot_ace_se, Boot_spce1_se, 
#                                          Boot_spce2_se, Boot_spcecure_se
# 
#  The script by parallel for one setting may take ~72 hours to execute. 
#  To save time, the results for both scenarios have been saved 
#  in the folder "bootstrap_intermediate_results/..":
#  1.-Table 3 for weak association in the manuscript
#      (1) boots_weak_alp0_PE20_PSBeran_n150.RData
#      (2) boots_weak_alp0_PE20_PSBeran_n300.RData
#      (3) boots_weak_alp0_PE20_PSBeran_n500.RData
#      (4) boots_weak_alp0_PE35_PSBeran_n150.RData
#      (5) boots_weak_alp0_PE35_PSBeran_n300.RData
#      (6) boots_weak_alp0_PE35_PSBeran_n500.RData
#      (7) boots_weak_alp0_PE50_PSBeran_n150.RData
#      (8) boots_weak_alp0_PE50_PSBeran_n300.RData
#      (9) boots_weak_alp0_PE50_PSBeran_n500.RData
#      ......
#  2.-Table 4 for strong association in the manuscript
#      (1) boots_strong_alp0_PE20_PSBeran_n150.RData
#      (2) boots_strong_alp0_PE20_PSBeran_n300.RData
#      (3) boots_strong_alp0_PE20_PSBeran_n500.RData
#      (4) boots_strong_alp0_PE35_PSBeran_n150.RData
#      (5) boots_strong_alp0_PE35_PSBeran_n300.RData
#      (6) boots_strong_alp0_PE35_PSBeran_n500.RData
#      (7) boots_strong_alp0_PE50_PSBeran_n150.RData
#      (8) boots_strong_alp0_PE50_PSBeran_n300.RData
#      (9) boots_strong_alp0_PE50_PSBeran_n500.RData
#      ......
#
#  The saved results can be loaded with load(). 
#  Example:load("bootstrap_intermediate_results/boots_weak_alp0_PE20_PSBeran_n150.RData")
#  
#  Using these intermediate results, Table 3 and Table 4 in the paper
#  can be represented using the master scripts "table3.R" and "table4.R"
##################################################################################

##################################################################################
#  Remark 1: The following code gives the results for nn = 150, gamma0_true = 1.5 (PE = 20%), 
#            alpha0_true = 0 and alpha11_true = beta11_true = log(1.5)    
#            (boots_weak_alp0_PE20_PSBeran_n150.RData)
#
#            Change nn = 300, 500; gamma0_true = 1, 0 and alpha0_true = 0, 0.5, 1 in lines 250 - 295 
#            to obtain the other results for weak association:
#            boots_weak_alp0_PE35_PSBeran_n300.RData 
#            boots_weak_alp0_PE35_PSBeran_n500.RData 
#            boots_weak_alp0_PE50_PSBeran_n300.RData
#            boots_weak_alp0_PE50_PSBeran_n500.RData
#            ......
##################################################################################
#  Remark 2: To obtain the results for strong association 
#            (boots_strong_alp0_PE20_PSBeran_n150.RData)
#            In this script, change the values of alpha11_true and beta11_true in lines 250 - 295
# 
#            Change nn = 300, 500; gamma0_true = 1, 0 and alpha0_true = 0, 0.5, 1 to obtain the results for 
#            boots_strong_alp0_PE35_PSBeran_n300.RData 
#            boots_strong_alp0_PE35_PSBeran_n500.RData 
#            boots_strong_alp0_PE50_PSBeran_n300.RData
#            boots_strong_alp0_PE50_PSBeran_n500.RData
#            ......
##################################################################################

# Source functions and load true value ----------------------------------------- 
source("functions/function_definitions.R")
load("true/True_value.RData")

# Necessary packages -----------------------------------------------------------
packagelist = c('condSURV', 'npcure', 'foreign', 'parallel', 'doParallel')
lapply(packagelist, function(x) do.call("require", list(x)))
timestart <- Sys.time()

Getace <- function (ii){
  set.seed(30 * ii)
  
  # generate survival data
  data_ran <- getdata()
  
  PS_Beran <- function(data){
    ace <- spce1 <- spce2 <- spce_cure <- rep(0, J)
    tgrid <- seq(min(data$obse_time), max(data$obse_time), length = 101)[-1]
    
    for(j in 1:J){
      # stratification using propensity score
      DATA <- datastratify(data)[[j]]
      DATA1 <- DATA[DATA$z == 1,]
      DATA0 <- DATA[DATA$z == 0,]
      
      # calculate the bandwidth
      hcv <- berancv(pscore, obse_time, delta, DATA, x0 = DATA$pscore)
      hncure <- probcurehboot(pscore, obse_time, delta, DATA, x0 = DATA$pscore)
      
      # Beran estimator hat_S(t|Xi)
      surv1_t <- surv0_t <- matrix(0, nrow(DATA), length(tgrid))
      for(ix in 1:nrow(DATA)){
        for(it in 1:length(tgrid)){
          surv1_t[ix, it] <- Beran(time = DATA1$obse_time, status = DATA1$delta,
                                   covariate = DATA1$pscore, x = hcv$x0[ix], y = tgrid[it], bw = hcv$h[ix] )
          surv0_t[ix, it] <- Beran(time = DATA0$obse_time, status = DATA0$delta,
                                   covariate = DATA0$pscore, x = hcv$x0[ix], y = tgrid[it], bw = hcv$h[ix] )
        }
      }
      
      # Estimate the difference in area under the survival curve
      # Calculate the area using trapezoidal integration
      # Estimate the RACE
      ace[j] <- trapezoid(tgrid, apply(surv1_t, 2, mean)) - trapezoid(tgrid, apply(surv0_t, 2, mean))
      
      # Estimate the survival probability at time 0.8 and 2.4
      surv1_08 <- surv0_08 <- surv1_24 <- surv0_24 <- rep(0, nrow(DATA))
      for (ix in 1:nrow(DATA)) {
        surv1_08[ix] <- Beran(time = DATA1$obse_time, status = DATA1$delta,
                              covariate = DATA1$pscore, x = hcv$x0[ix], y = 0.8, bw = hcv$h[ix] )
        surv0_08[ix] <- Beran(time = DATA0$obse_time, status = DATA0$delta,
                              covariate = DATA0$pscore, x = hcv$x0[ix], y = 0.8, bw = hcv$h[ix] )
        surv1_24[ix] <- Beran(time = DATA1$obse_time, status = DATA1$delta,
                              covariate = DATA1$pscore, x = hcv$x0[ix], y = 2.4, bw = hcv$h[ix] )
        surv0_24[ix] <- Beran(time = DATA0$obse_time, status = DATA0$delta,
                              covariate = DATA0$pscore, x = hcv$x0[ix], y = 2.4, bw = hcv$h[ix] )
      }
      # Estimate the SPCE at the time 0.8 and 2.4
      spce1[j] <- mean(surv1_08) - mean(surv0_08)
      spce2[j] <- mean(surv1_24) - mean(surv0_24)
      
      # Estimate the cure probability
      cure1_x <- probcure(pscore, obse_time, delta, DATA1, x0 = hncure$x0, h = hncure$h)$q
      cure0_x <- probcure(pscore, obse_time, delta, DATA0, x0 = hncure$x0, h = hncure$h)$q
      
      # Estimate the SPCE about the cure rate
      spce_cure[j] <- mean(cure1_x) - mean(cure0_x)
    }
    
    ACE <- mean(ace)
    SPCE1 <- mean(spce1)
    SPCE2 <- mean(spce2)
    SPCE_cure <- mean(spce_cure)
    
    return(list(ACE_est = ACE, SPCE1_est = SPCE1, SPCE2_est = SPCE2, SPCE_cure_est = SPCE_cure))
  }
  
  # Calculate the bootstrap sample standard error
  PS_Beran_boot <- function(data){
    B <- 50
    Boot_spce_cure <- Boot_spce1 <- Boot_spce2 <- Boot_ace <- rep(0,B)
    b <- 1
    jj <- 0
    #estimate PS
    ps.model <- glm(z ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10, data = data,
                    family = binomial(link = "logit"), na.action = na.pass  )   
    Xdat <- predict(ps.model, newdata = data, type = "response")  
    Zdat <- cbind(Xdat)
    nZ <- nrow(Zdat)
    Xord <- sort(as.numeric(Xdat)) 
    
    # computes the pilot bandwidth g
    gpilot <- numeric(nn)
    for (ix in 1:nn)
    { 
      x_fix <- Xdat[ix]
      kneighbour <- round(length(Xdat)/4)
      gpilot[ix] <- gpilotfun(Xdata = Xdat, Xord = Xord, x_fix = x_fix, 
                              side = side, kneighbour = kneighbour)
    }
    
    NWweight <- function(Z, hn, lam = rep(1, nZ)){
      Z0 <- Z <- as.data.frame(Z)
      p1 <- dim(Z)[2]
      
      MZ0 <- array(apply(Z0, 2, function(x) rep(x, nZ)), dim = c(nZ, nZ, p1))
      MZ <- array(apply(Z, 2, function(x) rep(x, nZ)), dim = c(nZ, nZ, p1))  
      # each fixed p1, dim is n*n in which each column denotes a person

      tMZ <- aperm(MZ, c(2, 1, 3))   # each fixed p1, dim is n*n in which each row denotes a person
      dMZ <- MZ0 - tMZ # each fixed p1, dim is n*n, the i row is the i person, the k column detnotes a person's Z minus the k person's Z, namely Z-Z.k 
      
      Kp <- array(NA, dim = c(nZ, nZ, p1))
      for(i in 1:p1){
        
        ## Kp[,,i] denotes for fixed Z.i, in the the n*n matrix, the k column denotes a person's Z minus the k person's Z
        if(sum(dMZ[, , i] == 1 | dMZ[, , i] == -1 | dMZ[, , i] == 0) == (nZ * nZ) & any(dMZ[, , i] != 0)){
          Kp[, , i] <- Kernal_b(dMZ[, , i]) # if Z is a binary covariate, use the indicator function
        }
        else{
          Kp[, , i] <- Kernal_c(dMZ[, , i], hn) # if Z is a continuous covariate, use a kernel function
        }
      }
      Lam <- matrix(lam, nrow = nZ, ncol = nZ, byrow = T)
      prodKp <- apply(Kp, c(1, 2), prod)
      lamKp <- prodKp * Lam  

      return(lamKp)  
    }
    
    # calculate the NW weight
    w_matrix <- NWweight(Z = Zdat, hn = gpilot, lam = rep(1, nZ))
    repeat{
      databoot <- as.data.frame(matrix(0, nn, ncol(data)))
      databoot[, -c(1:2)] <- data[, -c(1:2)] 
      
      for (ij in 1:nn) {
        databoot[ij, c(1,2)] <- data[sample(c(1:nn), size = 1, prob = w_matrix[ij, ]), c(1,2)]
      }
      colnames(databoot) <- c("obse_time", "delta", "x1", "x2","x3", "x4","x5", "x6",
                              "x7", "x8","x9", "x10","z")
      
      Boot_result <- PS_Beran(databoot)
      Boot_ace[b] <- Boot_result$ACE_est
      Boot_spce1[b] <- Boot_result$SPCE1_est
      Boot_spce2[b] <- Boot_result$SPCE2_est
      Boot_spce_cure[b] <- Boot_result$SPCE_cure_est
      
      b <- b + 1
      if (b > B)
        break
    }
    # calculate the sample standard error
    Boot_ace_se <- sd(Boot_ace)
    Boot_spce1_se <- sd(Boot_spce1)
    Boot_spce2_se <- sd (Boot_spce2)
    Boot_spcecure_se <- sd(Boot_spce_cure)
    
    return(list(Boot_ace_se = Boot_ace_se, Boot_spce1_se = Boot_spce1_se,
                Boot_spce2_se = Boot_spce2_se, Boot_spcecure_se = Boot_spcecure_se))
  }
  
  summ_est <- PS_Beran(data_ran)
  boot_se <- PS_Beran_boot(data_ran)
  return(list(summ_est = summ_est, boot_se = boot_se))
}

## ----------- SET -------------------------------------------------------------
J <- 5    # number of the stratification
u1 <- 5   # upper bound of censored time, Tau_C
tau <- u1 - 0.5
nn <- 150         # sample size n=150, 300, 500
repli <- 500      # the number of replications

mu1 <- 2     # used in the baseline hazard function under treat z=1
mu0 <- 6     # used in the baseline hazard function under treat z=0
lamb <- 0.5  # used in the baseline hazard function
p <- 4       # the number of results obtained

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

summ_true <- summ_weak_alpha0 # true value

## ----------- Parallel --------------------------------------------------------
no_cores <- detectCores(logical = FALSE) 
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

MyRsult <- foreach(ii = 1:repli, .combine = rbind, .multicombine = TRUE,
                   .packages = packagelist)  %dopar%  Getace(ii)

MySummary1 <- Getsummary()
stopImplicitCluster()


time1 <- Sys.time()
runningave <- time1 - timestart

#save(MyRsult, file = "simulation_study/bootstrap_results/boots_weak_alp0_PE20_PSBeran.RData")

