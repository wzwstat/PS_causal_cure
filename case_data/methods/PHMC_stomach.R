###############################################################################
#
#  Title:         Estimating causal effects in observational studies for survival 
#                 data with a cure fraction using propensity score adjustment
#
#  Script:        PHMC_stomach.R
#
#  Authors:       Ziwen Wang | Chenguang Wang | Xiaoguang Wang
#
#  Purpose:       Reproducibility of the results of the PHMC method 
#                 is shown in Table 7 in the paper 
#
#  Note:          The script may take ~10 mins to execute. 
#
###############################################################################
# Necessary packages 
library(tidyverse)
library(smcure)

# Source functions and load data
mydata <- read.csv("case_data/data/stomach.csv")
source("case_data/function/function_stomach.R")
load("case_data/intermediate_results/bandwidth_Beran.Rdata")

##===== Cleaning data ====================================================
mydata <- mydata %>%
  dplyr::select(time, status, age, sex, size, surg) %>%
  mutate(status = case_when(status == "Alive" ~ 0,
                            status == "Dead" ~ 1)) %>%
  mutate(surg = case_when(surg == "Surgery performed" ~ 1,
                          surg == "Not recommended" ~ 0)) %>%
  mutate(sex = case_when(sex == "Male" ~ 1,
                         sex == "Female" ~ 0))

data <- mydata
data$age <- as.numeric(scale(data$age))
data$size <- as.numeric(scale(data$size))

##====== Fit PS_KM ============================================
j <- J <- 1
tau <- max(data$time)

PHMC_est <- function(data){

  DATA <- datastratify_stomach(data)[[j]]
  phmc <- smcure( Surv(time, status) ~ age + sex + size + surg,
                  cureform = ~ age + sex + size + surg,
                  data = DATA, model = "ph", Var = FALSE)

  #---Estimate the survival probability under treat Z=1 ------------------
  Z1 <- cbind(DATA$age, DATA$sex, DATA$size, 1)
  pred1 <- predictsmcure(phmc, newX = Z1, newZ = Z1, model = "ph")
  
  id <- dim(pred1$prediction)[2]
  t1 <- as.numeric(pred1$prediction[, id])   
  surv1_t <- pred1$prediction[, -id]
  Cox_S1_t <- as.numeric(apply(surv1_t, 1, mean))
  
  #---Estimate the survival probability under treat Z=0 ------------------
  Z0 <- cbind(DATA$age,  DATA$sex,  DATA$size,  0)
  pred0 <- predictsmcure(phmc, newX = Z0, newZ = Z0, model = "ph")
  
  id <- dim(pred0$prediction)[2]
  t0 <- as.numeric(pred0$prediction[, id])   
  surv0_t <- pred0$prediction[, -id]
  Cox_S0_t <- as.numeric(apply(surv0_t, 1, mean))
  
  #---Caculate the area under the survival curves ------------------
  surv1_area <- trapezoid(t1, Cox_S1_t)
  surv0_area <- trapezoid(t0, Cox_S0_t)
  
  #---caculate the RACE and SPCE for cure rate---------------------------------
  ace <- surv1_area - surv0_area
  spce <- min(Cox_S1_t[Cox_S1_t != 0]) - min(Cox_S0_t[Cox_S0_t != 0])
  
  return(list(ACE_est = ace, SPCE_cure_est = spce))
}

PHMC_boot <- function(data){
  B <- 50
  Boot_spce_cure <- Boot_ace <- rep(0, B)
  b <- 1
  jj <- 0
  nn <- nrow(data)
  ps_model <- glm(surg ~ age + size + sex, data = data, 
                  family = binomial(link="logit"), na.action = na.pass)    #estimate PS
  
  Xdat <- predict(ps_model, newdata = data, type = "response")  
  Zdat <- cbind(Xdat) 
  nZ <- nrow(Zdat)
  Xord <- sort(as.numeric(Xdat))
  gpilot <- numeric(nn)
  for(ix in 1:nn)
  { 
    x_fix <- Xdat[ix]
    kneighbour <- round(length(Xdat)/4)
    gpilot[ix] <- gpilotfun(Xdata = Xdat, Xord = Xord, x_fix = x_fix, 
                            side = side, kneighbour = kneighbour)
  }
  
  NWweight <- function(Z,hn,lam=rep(1,nZ)){
    Z0 <- Z <- as.data.frame(Z)
    p1 <- dim(Z)[2]
    MZ0 <- array(apply(Z0, 2, function(x) rep(x,nZ)), dim = c(nZ,nZ,p1))  # each fixed p1(Z is two dims, p1=2), dim is n*n in which each column denotes a person
    MZ <- array(apply(Z, 2, function(x) rep(x,nZ)), dim = c(nZ,nZ,p1))  # each fixed p1(Z is two dims, p1=2), dim is n*n in which each column denotes a person
    
    tMZ <- aperm(MZ, c(2,1,3))   # each fixed p1(Z is two dims, p1=2), dim is n*n in which each row denotes a person
    dMZ <- MZ0 - tMZ # each fixed p1(Z is two dims, p1=2), dim is n*n, the i row is the i person, the k column detnotes a person's Z minus the k person's Z, namely Z-Z.k 
    
    Kp <- array(NA,dim=c(nZ,nZ,p1))
    for(i in 1:p1){
       if(sum(dMZ[,,i]==1|dMZ[,,i]==-1|dMZ[,,i]==0)==(nZ*nZ) & any(dMZ[,,i]!=0)){
        Kp[,,i] <- Kernal_b(dMZ[,,i]) # if Z is a binay covariate, use the Indicator function
      }
      else{
        Kp[,,i] <- Kernal_c(dMZ[,,i],hn) # if Z is a continous covariate, use a kernel function
      }
    }
    Lam <- matrix(lam, nrow = nZ, ncol = nZ, byrow = T)
    prodKp <- apply(Kp, c(1, 2), prod)  # K(x) = K(x1) * K(X2), since Z is two dims. K(x) is n*n dim.
    lamKp <- prodKp * Lam  # each Z.k, times a disturbe value !!!
    
    return(lamKp)  
  }
  w_matrix <- NWweight(Z = Zdat, hn = gpilot, lam = rep(1, nZ))
  repeat {
    set.seed(5 * b)
    databoot <- as.data.frame(matrix(0, nn, ncol(data)))
    databoot[, -c(1:2)] <- data[, -c(1:2)]
    
    for (ij in 1:nn) {
      databoot[ij, c(1, 2)] <- data[sample(c(1:nn), size = 1, prob = w_matrix[ij,]), c(1, 2)]
    }
    colnames(databoot) <- c("time", "status", "age", "sex", "size", "surg")
    
    #data<-databoot
    PHMC_BOOT <- PHMC_est(databoot)
    Boot_ace[b] <- PHMC_BOOT$ACE_est
    Boot_spce_cure[b] <- PHMC_BOOT$SPCE_est
    
    b <- b + 1
    if (b > B)
      break
  }
  
  Boot_ace_se <- sd(Boot_ace)
  Boot_spcecure_se <- sd(Boot_spce_cure)
  
  return(list(Boot_ace_se = Boot_ace_se, Boot_spcecure_se = Boot_spcecure_se))
}

summ_est <- PHMC_est(data)
boot_se <- PHMC_boot(data)

Summary_table <- Get_table()

#save.image("case_data/intermediate_results/boots_PHMC_stomach.RData")
