###############################################################################
#
#  Title:         Estimating causal effects in observational studies for survival 
#                 data with a cure fraction using propensity score adjustment
#
#  Script:        PSKM_stomach.R
#
#  Authors:       Ziwen Wang | Chenguang Wang | Xiaoguang Wang
#
#  Purpose:       Reproducibility of the results of the PSKM method 
#                 is shown in Table 7 in the paper 
#
#  Note:          The script may take ~2 mins to execute.  
# 
###############################################################################
# Necessary packages 
library(tidyverse)

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
                          surg == "Not recommended" ~ 0))

data <- mydata
data$age <- as.numeric(scale(data$age))
data$size <- as.numeric(scale(data$size))

##====== Fit PS_KM ============================================
J <- 5
tau <- max(data$time)

PS_KM <- function(data){
  
  spce <- ace <- rep(0,J)
  for(j in 1:J){
    DATA <- datastratify_stomach(data)[[j]]
    DATA1 <- DATA[DATA$surg == 1,]
    DATA0 <- DATA[DATA$surg == 0,]
    
    kmfit1 <- survfit(Surv(DATA1$time, DATA1$status==1)~ 1)
    t1 <- kmfit1$time
    surv1_t <- kmfit1$surv
    
    kmfit0 <- survfit(Surv(DATA0$time, DATA0$status==1)~ 1)
    t0 <- kmfit0$time
    surv0_t <- kmfit0$surv
    
    t0 <- c(t0, seq(max(t0), max(data$time) , length.out = 2))
    su0 <- rep(min(surv0_t[surv0_t != 0]), 2)
    surv0_t <- c(surv0_t, su0)
    
    surv1_area <- trapezoid(t1, surv1_t)
    surv0_area <- trapezoid(t0, surv0_t)
    
    ace[j] <- surv1_area - surv0_area
    spce[j] <- min(surv1_t[surv1_t != 0]) - min(surv0_t[surv0_t != 0])
  }
  
  ACE <- mean(ace)
  SPCE <- mean(spce)
  
  return(list(ACE_est = ACE, SPCE_cure_est = SPCE))
}

PSKM_boot <- function(data){
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
    MZ0 <- array(apply(Z0, 2, function(x) rep(x, nZ)), dim = c(nZ, nZ, p1))  # each fixed p1(Z is two dims, p1=2), dim is n*n in which each column denotes a person
    MZ <- array(apply(Z, 2, function(x) rep(x, nZ)), dim = c(nZ, nZ, p1))  # each fixed p1(Z is two dims, p1=2), dim is n*n in which each column denotes a person
    
    tMZ <- aperm(MZ, c(2, 1, 3))   # each fixed p1(Z is two dims, p1=2), dim is n*n in which each row denotes a person
    dMZ <- MZ0 - tMZ # each fixed p1(Z is two dims, p1=2), dim is n*n, the i row is the i person, the k column detnotes a person's Z minus the k person's Z, namely Z-Z.k 
    
    Kp <- array(NA, dim = c(nZ, nZ, p1))
    for(i in 1:p1){
      ## Kp[,,i] denotes for fixed Z.i, in the the n*n matrix, the k column denotes a person's Z minus the k person's Z
      # Kp[,,i]=KK(dMZ[,,i],hn[i])  # here, i deonts Z.i, i=1,2,...,p1'dim
      if(sum(dMZ[,,i]==1|dMZ[,,i]==-1|dMZ[,,i]==0)==(nZ*nZ) & any(dMZ[,,i]!=0)){
        Kp[,,i] <- Kernal_b(dMZ[,,i]) # if Z is a binay covariate, use the Indicator function
      }
      else{
        Kp[,,i] <- Kernal_c(dMZ[,,i],hn) # if Z is a continous covariate, use a kernel function
      }
    }
    Lam <- matrix(lam, nrow = nZ, ncol = nZ, byrow = T)  
    prodKp <- apply(Kp, c(1,2), prod)  # K(x) = K(x1) * K(X2), since Z is two dims. K(x) is n*n dim.
    lamKp <- prodKp * Lam  # each Z.k, times a disturbe value !!!
    return(lamKp)  
  }
  w_matrix <- NWweight(Z = Zdat, hn = gpilot, lam = rep(1, nZ))
  repeat{
    set.seed(65 * b)
    databoot <- as.data.frame(matrix(0, nn, ncol(data)))
    databoot[, -c(1:2)] <- data[, -c(1:2)]
    
    for (ij in 1:nn) {
      databoot[ij, c(1, 2)] <- data[sample(c(1:nn), size = 1, prob = w_matrix[ij,]), c(1, 2)]
    }
    colnames(databoot) <- c("time", "status", "age", "sex", "size", "surg")
    
    PSKM_BOOT <- PS_KM(databoot)
    Boot_ace[b] <- PSKM_BOOT$ACE_est
    Boot_spce_cure[b] <- PSKM_BOOT$SPCE_cure_est
    
    b <- b + 1
    if (b > B)
      break
  }
  Boot_ace_se <- sd(Boot_ace)
  Boot_spcecure_se <- sd(Boot_spce_cure)
  
  return(list(Boot_ace_se = Boot_ace_se, Boot_spcecure_se = Boot_spcecure_se))
}

summ_est <- PS_KM(data)
boot_se <- PSKM_boot(data)

Summary_table <- Get_table()

#save.image("case_data/intermediate_results/boots_PSKM_stomach.RData")
