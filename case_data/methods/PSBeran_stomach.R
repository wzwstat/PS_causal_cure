###############################################################################
#
#  Title:         Estimating causal effects in observational studies for survival 
#                 data with a cure fraction using propensity score adjustment
#
#  Script:        PSBeran_stomach.R
#
#  Authors:       Ziwen Wang | Chenguang Wang | Xiaoguang Wang
#
#  Purpose:       The reproducibility of the results of the PS-Beran method 
#                 is shown in Table 7 in the paper.

#  Note:          The script may take ~25 mins to execute. 
# 
###############################################################################
# Necessary packages 
library(tidyverse)
library(npcure)
library(condSURV)

# Source functions and load data
mydata <- read.csv("case_data/data/stomach.csv")
source("case_data/function/function_stomach.R")
load("case_data/intermediate_results/bandwidth.Rdata")

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

##====== Fit PS_Beran ============================================
J <- 5
tau <- max(data$time)

PS_Beran <- function(data){
  ace <- spce_cure <- rep(0, J)
  tgrid <- seq(min(data$time), max(data$time), length = 101)[-1]
  
  for(j in 1:J){
    
    DATA <- datastratify_stomach(data)[[j]]
    DATA1 <- DATA[DATA$surg == 1,]
    DATA0 <- DATA[DATA$surg == 0,]
    #hcv <- berancv(pscore, time, status, DATA, x0=DATA$pscore)
    #hncure <- probcurehboot(pscore, time, status, DATA, x0=DATA$pscore)
    hcv <- bandwid_hcv[[j]]
    hncure <- bandwid_cure[[j]]
    
    surv1_t <- surv0_t <- matrix(0, nrow(DATA), length(tgrid)) #hat S(t|Xi)
    for(ix in 1:nrow(DATA)){
      for(it in 1:length(tgrid)){
        surv1_t[ix,it]<- Beran(time = DATA1$time, status = DATA1$status,
                               covariate = DATA1$pscore, x = hcv$x0[ix], y =tgrid[it], bw = hcv$h[ix] )
        surv0_t[ix,it]<- Beran(time = DATA0$time, status = DATA0$status,
                               covariate = DATA0$pscore, x = hcv$x0[ix], y =tgrid[it], bw = hcv$h[ix] )
      }
    }
    
    surv1_area <- surv0_area <- rep(0,length(tgrid))
    for(ix in 1:nrow(DATA) ){
      surv1_area[ix] <- trapezoid(tgrid, surv1_t[ix,])
      surv0_area[ix] <- trapezoid(tgrid, surv0_t[ix,])
    }
    ace[j] <- mean(surv1_area - surv0_area)

    cure1_x <- probcure(pscore, time, status, DATA1, x0 = hncure$x0, h=hncure$h)$q
    cure0_x <- probcure(pscore, time, status, DATA0, x0 = hncure$x0, h=hncure$h)$q
    spce_cure[j] <- mean(cure1_x) - mean(cure0_x)
  }
  
  ACE <- mean(ace)
  SPCE_cure <- mean(spce_cure)
  
  return(list(ACE_est = ACE, SPCE_cure_est = SPCE_cure))
}

PSBeran_boot<-function(data){
  B <- 50
  Boot_spce_cure <- Boot_ace <- rep(0, B)
  b <- 1
  jj <- 0
  nn <- nrow(data)
  ps_model <- glm(surg ~ age + size +sex, data = data, 
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
  
  NWweight <- function(Z, hn, lam = rep(1, nZ)){
    Z0 <- Z <- as.data.frame(Z)
    p1 <- dim(Z)[2]
    MZ0 <- array(apply(Z0, 2, function(x) rep(x,nZ)), dim = c(nZ,nZ,p1))  # each fixed p1(Z is two dims, p1=2), dim is n*n in which each column denotes a person
    MZ <- array(apply(Z, 2, function(x) rep(x,nZ)), dim = c(nZ,nZ,p1))  # each fixed p1(Z is two dims, p1=2), dim is n*n in which each column denotes a person
    
    tMZ <- aperm(MZ, c(2, 1, 3))   # each fixed p1(Z is two dims, p1=2), dim is n*n in which each row denotes a person
    dMZ <- MZ0 - tMZ # each fixed p1(Z is two dims, p1=2), dim is n*n, the i row is the i person, the k column detnotes a person's Z minus the k person's Z, namely Z-Z.k 
    
    Kp <- array(NA, dim = c(nZ, nZ, p1))
    for(i in 1:p1){
      if(sum(dMZ[,,i]==1|dMZ[,,i]==-1|dMZ[,,i]==0)==(nZ*nZ) & any(dMZ[,,i]!=0)){
        Kp[,,i] <- Kernal_b(dMZ[,,i]) # if Z is a binay covariate, use the Indicator function
      }
      else{
        Kp[,,i] <- Kernal_c(dMZ[,,i], hn) # if Z is a continous covariate, use a kernel function
      }
    }
    Lam <- matrix(lam, nrow = nZ, ncol = nZ, byrow = T)  
    prodKp <- apply(Kp, c(1,2), prod)  # K(x) = K(x1) * K(X2), since Z is two dims. K(x) is n*n dim.
    lamKp <- prodKp * Lam  # each Z.k, times a disturbe value !!!
    return(lamKp)  
  }
  w_matrix <- NWweight(Z = Zdat, hn = gpilot, lam = rep(1, nZ))
  repeat{
    set.seed(5 * b)
    databoot <- as.data.frame(matrix(0, nn, ncol(data)))
    databoot[,-c(1:2)] <- data[,-c(1:2)]
    
    for (ij in 1:nn) {
      databoot[ij, c(1,2)] <- data[sample(c(1:nn), size = 1, prob = w_matrix[ij, ]), c(1,2)]
    }
    colnames(databoot) <- c("time", "status", "age", "sex","size", "surg")
    
    PSBeran_BOOT <- PS_Beran(databoot)
    Boot_ace[b] <- PSBeran_BOOT$ACE_est
    Boot_spce_cure[b] <- PSBeran_BOOT$SPCE_cure_est
    
    b <- b + 1
    if (b > B)
      break
    
  }
  Boot_ace_se <- sd(Boot_ace)
  Boot_spcecure_se <- sd(Boot_spce_cure)
  
  return(list(Boot_ace_se = Boot_ace_se, Boot_spcecure_se = Boot_spcecure_se))
}

summ_est <- PS_Beran(data)
boot_se <- PSBeran_boot(data)

Summary_table <- Get_table()

# save.image("case_data/intermediate_results/boots_PSBeran_stomach.RData")
