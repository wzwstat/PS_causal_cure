#################################################################################
#
#  Title:         Estimating causal effects in observational studies for survival 
#                 data with a cure fraction using propensity score adjustment
#
#  Script:        Plot_methods_figure3.R
#
#  Authors:       Ziwen Wang | Chenguang Wang | Xiaoguang Wang
#
#  Purpose:       Create Figure 3 of the paper
#
#                 Automatically saved in the folder "case_data/results"
# 
#################################################################################

# Necessary packages -----------------------------------------------------
library(tidyverse)
library(condSURV)
library(smcure)

# Source functions and load data -----------------------------------------
source("case_data/function/function_stomach.R")
load("case_data/intermediate_results/bandwidth_Beran.Rdata")
load("case_data/intermediate_results/bandwidth.Rdata")
mydata <- read.csv("case_data/data/stomach.csv")

## ===== Cleaning data ===================================================
mydata <- mydata %>%
  dplyr::select(time, status, age, sex, size, surg) %>%
  mutate(status = case_when(status == "Alive" ~ 0,
                            status == "Dead" ~ 1)) %>%
  mutate(surg = case_when(surg == "Surgery performed" ~ 1,
                          surg == "Not recommended" ~ 0))  %>%
  mutate(sex = case_when(sex == "Male" ~ 1,
                         sex == "Female" ~ 0))

data <- mydata
data$age <- as.numeric(scale(data$age))
data$size <- as.numeric(scale(data$size))

## ===== PS_Beran method =================================================

PS_Beran_plot <- function(data){
  SUR0 <- SUR1 <- matrix(0, nrow = J,ncol = 100)
  
  for(j in 1:J){
    DATA <- datastratify_stomach(data)[[j]]
    DATA1 <- DATA[DATA$surg == 1,]
    DATA0 <- DATA[DATA$surg == 0,]
    
    hcv <- bandwid_hcv[[j]]
    hncure <- bandwid_cure[[j]]
    
    surv1_t <- surv0_t <- matrix(0, nrow(DATA), length(tgrid)) #hat S(t|Xi)
    for(ix in 1:nrow(DATA)){
      for(it in 1:length(tgrid)){
        surv1_t[ix, it] <- Beran(time = DATA1$time, status = DATA1$status,
                               covariate = DATA1$pscore, x = hcv$x0[ix], y = tgrid[it], bw = hcv$h[ix] )
        surv0_t[ix, it] <- Beran(time = DATA0$time, status = DATA0$status,
                               covariate = DATA0$pscore, x = hcv$x0[ix], y = tgrid[it], bw = hcv$h[ix] )
      }
    }
    SUR1[j, ] <- apply(surv1_t, 2, mean)
    SUR0[j, ] <- apply(surv0_t, 2, mean)
  }
  
  surv1_mean <- apply(SUR1, 2, mean)
  surv0_mean <- apply(SUR0, 2, mean)
  
  return(list(surv1_mean = surv1_mean, surv0_mean = surv0_mean))
}

J <- 5
tgrid <- seq(min(data$time), max(data$time), length = 101)[-1]
tau <- max(data$time)
pre_PSberan <- PS_Beran_plot(data)

pre_S1_psberan <- pre_PSberan$surv1_mean
pre_S0_psberan <- pre_PSberan$surv0_mean


## ===== Beran method ====================================================

Beran_plot <- function(data){
  SUR0 <- SUR1 <- matrix(0, nrow = J, ncol = 100)
  
  for(j in 1:J){
    DATA <- datastratify_stomach(data)[[j]]
    DATA1 <- DATA[DATA$surg == 1,]
    DATA0 <- DATA[DATA$surg == 0,]
    
    hcv <- hncure <- hcv1
    
    surv1_t <- surv0_t <- matrix(0, nrow(DATA), length(tgrid)) #hat S(t|Xi)
    for(ix in 1:nrow(DATA)){
      for(it in 1:length(tgrid)){
        surv1_t[ix, it] <- Beran(time = DATA1$time, status = DATA1$status,
                               covariate = DATA1$pscore, x = hcv1$x0[ix], y = tgrid[it], bw = hcv1$h[ix] )
        surv0_t[ix, it] <- Beran(time = DATA0$time, status = DATA0$status,
                               covariate = DATA0$pscore, x = hcv1$x0[ix], y = tgrid[it], bw = hcv1$h[ix] )
        }
    }
    SUR1[j, ] <- apply(surv1_t, 2, mean)
    SUR0[j, ] <- apply(surv0_t, 2, mean)
  }
  
  surv1_mean <- apply(SUR1, 2, mean)
  surv0_mean <- apply(SUR0, 2, mean)
  
  return(list(surv1_mean = surv1_mean, surv0_mean = surv0_mean))
}

J <- 1
tgrid <- seq( min(data$time), max(data$time), length=101)[-1]
tau <- max(data$time)
pre_Beran <- Beran_plot(data)

pre_S1_beran <- pre_Beran$surv1_mean
pre_S0_beran <- pre_Beran$surv0_mean

## ===== PHMC method =====================================================

PHMC_plot <- function(data){
  
  DATA <- datastratify_stomach(data)[[j]]
    
  phmc <- smcure( Surv(time, status) ~ age + as.factor(sex) + size + surg,
                  cureform = ~ age + as.factor(sex) + size + surg,
                  data = DATA, model = "ph", Var = FALSE)

  #--- Estimate the survival probability under treat Z=1 -----------------
  Z1 <- cbind(DATA$age, DATA$sex, DATA$size, 1)
  pred1 <- predictsmcure(phmc, newX = Z1, newZ = Z1, model = "ph")
  
  id <- dim(pred1$prediction)[2]
  t1 <- as.numeric(pred1$prediction[, id])   
  surv1_t <- pred1$prediction[, -id]
  Cox_S1_t <- as.numeric(apply(surv1_t, 1, mean))
  
  #--- Estimate the survival probability under treat Z=0 -----------------
  Z0 <- cbind(DATA$age,  DATA$sex,  DATA$size,  0)
  pred0 <- predictsmcure(phmc, newX = Z0, newZ = Z0, model = "ph")
  
  id <- dim(pred0$prediction)[2]
  t0 <- as.numeric(pred0$prediction[, id])   
  surv0_t <- pred0$prediction[, -id]
  Cox_S0_t <- as.numeric(apply(surv0_t, 1, mean))
    
  return(list(Cox_S1_t = Cox_S1_t, t1 = t1, Cox_S0_t = Cox_S0_t, t0 = t0))
}

j <- J <- 1
tau <- max(data$time)
pre_PHMC <- PHMC_plot(data)

pre_S1_PHMC <- pre_PHMC$Cox_S1_t
t1_PHMC <- pre_PHMC$t1
pre_S0_PHMC <- pre_PHMC$Cox_S0_t
t0_PHMC <- pre_PHMC$t0

## ===== PSKM method =====================================================

PSKM_plot <- function(data){
  
  SUR1_KM <- t1_sum <- matrix(0, nrow = J,ncol = 60)
  SUR0_KM <- t0_sum <- matrix(0, nrow = J,ncol = 38)
  for(j in 1:J){
    
    DATA <- datastratify_stomach(data)[[j]]
    DATA1 <- DATA[DATA$surg==1, ]
    DATA0 <- DATA[DATA$surg==0, ]
    
    kmfit1 <- survfit(Surv(DATA1$time, DATA1$status == 1) ~ 1)
    t1_sum[j, ] <- kmfit1$time
    SUR1_KM[j, ] <- kmfit1$surv
    
    kmfit0 <- survfit(Surv(DATA0$time, DATA0$status == 1) ~ 1)
    t0 <- kmfit0$time
    surv0_t <- kmfit0$surv
    
    t0 <- c(t0, seq(max(t0), max(data$time), length.out = 2))
    surv0_t <- c(surv0_t, rep(min(surv0_t[surv0_t != 0]), 2))
    t0_sum[j, ] <- t0[1:38]
    SUR0_KM[j, ] <- surv0_t[1:38]
  }
  
  surv1_KM <- apply(SUR1_KM, 2, mean)
  surv0_KM <- apply(SUR0_KM, 2, mean)
  t1 <- apply(t1_sum, 2, mean)
  t0 <- apply(t0_sum, 2, mean)
  
  return(list(surv1_KM = surv1_KM, t1=t1, surv0_KM = surv0_KM, t0 = t0))
}
J <- 5
tau <- max(data$time)
pre_PSKM <- PSKM_plot(data)

pre_S1_PSKM <- pre_PSKM$surv1_KM
pre_S0_PSKM <- pre_PSKM$surv0_KM
t1_PSKM <- pre_PSKM$t1
t0_PSKM <- pre_PSKM$t0

## =====================================================================
## plot Figure 3
## =====================================================================

#postscript("case_data/results/Figure3.eps") # save as .eps format
postscript("case_data/results/Figure3.eps", width = 12, height = 6) # save as .eps format
par(mfrow = c(1, 2))
## ===== for Non-surgery ===============================================

plot(pre_S1_psberan ~ tgrid, type = 's',ylim = c(0.5, 1), xlim = c(0, tau),
     lwd = 2, xlab = "Months", ylab = "Survival curves")
points(pre_S1_beran ~ tgrid, type = 's', lty = 2, lwd = 2)
points(pre_S1_PHMC ~ t1_PHMC, type = 's', lty = 3, lwd = 2)
points(pre_S1_PSKM ~ t1_PSKM, type = 's', lty = 4, lwd = 2)

legend("bottomleft", c("Surgery, PS-Beran",
                       "Surgery, Beran",
                       "Surgery, PHMC",
                       "Surgery, PS-KM"),
       lty = c(1, 2, 3, 4), lwd = c(2, 2, 2, 2),
       cex = 1)

## ===== for Non-surgery ================================================

plot(pre_S0_psberan ~ tgrid, type = 's', ylim = c(0.5, 1), xlim = c(0, tau),
     lwd = 2, xlab = "Months", ylab = "Survival curves")
points(pre_S0_beran ~ tgrid, type = 's', lty = 2, lwd = 2)
points(pre_S0_PHMC ~ t0_PHMC, type = 's', lty = 3, lwd = 2)
points(pre_S0_PSKM ~ t0_PSKM, type = 's', lty = 4, lwd = 2)

legend("bottomleft", c("Non-surgery, PS-Beran",
                       "Non-surgery, Beran",
                       "Non-surgery, PHMC",
                       "Non-surgery, PS-KM"),
       lty = c(1, 2, 3, 4), lwd = c(2, 2, 2, 2), 
       cex = 0.8 )

dev.off()

