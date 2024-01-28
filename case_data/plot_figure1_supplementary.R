#################################################################################
#
#  Title:         Estimating causal effects in observational studies for survival 
#                 data with a cure fraction using propensity score adjustment
#
#  Script:        plot_figure1_supplementary.R
#
#  Authors:       Ziwen Wang | Chenguang Wang | Xiaoguang Wang
#
#  Purpose:       Create Figure 1 of Supplementary material 
#
#                 Automatically saved in the folder "case_data/results" 
# 
#################################################################################

# Necessary packages -----------------------------------------------------
library(tidyverse)
library(npcure)
library(condSURV)

# Source functions and load data------------------------------------------
source("case_data/function/function_stomach.R")
load("case_data/intermediate_results/bandwidth.Rdata")
mydata <- read.csv("case_data/data/stomach.csv")

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

##====== First Fit PS_Beran ==============================================

J <- 5
tgrid <- seq(min(data$time), max(data$time), length = 101)[-1]
tau <- max(data$time)

PS_stratify <- function(data){
  surv1_mean <- surv0_mean <- matrix(0, J, 100)

  for(j in 1:J){
    
    DATA <- datastratify_stomach(data)[[j]]
    DATA1 <- DATA[DATA$surg == 1,]
    DATA0 <- DATA[DATA$surg == 0,]
    
    hcv <- bandwid_hcv[[j]]
    hncure <- bandwid_cure[[j]]
    
    surv1_t <- surv0_t <- matrix(0, nrow(DATA), length(tgrid)) #hat S(t|Xi)
    for(ix in 1:nrow(DATA)){
      for(it in 1:length(tgrid)){
        surv1_t[ix,it] <- Beran(time = DATA1$time, status = DATA1$status,
                               covariate = DATA1$pscore, x = hcv$x0[ix], y =tgrid[it], bw = hcv$h[ix] )
        surv0_t[ix,it] <- Beran(time = DATA0$time, status = DATA0$status,
                               covariate = DATA0$pscore, x = hcv$x0[ix], y =tgrid[it], bw = hcv$h[ix] )
      }
    }
    surv1_mean[j, ] <- apply(surv1_t, 2, mean)
    surv0_mean[j, ] <- apply(surv0_t, 2, mean)
  }
  
  return(list(surv1_mean = surv1_mean, surv0_mean = surv0_mean))
}

surv_fit <- PS_stratify(data)

##------------------------------------------------
## Then Plot Figure 1 in Supplementary
##------------------------------------------------

postscript("case_data/results/Figure1_supplementary.eps", width = 14, height = 6, horizontal = FALSE) # save as .eps format
par(mfrow = c(2, 3))

plot(
  surv_fit$surv1_mean[1, ] ~ tgrid,
  type = 's',
  lwd = 2,
  ylim = c(0.3, 1),
  xlim = c(0, tau),
  xlab = "Months",
  ylab = "Survival probability",
  main = "Stratum 1"
)
points(surv_fit$surv0_mean[1, ],
       type = 's',
       lwd = 2,
       lty = 2)
legend(
  "bottomleft",
  c("Surgery", "Non-surgery", "RACE: 11.07", "SPCE:  0.21"),
  lty = c(1, 2, 0, 0),
  lwd = c(2),
  cex = 1
)

plot(
  surv_fit$surv1_mean[2, ] ~ tgrid,
  type = 's',
  lwd = 2,
  ylim = c(0.3, 1),
  xlim = c(0, tau),
  xlab = "Months",
  ylab = "Survival probability",
  main = "Stratum 2"
)
points(surv_fit$surv0_mean[2, ],
       type = 's',
       lwd = 2,
       lty = 2)
legend(
  "bottomleft",
  c("Surgery", "Non-surgery",  "RACE: 11.17", "SPCE:  0.22"),
  lty = c(1, 2, 0, 0),
  lwd = c(2),
  cex = 1
)

plot(
  surv_fit$surv1_mean[3, ] ~ tgrid,
  type = 's',
  lwd = 2,
  ylim = c(0.3, 1),
  xlim = c(0, tau),
  xlab = "Months",
  ylab = "Survival probability",
  main = "Stratum 3"
)
points(surv_fit$surv0_mean[3, ],
       type = 's',
       lwd = 2,
       lty = 2)
legend(
  "bottomleft",
  c("Surgery", "Non-surgery", "RACE: 12.21", "SPCE:  0.34"),
  lty = c(1, 2, 0, 0),
  lwd = c(2),
  cex = 1
)

plot(
  surv_fit$surv1_mean[4, ] ~ tgrid,
  type = 's',
  lwd = 2,
  ylim = c(0.3, 1),
  xlim = c(0, tau),
  xlab = "Months",
  ylab = "Survival probability",
  main = "Stratum 4"
)
points(surv_fit$surv0_mean[4, ],
       type = 's',
       lwd = 2,
       lty = 2)
legend(
  "bottomleft",
  c("Surgery", "Non-surgery", "RACE: 7.78", "SPCE: 0.10"),
  lty = c(1, 2, 0, 0),
  lwd = c(2),
  cex = 1
)


plot(
  surv_fit$surv1_mean[5, ] ~ tgrid,
  type = 's',
  lwd = 2,
  ylim = c(0.3, 1),
  xlim = c(0, tau),
  xlab = "Months",
  ylab = "Survival probability",
  main = "Stratum 5"
)
points(surv_fit$surv0_mean[5, ],
       type = 's',
       lwd = 2,
       lty = 2)
legend(
  "bottomleft",
  c("Surgery", "Non-surgery", "RACE: 7.16", "SPCE: 0.13"),
  lty = c(1, 2, 0, 0),
  lwd = c(2),
  cex = 1
)

dev.off()


