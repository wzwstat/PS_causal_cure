##==================================================================================
## Calculate and Save the bandwidth of Beran for each stratum
##
## This is an intermediate process to get bandwidth, 
## the result obtained have been saved in the folder:
## "case_data/intermediate_results/bandwidth"
##==================================================================================
# Necessary packages 
library(tidyverse)
library(npcure)
library(condSURV)

# Source functions and load data
mydata <- read.csv("case_data/data/stomach.csv")
source("case_data/function/function_stomach.R")

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

##====== Fit PS_Beran and calculate the bandwidth ============================================
J <- 5
tau <- max(data$time)

DATA1 <- datastratify_stomach(data)[[1]]
DATA2 <- datastratify_stomach(data)[[2]]
DATA3 <- datastratify_stomach(data)[[3]]
DATA4 <- datastratify_stomach(data)[[4]]
DATA5 <- datastratify_stomach(data)[[5]]

hcv1 <- berancv(pscore, time, status, DATA1, x0 = DATA1$pscore)
hcv2 <- berancv(pscore, time, status, DATA2, x0 = DATA2$pscore)
hcv3 <- berancv(pscore, time, status, DATA3, x0 = DATA3$pscore)
hcv4 <- berancv(pscore, time, status, DATA4, x0 = DATA4$pscore)
hcv5 <- berancv(pscore, time, status, DATA5, x0 = DATA5$pscore)
bandwid_hcv <- list(hcv1, hcv2, hcv3, hcv4, hcv5)

hncure1 <- probcurehboot(pscore, time, status, DATA1, x0 = DATA1$pscore)
hncure2 <- probcurehboot(pscore, time, status, DATA2, x0 = DATA2$pscore)
hncure3 <- probcurehboot(pscore, time, status, DATA3, x0 = DATA3$pscore)
hncure4 <- probcurehboot(pscore, time, status, DATA4, x0 = DATA4$pscore)
hncure5 <- probcurehboot(pscore, time, status, DATA5, x0 = DATA5$pscore)
bandwid_cure <- list(hncure1, hncure2, hncure3, hncure4, hncure5)

# save(bandwid_hcv, bandwid_cure, file = "case_data/intermediate_results/bandwidth.Rdata")

