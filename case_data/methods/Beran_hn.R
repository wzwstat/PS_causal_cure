##==================================================================================
## Calculate and Save the bandwidth of Beran
##
## This is an intermediate process to get bandwidth of Beran estimator, 
## the result obtained have been saved in the folder:
## "case_data/intermediate_results/bandwidth_Beran"
##==================================================================================
# Necessary packages 
library(tidyverse)
library(npcure)
library(condSURV)

# Source functions and load data
source("case_data/function/function_stomach.R")
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

##====== Calculate the bandwidth of Beran estimator============================================
J <- 1
tau <- max(data$time)

DATA <- datastratify_stomach(data)[[J]]

hcv1 <- berancv(pscore, time, status, DATA1, x0 = DATA1$pscore)

# save(hcv1, file = "case_data/intermediate_results/bandwidth_Beran.Rdata")
