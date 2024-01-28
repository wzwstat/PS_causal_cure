###############################################################################
#
#  Title:         Estimating causal effects in observational studies for survival 
#                 data with a cure fraction using propensity score adjustment
#
#  Script:        table5&6_figure1&2_table2supplementary.R
#
#  Authors:       Ziwen Wang | Chenguang Wang | Xiaoguang Wang
#
#  Purpose:       Create results of Table 5 (print in line 55), 
#                 Table 6 (print in line 107),
#                 Figure 1 (print in lines 118-122), 
#                 Figure 2 (print in lines 136-228),
#                 and Table 2 in Supplementary Material (print in line 264).
#
#                 Automatically saved in the folder "case_data/results"
# 
###############################################################################

# Necessary packages -----------------------------------------------------
library(tidyverse)

# Source functions and load data------------------------------------------
source("case_data/function/function_stomach.R")
mydata <- read.csv("case_data/data/stomach.csv")

##===== Cleaning data ====================================================
mydata <- mydata %>%
  dplyr::select(time, status, age, sex, size, surg) %>%
  mutate(status = case_when(status == "Alive" ~ 0,
                            status == "Dead" ~ 1)) %>%
  mutate(surg = case_when(surg == "Surgery performed" ~ 1,
                          surg == "Not recommended" ~ 0)) 

##=========================================================================
## Table 5 Summary statistics of the covariates for the stomach data 
##=========================================================================
age_sum <- mydata %>%
  group_by(surg) %>% 
  summarise(n = n(), age_mean = mean(age), 
            Q1 = quantile(age, 0.25), Q2 = quantile(age, 0.75) )

size_sum <- mydata %>%
  group_by(surg) %>%
  summarise(n = n(), size_mean = mean(size), 
            Q1 = quantile(size, 0.25), Q2 = quantile(size, 0.75) )
  
sex_sum <- mydata %>%
  group_by(surg, sex) %>%
  summarise(n = n(), .groups = "drop")

##----- Table5: surg=1 means receive surgery ----------------------------------
table5_sum <- cbind(round(age_sum, 1), round(size_sum, 1),
                    sex_sum[sex_sum$sex=="Male",], sex_sum[sex_sum$sex=="Female",])
print(table5_sum)


##===== Scale data =============================================================
data <- mydata
data$age <- as.numeric(scale(data$age))
data$size <- as.numeric(scale(data$size))

##===== Propensity score estimation ==========================================
library(ResourceSelection)
ps_model <- glm(surg ~ age + size + sex , data =  data, 
                family = binomial(link="logit"), na.action=na.pass)    #estimate PS
summary(ps_model)
hoslem.test(data$surg, fitted(ps_model), g = 10)


##=========================================================================
## Table 6 summary statistics for the PS 
##=========================================================================
J <- 5
DATA1 <- datastratify_stomach(data)[[1]]
DATA2 <- datastratify_stomach(data)[[2]]
DATA3 <- datastratify_stomach(data)[[3]]
DATA4 <- datastratify_stomach(data)[[4]]
DATA5 <- datastratify_stomach(data)[[5]]

table6 <- cbind(
  DATA1 %>%
    group_by(surg) %>% summarise(d1 = n()),
  DATA2 %>%
    group_by(surg) %>% summarise(d2 = n()),
  DATA3 %>%
    group_by(surg) %>% summarise(d3 = n()),
  DATA4 %>%
    group_by(surg) %>% summarise(d4 = n()),
  DATA5 %>%
    group_by(surg) %>% summarise(d5 = n())
)

pre_ps <- predict(ps_model, newdata = data, type = "response")
orderps <- sort(pre_ps)    #sort data based on PS
q <- quantile(orderps, probs = seq(0, 1, 1 / J))
q[1] <- 0
q[J+1] <- 1

##----- Table6 -------------------------------------------------------
table6 <- cbind(table6$d1, table6$d2, table6$d3, table6$d4, table6$d5)
rownames(table6) <- c("Non_surgery", "Surgery")
surg_count <- data %>% group_by(surg) %>% summarise(n = n()) 

table6_sum <- cbind(rbind(surg_count, sum(surg_count$n)),
                    rbind(table6, table(cut(orderps, breaks = q))) )
print(table6_sum)


##=========================================================================
## Figure 1 KM curves
##=========================================================================
tau.star <- max(data$time[which(data$status == 1)])

library(survival)
kmfit <- survfit(Surv(data$time, data$status) ~ data$surg)

postscript("case_data/results/Figure1.eps", width = 6, height = 6, horizontal = FALSE) # save as .eps format
plot(kmfit, xlim = c(0,max(data$time)), lty = 2:1, 
     ylab = "Survival curvrs", xlab = "Months", ylim = c(0.5,1))
legend("bottomleft",legend = c("Surgery","Non-Surgery"),lty = c(1,2), cex = 0.8)
dev.off()

##===== Maller-Zhou test ================================================
library(npcure)
testmz(time, status, data)

data1 <- data[data$surg == 1, ]
data0 <- data[data$surg == 0, ]
testmz(time, status, data1)
testmz(time, status, data0)

##=========================================================================
## Figure 2 Balance checking
##=========================================================================
postscript("case_data/results/Figure2.eps", width = 12, height = 7)
par(mfrow = c(3, 5))

#------- age --------------------------------------------------------------
plot(density(DATA1$age[DATA1$surg == 0]), ylab = "Age", 
     main = "Stratum 1", xlab=" ", xlim = c(-2, 2), ylim = c(0, 1.5) )
points(density(DATA1$age[DATA1$surg == 1]), type = 'l', lty = 2)
legend("topleft", legend = c("Non-surgery", "Surgery"), lty = c(1, 2), cex = 0.6)

plot(density(DATA2$age[DATA2$surg == 0]), ylab = " ", main = "Stratum 2",
     xlab = "", xlim = c(-2, 2), ylim = c(0, 1.5) )
points(density(DATA2$age[DATA2$surg == 1]), type = 'l', lty = 2)

plot(density(DATA3$age[DATA3$surg == 0]),ylab = "",main = "Stratum 3",
     xlab = "", xlim = c(-2, 2), ylim = c(0, 1.5) )
points(density(DATA3$age[DATA3$surg == 1]),type='l', lty = 2)

plot(density(DATA4$age[DATA4$surg == 0]),  ylab = "", main = "Stratum 4",
     xlab = "", xlim = c(-3, 2), ylim = c(0, 1.5) )
points(density(DATA4$age[DATA4$surg == 1]), type = 'l', lty = 2)

plot(density(DATA5$age[DATA5$surg == 0]),  ylab = "", main = "Stratum 5",
     xlab = "", xlim = c(-5, 1), ylim = c(0, 1.5) )
points(density(DATA5$age[DATA5$surg == 1]), type = 'l', lty = 2)

#------- size -------------------------------------------------------------
plot(density(DATA1$size[DATA1$surg == 0]), ylab = "Tumor size", main = "", 
     xlab = "",  xlim = c(-1.5, 4), ylim = c(0, 1.6) )
points(density(DATA1$size[DATA1$surg == 1]), type = 'l', lty = 2)
legend("topleft",legend = c("Non-surgery", "Surgery"), lty = c(1, 2), cex = 0.6)


plot(density(DATA2$size[DATA2$surg == 0]), ylab = "", main = "",
     xlab = "", xlim = c(-1.5, 3), ylim = c(0, 1.6) )
points(density(DATA2$size[DATA2$surg == 1]), type = 'l', lty = 2)

plot(density(DATA3$size[DATA3$surg == 0]), ylab = "", main = "",
     xlab = "", xlim = c(-1.5, 3), ylim = c(0, 1.6) )
points(density(DATA3$size[DATA3$surg == 1]), type = 'l', lty = 2)

plot(density(DATA4$size[DATA4$surg == 0]), ylab = "", main = "",
     xlab = "", xlim = c(-1.5, 3), ylim = c(0, 1.6) )
points(density(DATA4$size[DATA4$surg == 1]), type = 'l', lty = 2)

plot(density(DATA5$size[DATA5$surg == 0]), ylab = "", main = "",
     xlab = "", xlim = c(-1.5, 3), ylim = c(0, 1.6) )
points(density(DATA5$size[DATA5$surg == 1]), type = 'l', lty = 2)

#------ gender ------------------------------------------------------------
H1.sex <- c(sum(DATA1$sex[DATA1$surg == 0] == "Male"),
            sum(DATA1$sex[DATA1$surg == 0] == "Female"),
            sum(DATA1$sex[DATA1$surg == 1] == "Male"),
            sum(DATA1$sex[DATA1$surg == 1] == "Female") )
barplot(H1.sex, names.arg = c("M", "F", "M", "F"), ylab = "Gender", 
        col = c("darkgray", "darkgray", "light gray", "light gray"), 
        main = "", border = "black")
legend("topleft", c("Non-surgery", "Surgery"),
       col = c("darkgray", "light gray"),
       pch = 15, bty = "o", cex = 0.6)

H2.sex<-c(sum(DATA2$sex[DATA2$surg == 0] == "Male"),
          sum(DATA2$sex[DATA2$surg == 0] == "Female"),
          sum(DATA2$sex[DATA2$surg == 1] == "Male"),
          sum(DATA2$sex[DATA2$surg == 1] == "Female"))
barplot(H2.sex,names.arg=c("M","F","M","F"),ylab="", 
        col=c("darkgray","darkgray","light gray","light gray"),
        main="",border="black")

H3.sex<-c(sum(DATA3$sex[DATA3$surg==0]=="Male"),
          sum(DATA3$sex[DATA3$surg==0]=="Female"),
          sum(DATA3$sex[DATA3$surg==1]=="Male"),
          sum(DATA3$sex[DATA3$surg==1]=="Female"))
barplot(H3.sex, names.arg = c("M", "F", "M", "F"), ylab = "", 
        col = c("darkgray", "darkgray", "light gray", "light gray"),
        main = "", border = "black")


H4.sex<-c(sum(DATA4$sex[DATA4$surg == 0] == "Male"),
          sum(DATA4$sex[DATA4$surg == 0] == "Female"),
          sum(DATA4$sex[DATA4$surg == 1] == "Male"),
          sum(DATA4$sex[DATA4$surg == 1] == "Female") )
barplot(H4.sex, names.arg = c("M", "F", "M", "F"), ylab = "",
        col = c("darkgray", "darkgray", "light gray", "light gray"),
        main = "", border = "black")

H5.sex<-c(sum(DATA5$sex[DATA5$surg == 0] == "Male"),
          sum(DATA5$sex[DATA5$surg == 0] == "Female"),
          sum(DATA5$sex[DATA5$surg == 1] == "Male"),
          sum(DATA5$sex[DATA5$surg == 1] == "Female") )
barplot(H5.sex, names.arg = c("M", "F", "M", "F"), ylab = "", 
        col = c("darkgray", "darkgray", "light gray", "light gray"),
        main = "", border = "black")
dev.off()
##=========================================================================
## Table 2 in Supplementary 
## test between surgery group and non-surgery group after stratification
##=========================================================================
age_S1 <- ks.test(DATA1$age[DATA1$surg == 0], DATA1$age[DATA1$surg == 1], exact = TRUE)
age_S2 <- ks.test(DATA2$age[DATA2$surg == 0], DATA2$age[DATA2$surg == 1], exact = TRUE)
age_S3 <- ks.test(DATA3$age[DATA3$surg == 0], DATA3$age[DATA3$surg == 1], exact = TRUE)
age_S4 <- ks.test(DATA4$age[DATA4$surg == 0], DATA4$age[DATA4$surg == 1], exact = TRUE)
age_S5 <- ks.test(DATA5$age[DATA5$surg == 0], DATA5$age[DATA5$surg == 1], exact = TRUE)
age_KS_test <- c(age_S1$p.value, age_S2$p.value, age_S3$p.value, age_S4$p.value, age_S5$p.value)

size_S1 <- ks.test(DATA1$size[DATA1$surg == 0], DATA1$size[DATA1$surg == 1], exact = TRUE)
size_S2 <- ks.test(DATA2$size[DATA2$surg == 0], DATA2$size[DATA2$surg == 1], exact = TRUE)
size_S3 <- ks.test(DATA3$size[DATA3$surg == 0], DATA3$size[DATA3$surg == 1], exact = TRUE)
size_S4 <- ks.test(DATA4$size[DATA4$surg == 0], DATA4$size[DATA4$surg == 1], exact = TRUE)
size_S5 <- ks.test(DATA5$size[DATA5$surg == 0], DATA5$size[DATA5$surg == 1], exact = TRUE)
size_KS_test <- c(size_S1$p.value, size_S2$p.value, size_S3$p.value, size_S4$p.value, size_S5$p.value)

age_W1 <- wilcox.test(DATA1$age[DATA1$surg == 0], DATA1$age[DATA1$surg == 1])
age_W2 <- wilcox.test(DATA2$age[DATA1$surg == 0], DATA2$age[DATA1$surg == 1])
age_W3 <- wilcox.test(DATA3$age[DATA1$surg == 0], DATA3$age[DATA1$surg == 1])
age_W4 <- wilcox.test(DATA4$age[DATA1$surg == 0], DATA4$age[DATA1$surg == 1])
age_W5 <- wilcox.test(DATA5$age[DATA1$surg == 0], DATA5$age[DATA1$surg == 1])
age_WR_test <- c(age_W1$p.value, age_W2$p.value, age_W3$p.value, age_W4$p.value, age_W5$p.value)

size_W1 <- wilcox.test(DATA1$size[DATA1$surg == 0], DATA1$size[DATA1$surg == 1])
size_W2 <- wilcox.test(DATA2$size[DATA1$surg == 0], DATA2$size[DATA1$surg == 1])
size_W3 <- wilcox.test(DATA3$size[DATA1$surg == 0], DATA3$size[DATA1$surg == 1])
size_W4 <- wilcox.test(DATA4$size[DATA1$surg == 0], DATA4$size[DATA1$surg == 1])
size_W5 <- wilcox.test(DATA5$size[DATA1$surg == 0], DATA5$size[DATA1$surg == 1])
size_WR_test <- c(size_W1$p.value, size_W2$p.value, size_W3$p.value, size_W4$p.value, size_W5$p.value)

table2_supplemetary <- rbind(age_KS_test, size_KS_test, age_WR_test, size_WR_test)
colnames(table2_supplemetary) <- c("d1", "d2", "d3", "d4", "d5")

print(table2_supplemetary)


## ------ save results --------------------------------------------------------
save(table5_sum, file = "case_data/results/table5.RData")
save(table6_sum, file = "case_data/results/table6.RData")
save(table2_supplemetary, file = "case_data/results/table2_supplementary.RData")


