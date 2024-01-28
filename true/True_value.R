###------Monte Carlo Method for caculating true RACE and SPCE---------------------
remove(list = ls())
set.seed(12345)

###------------------------------------------------------------------------
### Set alpha0=0 and alpha*=beta*=log(1.5) (Weak association)
###------------------------------------------------------------------------

u1 <- 5   # upper bound of censored time, Tau_C
tau <- u1 - 0.5
nn <- 500         # sample size
repli <- 500       # the number of replications

m1 <- 2  #used in the baseline hazard function under treat z=1
m0 <- 6  #used in the baseline hazard function under treat z=0
lam <- 0.5  #used in the baseline hazard function
n <- 1000000  #sample size

true_calculate <- function(alpha0_true, alpha11_true, beta11_true){
  
  alpha1_true <- log(1.5)
  alpha2_true <- -log(1.5)
  alpha3_true <- log(1.5)
  alpha4_true <- -log(1.5)
  alpha5_true <- log(2)
  alpha6_true <- log(2)
  alpha7_true <- log(1.5)
  alpha8_true <- log(1.5)
  alpha9_true <- log(1.25)
  alpha10_true <- log(1.25)
  alp <- c(alpha0_true, alpha1_true, alpha2_true, alpha3_true, alpha4_true, 
           alpha5_true, alpha6_true, alpha7_true, alpha8_true, alpha9_true,
           alpha10_true, alpha11_true) 
  
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
  bet <- c(beta1_true, beta2_true, beta3_true, beta4_true, 
           beta5_true, beta6_true, beta7_true, beta8_true, 
           beta9_true, beta10_true, beta11_true)
  
  
  #===============generate dataset=====================
  
  ran_prop <- runif(n, min = 0, max = 1)
  x1 <- rbinom(n, 1, 0.5)
  x2 <- rbinom(n, 1, 0.5)
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rbinom(n, 1, 0.5)
  x5 <- runif(n, min = -1, max = 1)
  x6 <- runif(n, min = -1, max = 1)
  x7 <- runif(n, min = -1, max = 1)
  x8 <- runif(n, min = -1, max = 1)
  x9 <- runif(n, min = -1, max = 1)
  x10 <- runif(n, min = -1, max = 1)
  
  # When all the subjects were exposed z=1 ----------------------------------------------------
  z <- 1
  zz_1 <- cbind(1,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,z)
  ucure_prob1 <- exp(zz_1%*%alp)/( 1+exp(zz_1%*%alp) )
  ucure_indc1 <- rbinom(n, 1, ucure_prob1)  #generate uncure indcator 
  
  fail_time1 <- rep(0,n)
  for (i in 1:n) {
    temp <- (-log(ran_prop[i])/(lam*exp(beta1_true*x1[i]+beta2_true*x2[i]+beta3_true*x3[i]
                                        +beta4_true*x4[i]+beta5_true*x5[i]+beta6_true*x6[i]
                                        +beta7_true*x7[i]+beta8_true*x8[i]+beta9_true*x9[i]
                                        +beta10_true*x10[i]+beta11_true*z)))^(1/m1)
    while( ( temp > tau )==TRUE ){
      a <- runif(1, min = 0, max = 1)
      temp <- (-log(a)/(lam*exp(beta1_true*x1[i]+beta2_true*x2[i]+beta3_true*x3[i]
                                +beta4_true*x4[i]+beta5_true*x5[i]+beta6_true*x6[i]
                                +beta7_true*x7[i]+beta8_true*x8[i]+beta9_true*x9[i]
                                +beta10_true*x10[i]+beta11_true*z)))^(1/m1)
    }
    fail_time1[i] <- temp
  }
  fail_time1[ucure_indc1==0] <- u1 #generate survival time under the treatment z=1
  race_true1 <- mean(fail_time1)    #true race under z=1 
  cure_rate1 <- mean(1-ucure_indc1) #true cure rate under z=1
  
  # spce at t=0.8 under z=1
  tt1 <- 0.8
  Su1 <- exp(-lam*tt1^m1*exp(beta1_true*x1+beta2_true*x2+beta3_true*x3
                             +beta4_true*x4+beta5_true*x5+beta6_true*x6+beta7_true*x7
                             +beta8_true*x8+beta9_true*x9+beta10_true*x10+beta11_true*z))
  S1_tt1 <- mean(Su1*ucure_indc1+1-ucure_indc1)
  
  # spce at t=2.4 under z=1
  tt2<-2.4
  Su2<-exp(-lam*tt2^m1*exp(beta1_true*x1+beta2_true*x2+beta3_true*x3
                           +beta4_true*x4+beta5_true*x5+beta6_true*x6
                           +beta7_true*x7+beta8_true*x8+beta9_true*x9
                           +beta10_true*x10+beta11_true*z))
  S1_tt2<- mean(Su2*ucure_indc1+1-ucure_indc1)
  
  
  ##----  When all the subjects were exposed z=1 ----------------------------------------------------------------------
  z<-0
  zz_0<-cbind(1,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,z)
  ucure_prob0<-exp( zz_0%*%alp)/( 1+exp(zz_0%*%alp) )
  ucure_indc0<-rbinom(n,1,ucure_prob0)  
  
  fail_time0 <-rep(0,n)
  for (i in 1:n)   {
    temp<-(-log(ran_prop[i])/(lam*exp(beta1_true*x1[i]+beta2_true*x2[i]+beta3_true*x3[i]
                                      +beta4_true*x4[i]+beta5_true*x5[i]+beta6_true*x6[i]
                                      +beta7_true*x7[i]+beta8_true*x8[i]+beta9_true*x9[i]
                                      +beta10_true*x10[i]+beta11_true*z)))^(1/m0)
    while( ( temp > tau )==TRUE ){
      a<- runif(1, min = 0, max = 1)
      temp<-(-log(a)/(lam*exp(beta1_true*x1[i]+beta2_true*x2[i]+beta3_true*x3[i]+
                                beta4_true*x4[i]+beta5_true*x5[i]+beta6_true*x6[i]
                              +beta7_true*x7[i]+beta8_true*x8[i]+beta9_true*x9[i]
                              +beta10_true*x10[i]+beta11_true*z)))^(1/m0)
    }
    fail_time0[i]<-temp
  }
  
  fail_time0[ucure_indc0==0]<-u1  # generate survival time T under z=0
  race_true0 <- mean( fail_time0 )  # true race under z=0
  cure_rate0 <- 1 - mean( ucure_indc0 ) # cure rate under z=0
  
  # spce at t=0.8 under z=0
  tt1 <- 0.8
  Su01 <- exp(-lam*tt1^m0*exp(beta1_true*x1+beta2_true*x2+beta3_true*x3+beta4_true*x4+beta5_true*x5+beta6_true*x6+beta7_true*x7+beta8_true*x8+beta9_true*x9+beta10_true*x10+beta11_true*z))
  S0_tt1 <- mean(Su01*ucure_indc0+1-ucure_indc0) 
  
  # spce at t=0.8 under z=0 
  tt2 <- 2.4
  Su02 <- exp(-lam*tt2^m0*exp(beta1_true*x1+beta2_true*x2+beta3_true*x3+beta4_true*x4+beta5_true*x5+beta6_true*x6+beta7_true*x7+beta8_true*x8+beta9_true*x9+beta10_true*x10+beta11_true*z))
  S0_tt2 <- mean(Su02*ucure_indc0+1-ucure_indc0)
  
  # computate ACE and SPCE -----------------------------------------------------------
  race_true <- race_true1 - race_true0      # race mean(T1)-mean(T0)
  spce_1 <- S1_tt1 - S0_tt1                   # spce at t=0.8
  spce_2 <- S1_tt2 - S0_tt2                   # spce at t=2.4
  spce_cure <- cure_rate1 - cure_rate0      # spce for cure rate 
  
  return(c(race_true, spce_1, spce_2, spce_cure))
  
}

summ_weak_alpha0 <- true_calculate(alpha0_true = 0, alpha11_true = log(1.5), beta11_true = log(1.5))
summ_weak_alpha05 <- true_calculate(alpha0_true  = 0.5, alpha11_true = log(1.5), beta11_true = log(1.5))
summ_weak_alpha1 <- true_calculate(alpha0_true  = 1, alpha11_true = log(1.5), beta11_true = log(1.5))

summ_weak_true <- round(rbind(summ_weak_alpha0, summ_weak_alpha05, summ_weak_alpha1),4)
colnames(summ_weak_true) <- c("race","spce0.8","spce2.4","cure")


summ_strong_alpha0 <- true_calculate(alpha0_true = 0, alpha11_true = log(2), beta11_true = log(2))
summ_strong_alpha05 <- true_calculate(alpha0_true  = 0.5, alpha11_true = log(2), beta11_true = log(2))
summ_strong_alpha1 <- true_calculate(alpha0_true  = 1, alpha11_true = log(2), beta11_true = log(2))

summ_strong_true <- round(rbind(summ_strong_alpha0, summ_strong_alpha05, summ_strong_alpha1),4)
colnames(summ_strong_true) <- c("race","spce0.8","spce2.4","cure")


#save.image("~/true/True_value.RData")


save(MyRsult,file = "C:/Users/zoewa/Desktop/simulation_ps/Code_Data_re/simulation_study/MSE_intermediate_results/weak/MSE_ACE_SPCEcure_weak_alp1_PE50_MPSBeran.RData")



