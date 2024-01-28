##===============================================================================
## Functions of the simulated scenarios
##===============================================================================

##--------------------------------------------------------------------------------
## Generate dataset
##--------------------------------------------------------------------------------
## the propensity score is generated from logistic function
getdata <- function(){ 
  
  ran_prop <- runif(nn, min = 0, max = 1)
  x1 <- rbinom(nn, 1, 0.5)
  x2 <- rbinom(nn, 1, 0.5)
  x3 <- rbinom(nn, 1, 0.5)
  x4 <- rbinom(nn, 1, 0.5)
  x5 <- runif(nn, min = -1, max = 1)
  x6 <- runif(nn, min = -1, max = 1)
  x7 <- runif(nn, min = -1, max = 1)
  x8 <- runif(nn, min = -1, max = 1)
  x9 <- runif(nn, min = -1, max = 1)
  x10 <- runif(nn, min = -1, max = 1)
  xx <- cbind(1,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
  
  treat_prob <- exp(xx %*% gam)/( 1 + exp(xx %*% gam) )
  z <- rbinom(nn, 1, treat_prob)  #generate treatment  
  zz <- cbind(xx,z)
  
  ucure_prob <- exp( zz %*%  alp)/( 1 + exp(zz %*% alp) ) #uncure rate
  ucure_indc <- rbinom(nn,1,ucure_prob)  
  
  T1 <- rep(0,nn)
  for (i in 1:nn)   {
    temp <- (-log(ran_prop[i])/(lamb*exp(beta1_true*x1[i]+beta2_true*x2[i]+beta3_true*x3[i]+beta4_true*x4[i]
                                      +beta5_true*x5[i]+beta6_true*x6[i]+beta7_true*x7[i]+beta8_true*x8[i]
                                      +beta9_true*x9[i]+beta10_true*x10[i]+beta11_true*z[i])))^(1/mu1)
    while( ( temp > tau )==TRUE ){
      a <- runif(1, min = 0, max = 1)
      temp <- (-log(a)/(lamb*exp(beta1_true*x1[i]+beta2_true*x2[i]+beta3_true*x3[i]+beta4_true*x4[i]
                              +beta5_true*x5[i]+beta6_true*x6[i]+beta7_true*x7[i]+beta8_true*x8[i]
                              +beta9_true*x9[i]+beta10_true*x10[i]+beta11_true*z[i])))^(1/mu1)
    }
    T1[i] <- temp
  }
  
  T0 <- rep(0,nn)
  for (i in 1:nn)   {
    temp <- (-log(ran_prop[i])/(lamb*exp(beta1_true*x1[i]+beta2_true*x2[i]+beta3_true*x3[i]+beta4_true*x4[i]
                                      +beta5_true*x5[i]+beta6_true*x6[i]+beta7_true*x7[i]+beta8_true*x8[i]
                                      +beta9_true*x9[i]+beta10_true*x10[i]+beta11_true*z[i])))^(1/mu0)
    while( ( temp > tau )==TRUE ){
      a <- runif(1, min = 0, max = 1)
      temp <- (-log(a)/(lamb*exp(beta1_true*x1[i]+beta2_true*x2[i]+beta3_true*x3[i]+beta4_true*x4[i]
                              +beta5_true*x5[i]+beta6_true*x6[i]+beta7_true*x7[i]+beta8_true*x8[i]
                              +beta9_true*x9[i]+beta10_true*x10[i]+beta11_true*z[i])))^(1/mu0)
    }
    T0[i] <- temp
  }
  
  fail_time <- obse_time <- NULL
  fail_time[z==1] <- obse_time[z==1] <- T1[z==1]  
  fail_time[z==0] <- obse_time[z==0] <- T0[z==0]  
  
  cens_time <- runif(nn, min = 0, max = u1) # generate censored time from U(0,5)
  
  delta <- (fail_time <= cens_time) & (ucure_indc==1)  #right censored 
  obse_time[!delta] <- cens_time[!delta]  #generate observed time
  
  data <- data.frame(
    obse_time=obse_time, delta=delta, x1=x1, x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,x8=x8, x9=x9,x10=x10,z=z)
  return(data)
}

## the propensity score is generated from probit function
getdata_MPS <- function(){ 
  
  ran_prop <- runif(nn, min = 0, max = 1)
  x1 <- rbinom(nn, 1, 0.5)
  x2 <- rbinom(nn, 1, 0.5)
  x3 <- rbinom(nn, 1, 0.5)
  x4 <- rbinom(nn, 1, 0.5)
  x5 <- runif(nn, min = -1, max = 1)
  x6 <- runif(nn, min = -1, max = 1)
  x7 <- runif(nn, min = -1, max = 1)
  x8 <- runif(nn, min = -1, max = 1)
  x9 <- runif(nn, min = -1, max = 1)
  x10 <- runif(nn, min = -1, max = 1)
  xx <- cbind(1,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
  
  esp_t<-rnorm(nn)
  
  treat_prob <-xx %*% gam +esp_t
  z <- ifelse(treat_prob<=0, 0, 1)  #generate treatment  
  zz <- cbind(xx,z)
  
  ucure_prob <- exp( zz %*%  alp)/( 1 + exp(zz %*% alp) ) #uncure rate
  ucure_indc <- rbinom(nn, 1, ucure_prob)  
  
  T1 <- rep(0,nn)
  for (i in 1:nn)   {
    temp <- (-log(ran_prop[i])/(lamb*exp(beta1_true*x1[i]+beta2_true*x2[i]+beta3_true*x3[i]+beta4_true*x4[i]
                                         +beta5_true*x5[i]+beta6_true*x6[i]+beta7_true*x7[i]+beta8_true*x8[i]
                                         +beta9_true*x9[i]+beta10_true*x10[i]+beta11_true*z[i])))^(1/mu1)
    while( ( temp > tau )==TRUE ){
      a <- runif(1, min = 0, max = 1)
      temp <- (-log(a)/(lamb*exp(beta1_true*x1[i]+beta2_true*x2[i]+beta3_true*x3[i]+beta4_true*x4[i]
                                 +beta5_true*x5[i]+beta6_true*x6[i]+beta7_true*x7[i]+beta8_true*x8[i]
                                 +beta9_true*x9[i]+beta10_true*x10[i]+beta11_true*z[i])))^(1/mu1)
    }
    T1[i] <- temp
  }
  
  T0 <- rep(0,nn)
  for (i in 1:nn)   {
    temp <- (-log(ran_prop[i])/(lamb*exp(beta1_true*x1[i]+beta2_true*x2[i]+beta3_true*x3[i]+beta4_true*x4[i]
                                         +beta5_true*x5[i]+beta6_true*x6[i]+beta7_true*x7[i]+beta8_true*x8[i]
                                         +beta9_true*x9[i]+beta10_true*x10[i]+beta11_true*z[i])))^(1/mu0)
    while( ( temp > tau )==TRUE ){
      a <- runif(1, min = 0, max = 1)
      temp <- (-log(a)/(lamb*exp(beta1_true*x1[i]+beta2_true*x2[i]+beta3_true*x3[i]+beta4_true*x4[i]
                                 +beta5_true*x5[i]+beta6_true*x6[i]+beta7_true*x7[i]+beta8_true*x8[i]
                                 +beta9_true*x9[i]+beta10_true*x10[i]+beta11_true*z[i])))^(1/mu0)
    }
    T0[i] <- temp
  }
  
  fail_time <- obse_time <- NULL
  fail_time[z == 1] <- obse_time[z == 1] <- T1[z == 1]
  fail_time[z == 0] <- obse_time[z == 0] <- T0[z == 0]
  
  cens_time <- runif(nn, min = 0, max = u1) # generate censored time from U(0,5)
  
  delta <- (fail_time <= cens_time) & (ucure_indc == 1)  #right censored 
  obse_time[!delta] <- cens_time[!delta]  #generate observed time
  
  data <- data.frame(
    obse_time=obse_time, delta=delta, x1=x1, x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,x8=x8, x9=x9,x10=x10,z=z)
  return(data)
}

##--------------------------------------------------------------------------------
## Stratification function by propensity score
##--------------------------------------------------------------------------------
datastratify <- function(data){ 
  
  ps_model <- glm(z ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, data = data,
                  family = binomial(link = "logit"),
                  na.action = na.pass  )    #estimate propensity score (PS)
  data$pscore <- predict(ps_model, newdata = data, type = "response")  
  
  orderps <- sort(data$pscore)    #sort data based on PS
  q <- quantile(orderps, probs = seq(0, 1, 1/J))  
  q[1] <- 0
  
  Data <- data_order <- NULL
  
  for(i in 1:J){
    Data[[i]] <- data[(data$pscore > q[i]) & (data$pscore <= q[i + 1]),]
    data_order[[i]] <- Data[[i]][order(Data[[i]]$obse_time, decreasing = F),]
  }
  
  return( data_order)
}

## For Beran estimator, no stratification, only J=1 with estimated PS 
datastratify_full <- function(data){ 
  
  ps_model <- glm(z ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10, data = data,
                  family = binomial(link = "logit"),
                  na.action = na.pass  )    #estimate propensity score (PS)
  data$pscore <- predict(ps_model, newdata = data, type = "response")  
  
  orderps <- sort(data$pscore)    #sort data based on PS
  q <- quantile(orderps, probs = seq(0, 1, 1/J))  
  q[1] <- 0
  
  Data <- data_order <- NULL
  
  for(i in 1:J){
    Data[[i]] <- data[(data$pscore > q[i]) & (data$pscore <= q[i + 1]),]
    data_order[[i]] <- Data[[i]][order(Data[[i]]$obse_time, decreasing = F),]
  }
  
  return( data_order)
}

## the propensity score is estimated by probit function
datastratify_probit <- function(data){ 
  
  ps_model <- glm(z ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, data = data,
                  family = binomial(link = "probit"),
                  na.action = na.pass  )    #estimate propensity score (PS)
  data$pscore <- predict(ps_model, newdata = data, type = "response")  
  
  orderps <- sort(data$pscore)    #sort data based on PS
  q <- quantile(orderps, probs = seq(0, 1, 1/J))  
  q[1] <- 0
  
  Data <- data_order<-NULL
  
  for(i in 1:J){
    Data[[i]] <- data[(data$pscore > q[i]) & (data$pscore <= q[i + 1]),]
    data_order[[i]] <- Data[[i]][order(Data[[i]]$obse_time, decreasing = F),]
  }
  
  return(data_order)
}

##--------------------------------------------------------------------------------
## Trapezoidal Rule for Integration
##--------------------------------------------------------------------------------
trapezoid <- function(x, y)
  sum(diff(x) * (y[-1] + y[-length(y)])) / 2


##--------------------------------------------------------------------------------
## Kernel function to compute the Nadaraya-Watson weights for continuous covariates
##--------------------------------------------------------------------------------
Kernal_c = function(x, hn) {
  f = 3 / 4 * (1 - (x / hn) ^ 2) * (abs(x / hn) <= 1)  # Epanechnikov kernel
  return(f)
}

##--------------------------------------------------------------------------------
## Kernel function to compute the Nadaraya-Watson weights for a binary covariate
##--------------------------------------------------------------------------------
Kernal_b = function(x) {
  f = ifelse(x == 0, 1, 0)
  return(f)
}

##---------------------------------------------------------------------------------
## This function computes the k-nearest neighbour. 
## It is used in the computation of the pilot bandwidth g in the gpilot.fun() below
##---------------------------------------------------------------------------------
dkneighbour <- function(x_fix = x_fix, Xord = Xord, side = side, kneighbour = kneighbour) {
  if (side == "left") {
    left_ord <- which(Xord < x_fix)
    ind <- left_ord[length(left_ord) - kneighbour + 1]
    neighbour_k <- Xord[ind]
    distk <- abs(neighbour_k - x_fix)
  } else
    if (side == "right") {
      right_ord <- which(Xord > x_fix)
      ind <- right_ord[kneighbour]
      neighbour_k <- Xord[ind]
      distk <- abs(neighbour_k - x_fix)
    }
  return(dkneighbour = distk)
}

##---------------------------------------------------------------------------------
## This function computes the pilot bandwidth g.
## It is used in the computation of the bootstrap sample standard error
## in the PS_Beran_boot() function of PSBeran_boots.R 
##---------------------------------------------------------------------------------
gpilotfun <- function(Xdata = Xdat, Xord = Xord, x_fix = x_fix, side = side, kneighbour = kneighbour) {
  n <- length(Xord)
  num_neigh_left <- sum(Xdata < x_fix)
  num_neigh_right <- sum(Xdata > x_fix)
  if ((num_neigh_left >= kneighbour) & (num_neigh_right >= kneighbour)) {
    neighbourk_right <- dkneighbour(x = x_fix, Xord = Xord, side = "right", kneighbour = kneighbour)
    neighbourk_left <- dkneighbour(x = x_fix, Xord = Xord, side = "left", kneighbour = kneighbour)
    g_pilot <- (neighbourk_right + neighbourk_left)/2 * (100^(1/9)) * (n^(-1/9))
  } else
    if (num_neigh_left >= kneighbour) {
      neighbourk_left <- dkneighbour(x = x_fix, Xord = Xord, side = "left", kneighbour = kneighbour)
      g_pilot <- (neighbourk_left * 2) * (100^(1/9)) * (n^(-1/9))
    } else
      if (num_neigh_right >= kneighbour) {
        neighbourk_right <- dkneighbour(x = x_fix, Xord = Xord, side = "right", kneighbour = kneighbour)
        g_pilot <- (neighbourk_right * 2) * (100^(1/9)) * (n^(-1/9))
      } else
        if ((num_neigh_left < kneighbour) & (num_neigh_right < kneighbour)) {
          g_pilot <- (max_x - min_x) * (100^(1/9)) * (n^(-1/9))
        }
  return(g_pilot)
}

##=================================================================================
## Function to obtain MSE about Table 1 and Table 2 
## in MSE_computations_table1&2.R 
##=================================================================================
Get_MSE <- function(){
  summ_est <- matrix(0, nrow = repli, ncol =2)
  
  for(i in 1:repli){
    summ_est[i,] <- as.numeric(MyRsult[[i]])
  }
  summ_est <- summ_est[is.finite(rowSums(summ_est)), ]
  summ_est_ave <- apply(na.omit(summ_est), 2, mean)
  summ_sd <- apply(na.omit(summ_est), 2, sd )    
  
  MSE <- NULL
  summ_True <- as.numeric(summ_true)
  
  MSE[1] <- mean((summ_est[, 1] - summ_True[1]) ^ 2)
  MSE[2] <- mean((summ_est[, 2] - summ_True[2]) ^ 2)
  pb <- (summ_est_ave - summ_True) / summ_True
  
  Summary_summ <- cbind(summ_True, summ_est_ave, summ_est_ave - summ_True, pb, summ_sd, MSE)
  row_name <- c("ACE", "SPCE")
  col_name <- c("Ture value", "Est", "Bias", "PB", "SD", "MSE")
  dimnames(Summary_summ) <- list(row_name, col_name)
  
  #print(round(Summary_summ, 4))
  return(list( Summary_summ = round(Summary_summ, 4) ))
}

##=================================================================================
## Function to obtain Table 3 and Table 4 in table3.R and table4.R 
## The Getsummary() computes the stimates, Bise, PB, SD, SE, CP and MSE.
##=================================================================================
Getsummary <- function(){
  summ_est <- matrix(0, nrow = repli, ncol =p)
  boot_se <- matrix(0, nrow = repli, ncol =p)
  
  for(i in 1:repli){
    summ_est[i, ] <- as.numeric(MyRsult[[i]])
    boot_se[i, ] <- as.numeric(MyRsult[[repli + i]])
  }
  
  summ_est_ave <- apply(summ_est, 2, mean)
  summ_sd <- apply(summ_est, 2, sd)
  summ_se_ave <- apply(boot_se, 2, mean)
  
  MSE <- summ_CP <- summ_CP2 <- NULL
  summ_True <- as.numeric(summ_true)
  
  for (ip in 1:p) {
    MSE[ip] <- mean((summ_est[, ip] - summ_True[ip]) ^ 2)
    summ_CP[ip] <- sum(((summ_est[, ip] - qnorm(0.975, 0.1) * boot_se[, ip]) <= summ_true[ip])
          & (summ_true[ip] <= (summ_est[, ip] + qnorm(0.975, 0.1) * boot_se[, ip]))) / repli
    summ_CP2[ip] <- sum(((summ_est[, ip] - 1.96 * boot_se[, ip]) <= summ_true[ip])
                       & (summ_true[ip] <= (summ_est[, ip] + 1.96 * boot_se[, ip]))) / repli
    
      }
  
  pb <- (summ_est_ave - summ_True) / summ_True
  
  Summary_result <- cbind(summ_True, summ_est_ave, summ_est_ave - summ_True,
                          pb, summ_sd, summ_se_ave, summ_CP, summ_CP2, MSE)
  row_name <- c("ACE", "SPCE_08", "SPCE_24", "SPCE_cure")
  col_name <- c("Ture value", "Est", "Bias", "PB", "SD", "SE", "CP", "CP2", "MSE")
  dimnames(Summary_result) <- list(row_name, col_name)
  
  #' timeend <- Sys.time()
  #' runningtime <- timeend - timestart
  #' cat('----------------------------------------------', '\n',
  #'     'Run Time is            :', runningtime, 'secs or minus', '\n',
  #'     #'Right-censoring rata is:', data$RC, '\n',
  #'     #'Cure rate is           :', data$CR, '\n',
  #'     '---------------------------------------------', '\n'
  #' )
  #' a <- print(runningtime)
  #' print(round(Summary_result, 4))
  
  return(list( Summary_result = round(Summary_result, 4) ))
}


