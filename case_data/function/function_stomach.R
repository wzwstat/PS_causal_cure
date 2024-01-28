##===============================================================================
## Functions of the stomach studies
##===============================================================================

##-------------------------------------------------------------------
## Stratification function by propensity score
##-------------------------------------------------------------------
datastratify_stomach <- function(data){ 
  ps_model <- glm(surg ~ age + sex + size, data = data, 
                  family = binomial(link="logit"), na.action = na.pass)    #estimate PS
  data$pscore <- predict(ps_model, newdata = data, type = "response")  
  orderps <- sort(data$pscore)    #sort data based on PS
  q <- quantile(orderps, probs = seq(0, 1, 1/J))  
  q[1] <- 0
  
  Data <- data_order <- NULL
  
  for (i in 1:J) {
    Data[[i]] <- data[(data$pscore > q[i]) & (data$pscore <= q[i + 1]) , ]
    data_order[[i]] <- Data[[i]][order(Data[[i]]$time, decreasing = F), ]
  }
  return(data_order = data_order)
}


##--------------------------------------------------------------------
## Trapezoidal Rule for Integration
##--------------------------------------------------------------------
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
## Kernel function for a binary covariate 
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
  }
  else
    if (side == "right") {
      right_ord <- which(Xord > x_fix)
      ind <- right_ord[kneighbour]
      neighbour_k <- Xord[ind]
      distk <- abs(neighbour_k - x_fix)
    }
  return(dkneighbour = distk)
}

##----------------------------------------------------------
## This function computes the pilot bandwidth g.
## It is used in the computation of the bootstrap sample standard error
##----------------------------------------------------------
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
## Function to obtain results about Table 7 
## in MSE_computations_strong_table1&2.R 
## and in table7.R 
##=================================================================================
Get_table= function(){
  
  p_value_ace <- 2*pt(abs(summ_est$ACE_est/boot_se$Boot_ace_se), df=Inf, lower.tail=FALSE )
  p_value_spce <- 2*pt(abs(summ_est$SPCE_est/boot_se$Boot_spcecure_se), df=Inf, lower.tail=FALSE )
  
  down_ace <- summ_est$ACE_est - qnorm(0.975,0,1) * boot_se$Boot_ace_se
  up_ace <- summ_est$ACE_est + qnorm(0.975,0,1) * boot_se$Boot_ace_se
  
  down_spce <- summ_est$SPCE_cure_est - qnorm(0.975,0,1) * boot_se$Boot_spcecure_se
  up_spce <- summ_est$SPCE_cure_est + qnorm(0.975,0,1) * boot_se$Boot_spcecure_se
  
  
  estimator <- as.numeric(summ_est)
  se_est <- as.numeric(boot_se)
  lower <- c(down_ace, down_spce)
  upper <- c(up_ace, up_spce)
  p_value <- c(p_value_ace, p_value_spce)
  
  table_result <- cbind(estimator, se_est, lower, upper, p_value)
  row_name <- c("ACE", "SPCE")
  col_name <- c("Est", "SE", "Lower", "Upper", "p_value")
  dimnames(table_result) <- list(row_name, col_name)
  #' timeend <- Sys.time()
  #' runningtime <- timeend - timestart
  #' cat('----------------------------------------------', '\n',
  #'     'Run Time is            :', runningtime, 'secs or minus', '\n',
  #'     #'Right-censoring rata is:', data$RC, '\n', 
  #'     #'Cure rate is           :', data$CR, '\n',
  #'     '---------------------------------------------', '\n'
  #' )
  #' a <- print(runningtime)
  #' 
  #' print(round(table_result, 6))
  return(list(table_result = round(table_result, 6)))
}
