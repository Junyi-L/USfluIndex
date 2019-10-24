library(data.table)
library(xtable)
library(lubridate)
library(forecast)

load(file = "./Data/data_holidays.RData")
logit_FUN <- function(x){
  qlogis(x/100)
}

logistic_FUN <- function(x){
  plogis(x) * 100
}
#-----------------------------------------------------------------------------
train_data <- data[data$train == TRUE,]
train_data <- data.table(train_data)
sarima_cov_train <- train_data[, .(sin_InPeriod1, 
                                   cos_InPeriod1, 
                                   sin_InPeriod2, 
                                   cos_InPeriod2,
                                   sin_InPeriod3, 
                                   cos_InPeriod3,
                                   sin_InPeriod4, 
                                   cos_InPeriod4,
                                   sin_InPeriod5, 
                                   cos_InPeriod5
)]
sarima_cov_train <- as.matrix(sarima_cov_train)
sarima_cov_train2 <- train_data[, .( x, y)]
sarima_cov_train2 <- as.matrix(sarima_cov_train2)

logit_prediction_target <- logit_FUN(train_data[, weighted_ili])

missing <- is.na(logit_prediction_target)
firstnonmiss <- head(which(!missing), 1)

sarima_cov_train2 <- sarima_cov_train2[firstnonmiss : dim(sarima_cov_train2)[1], ]
sarima_cov_train <- sarima_cov_train[firstnonmiss : dim(sarima_cov_train)[1], ]
logit_prediction_target <- logit_prediction_target[firstnonmiss:length(logit_prediction_target)]

Try_ARIMA <- function(S, target, regressor, regressor2){
  idx <- c("x", "y")
  sarima_cov <- regressor[, 1:(2 * S)]
  sarima_cov2 <- regressor2[, idx]
  
  sarima_cov <- data.frame(sarima_cov, sarima_cov2)
  sarima_cov <- as.matrix(sarima_cov)
  res <- try(auto.arima(target, xreg = sarima_cov, seasonal = FALSE)
  )
  return(res)
}


i = 1
n <- sum(is.na(train_data$weighted_ili) != 1)
combi <- character()
LL <- numeric()
AIC <- numeric()
AICc <- numeric()
npar <- numeric()
term <- character()
for(h in 1:5){
  combi[i] <- paste0("$S = ", h, "$")
  model_h <- Try_ARIMA(S = h,
                       target =  logit_prediction_target,
                       regressor = sarima_cov_train,
                       regressor2 = sarima_cov_train2)
  if(is(model_h,"try-error")){
    npar[i] <- 0
    LL[i] <- 0
    AIC[i] <- 0
    AICc[i] <- 0
    term[i] <- 0
  }else{
    npar[i] <- sum(model_h$arma[c(1, 2)]) + 2 * h + 2 + 1
    LL[i] <- model_h$loglik
    AIC[i] <- model_h$aic
    AICc[i] <- model_h$aic + 2 * npar[i] * ( npar[i] + 1)/(n -  npar[i] - 1)
    term[i] <- paste(model_h$arma[c(1, 6, 2)], collapse = ",")
  }
  i = i + 1
  
}
Vali_table <- data.frame(term, 
                         npar,
                         LL,
                         AIC,
                         AICrank = rank(AIC),
                         AICc,
                         AICcrank = rank(AICc))

rownames(Vali_table) <- combi

Vali_table[Vali_table$AICcrank < 6,]

# $S = 1$ 5,1,0   10 337.3298 -654.6596       4 -654.3120        4
# $S = 2$ 5,1,0   12 341.5728 -659.1456       3 -658.6511        3
# $S = 3$ 5,1,0   14 348.1070 -668.2141       2 -667.5463        2
# $S = 4$ 5,1,0   16 350.8359 -669.6718       1 -668.8042        1
# $S = 5$     0    0   0.0000    0.0000       5    0.0000        5
