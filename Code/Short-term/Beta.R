library(betareg)
library(data.table)
library(lubridate)
library(here)

load(file = here("./Data/data_holidays.RData"))
source(file = here("./Code/Beta_forecast_p.R"))
lags <- 1

data <- data.table(data)
for (i in 1 : lags){
  data[, paste0("p", i) := logit_FUN(shift(weighted_ili,n = i, type = "lag"))]
}


Try_BetaReg5 <- function(S_mean, S_Precision, lags, data){
  mean_covar1 <- paste0(c(rep("sin_InPeriod",S_mean),rep("cos_InPeriod",S_mean)),1 : S_mean)
  mean_covar2 <- c("x", "y")
  mean_covar <- c(mean_covar1, mean_covar2)
  precision_covar <- paste0(c(rep("sin_InPeriod",S_Precision), rep("cos_InPeriod",S_Precision)), 1 : S_Precision)
  
  mean_AR <- paste("weighted_ili_org ~", paste(paste0("p", 1 : lags), collapse = " + "))
  mean_model <- paste(mean_AR, paste(mean_covar,collapse = " + "), sep = " + ")
  precision_model <- paste(precision_covar, collapse = " + ")
  full_model <- paste(mean_model, "|", precision_model)
  
  res <- try(betareg(full_model,
                     data = data))
  return(res)
}

train_data <- data[data$train == TRUE, ]
n <- sum(is.na(train_data$weighted_ili) != 1)
res <- list()
i = 1
combi <- character()
LL <- numeric()
AIC <- numeric()
AICc <- numeric()
npar <- numeric()
for(h in 1:5){
  for(k in 1:5){
    combi[i] <- paste0("S_v = ", h, "S_{\\phi} = ", k)
    model_hk_try <- Try_BetaReg5(S_mean = h, S_Precision = k, data = train_data, lags = lags)
    if(is(model_hk_try,"try-error")){
      npar[i] <- 0
      LL[i] <- 0
      AIC[i] <- 0
      AICc[i] <- 0
      res[[i]] <- list(S_mean = h, S_Precision = k, 0)
      
    }else{
      model_hk <- model_hk_try
      npar[i] <- h * 2 + k * 2 + 2 + 2 + lags
      LL[i] <- model_hk$loglik
      AIC[i] <- -2 * model_hk$loglik + 2 * npar[i]
      AICc[i] <- -2 * model_hk$loglik + 2 * npar[i] + 2 * npar[i] * ( npar[i] + 1)/(n -  npar[i] - 1)
      res[[i]] <- list(S_mean = h, S_Precision = k, model_hk)
      
    }
    i <- i + 1
    
  }
}

Vali_table <- data.frame(npar,LL,AIC,AICrank = rank(AIC), AICc, AICcrank = rank(AICc))
rownames(Vali_table) <- combi
Vali_table[Vali_table$AICcrank < 6,]


# npar       LL       AIC AICrank      AICc AICcrank
# S_v = 2S_{\\phi} = 4   17 3090.534 -6147.068       5 -6146.090        5
# S_v = 3S_{\\phi} = 3   17 3093.613 -6153.227       4 -6152.249        4
# S_v = 3S_{\\phi} = 4   19 3101.970 -6165.940       1 -6164.722        1
# S_v = 4S_{\\phi} = 3   19 3095.927 -6153.855       3 -6152.637        3
# S_v = 4S_{\\phi} = 4   21 3103.787 -6165.574       2 -6164.089        2
#-------------------------------------------------------- ph1-4 forecast--------------------------------------------------
data <- data.table(data)

ptm <- proc.time()

all_prediction_horizons <- 1:4
prediction_time_inds <- which(data$train == FALSE)
num_rows <- length(all_prediction_horizons) *
  length(prediction_time_inds)

results_row_ind <- 1L
count_error <- 0L
count_warning <- 0L
quantile_matrix <- matrix(nrow = length(prediction_time_inds), ncol = 99)

S_mean <- 3
S_Precision <- 4
lags <- 1

mean_covar1 <- paste0(c(rep("sin_InPeriod",S_mean),rep("cos_InPeriod",S_mean)),1 : S_mean)
mean_covar2 <- c("x", "y")
mean_covar <- c(mean_covar1, mean_covar2)

precisioin_covar <- paste0(c(rep("sin_InPeriod",S_Precision), rep("cos_InPeriod",S_Precision)), 1 : S_Precision)
mean_AR <- paste("weighted_ili_org ~", paste(paste0("p", 1 : lags), collapse = " + "))
mean_model <- paste(mean_AR, paste(mean_covar,collapse = " + "), sep = " + ")
precision_model <- paste(precisioin_covar, collapse = " + ")
full_model <- paste(mean_model, "|", precision_model)

max_S <- max(S_mean, S_Precision)
covar <- c(paste0(c(rep("sin_InPeriod",max_S), rep("cos_InPeriod",max_S)), 1 : max_S),mean_covar2)

data_set_results <- data.frame(data_set = data[data$train == FALSE,],
                               model = paste0("BetaReg", "(", S_mean, ",", S_Precision, ")"),
                               prediction_horizon = rep(NA_integer_, num_rows),
                               prediction_time = rep(NA, num_rows),
                               log_score = rep(NA_real_, num_rows),
                               pt_pred = rep(NA_real_, num_rows),
                               AE = rep(NA_real_, num_rows),
                               interval_pred_lb_95 = rep(NA_real_,  num_rows),
                               interval_pred_ub_95 = rep(NA_real_, num_rows),
                               interval_pred_lb_50 = rep(NA_real_, num_rows),
                               interval_pred_ub_50 = rep(NA_real_, num_rows),
                               precision = rep(NA_real_, num_rows),
                               shape1 = rep(NA_real_, num_rows),
                               shape2 = rep(NA_real_, num_rows),
                               DS_score = rep(NA_real_, num_rows),
                               stringsAsFactors = FALSE,
                               vari = rep(NA_real_, num_rows),
                               row.names = NULL)
class(data_set_results$prediction_time) <- class(data$time)
for(prediction_horizon in all_prediction_horizons) {
  for(prediction_time_ind in prediction_time_inds) {
    ## Set values describing case in data_set_results
    data_set_results$prediction_horizon[results_row_ind] <-
      prediction_horizon
    data_set_results$prediction_time[results_row_ind] <-
      data$time[prediction_time_ind]
    
    ## Get index of analysis time in data set
    ## (time from which we predict forward)
    analysis_time_ind <- prediction_time_ind - prediction_horizon
    
    ## Observed value at prediction time -- used in calculating log
    ## score and absolute error
    observed_prediction_target <-
      data[prediction_time_ind, weighted_ili_org]
    
    new_data <- data[
      seq_len(max(0, analysis_time_ind)), ]
    
    Update_BetaReg <- betareg(full_model,
                              data = new_data,
                              link = "logit")
    
    
    if (prediction_horizon == 1){
      predict_result <- 
        predict(Update_BetaReg, data[prediction_time_ind, ], type = "response")
      data_set_results$pt_pred[results_row_ind] <- predict_result * 100
      data_set_results$precision[results_row_ind] <- precision <-
        predict(Update_BetaReg, data[prediction_time_ind, ], type = "precision")
      
      data_set_results$shape1[results_row_ind] <- shape1 <- predict_result * precision
      data_set_results$shape2[results_row_ind] <- shape2 <- precision - shape1
      
      ## Compute log score of distribution prediction
      
      data_set_results$log_score[results_row_ind] <-
        dbeta(observed_prediction_target,
              shape1 = shape1,
              shape2 = shape2, 
              ncp = 0, 
              log = TRUE) - log(100)
      
      
      
      ## Compute absolute error of point prediction
      data_set_results$AE[results_row_ind] <-
        abs(data_set_results$pt_pred[results_row_ind] -
              observed_prediction_target * 100)
      
      ## Compute DSS of distribution prediction
      vari <- predict_result * (1 - predict_result)/(1 + precision)
      data_set_results$DS_score[results_row_ind] <- log(vari) + 
        (predict_result - observed_prediction_target)^2/vari + 
        4 * log(10)
      data_set_results$vari[results_row_ind] <- vari
      
      ## Compute prediction interval bounds
      data_set_results[results_row_ind,
                       c("interval_pred_lb_95",
                         "interval_pred_ub_95", 
                         "interval_pred_lb_50", 
                         "interval_pred_ub_50")] <-
        predict(Update_BetaReg, data[prediction_time_ind, ],
                type = "quantile", 
                at = c(0.05, 
                       0.95,
                       0.25, 
                       0.75)) * 100
      
      # save quantiles of predictive distribution for 1-week ahead prediction
      quantile_matrix[results_row_ind,] <-  qbeta(p = (1:99)/100,
                                                  shape1 = shape1,
                                                  shape2 = shape2, 
                                                  ncp = 0, 
                                                  log = FALSE) * 100
      
      
    }else{
      samples <- Beta_Forecast_p(object = Update_BetaReg,
                                 nsim = 10000, 
                                 seed = 1,
                                 ph = prediction_horizon,
                                 p = lags,
                                 start_value = data[(analysis_time_ind - lags + 1) : analysis_time_ind, weighted_ili],
                                 coefs_subset = data[analysis_time_ind + (1 : prediction_horizon),
                                                     ..covar])
      
      data_set_results[results_row_ind,
                       c("pt_pred",
                         "interval_pred_lb_95", 
                         "interval_pred_ub_95", 
                         "interval_pred_lb_50", 
                         "interval_pred_ub_50")] <-
        quantile(samples$pt[, prediction_horizon],
                 probs = c(0.5,
                           0.05,
                           0.95,
                           0.25, 
                           0.75))
      
      data_set_results$AE[results_row_ind] <-
        abs(data_set_results$pt_pred[results_row_ind] -
              observed_prediction_target * 100)
      
      data_set_results$vari[results_row_ind] <- var(samples$pt[, prediction_horizon])
      data_set_results$DS_score[results_row_ind] <- log(data_set_results$vari[results_row_ind]) + 
        (data_set_results$pt_pred[results_row_ind] - 
           observed_prediction_target * 100)^2/data_set_results$vari[results_row_ind]
      
      data_set_results$log_score[results_row_ind] <- cal_log_score(samples,
                                                                   (observed_prediction_target * 100),
                                                                   ph = prediction_horizon)
    }
    
    ## Increment results row
    results_row_ind <- results_row_ind + 1L
    print(results_row_ind)
  } # prediction_time_ind
}  # prediction_horizon
run_time <- proc.time() - ptm
BetaRegression <- list(result = data_set_results,
                       run_time = run_time,
                       quantile_matrix = quantile_matrix)
saveRDS(BetaRegression,
        file = here("./Results/Forecast_ph1-4/Beta.rds"))



