# SARIMA model after logit transform
library(lubridate)
library(forecast)
library(logitnorm)
load(file = "./Data/data_holidays.RData")

logit_FUN <- function(x){
  qlogis(x/100)
}

logistic_FUN <- function(x){
  plogis(x) * 100
}

prediction_target_var <- "weighted_ili"

# remove missing values at the beginning.
missing <- is.na(data$weighted_ili)
firstnonmiss <- head(which(!missing), 1)
data <- data[firstnonmiss : dim(data)[1], ]

data_train <- data[data$train == TRUE, ]

logit_prediction_target_train <- logit_FUN(data_train[, prediction_target_var])
logit_prediction_target_train <- ts(logit_prediction_target_train, frequency = 52)

sarima_fit <-
  auto.arima(logit_prediction_target_train,trace = TRUE, allowdrift = TRUE)

# Save fit result for peak week prediction.
saveRDS(sarima_fit,file = "./Results/Peak/sarima_fit.rds")
model_order <- sarima_fit$arma[c(1, 6, 2, 3, 7, 4, 5)]

# check whether Arima return same fitting with auto.arima
check_fit <- Arima(
  logit_prediction_target_train,
  order = model_order[1:3],
  seasonal = list(order = model_order[4:6], period = model_order[7]),
  include.mean = "intercept" %in% names(sarima_fit$coef),
  include.drift = "drift" %in% names(sarima_fit$coef),
  lambda = sarima_fit$lambda)

#######################################################################
ptm <- proc.time()
all_prediction_horizons <- 1 : 4

logit_prediction_target <- logit_FUN(data[, prediction_target_var])
logit_prediction_target <- ts(logit_prediction_target,frequency = 52)

## Row indices in data corresponding to times at which we want to make a prediction.
prediction_time_inds <- which(data$train == FALSE)
num_rows <- length(all_prediction_horizons) *
  length(prediction_time_inds)

quantile_matrix <- matrix(nrow = length(prediction_time_inds),ncol = 99)

data_set_results <- data.frame(data_set = data[data$train == FALSE,],
                               model = "SARIMA",
                               prediction_horizon = rep(NA_integer_, num_rows),
                               prediction_time = rep(NA, num_rows),
                               pt_pred = rep(NA_real_, num_rows),
                               AE = rep(NA_real_, num_rows),
                               interval_pred_lb_95 = rep(NA_real_,  num_rows),
                               interval_pred_ub_95 = rep(NA_real_, num_rows),
                               interval_pred_lb_50 = rep(NA_real_, num_rows),
                               interval_pred_ub_50 = rep(NA_real_, num_rows),
                               mu = rep(NA_real_, num_rows),
                               sigma = rep(NA_real_, num_rows),
                               var = rep(NA_real_, num_rows),
                               log_score = rep(NA_real_, num_rows),
                               DS_score = rep(NA_real_, num_rows),
                               stringsAsFactors = FALSE,row.names = NULL)


class(data_set_results$prediction_time) <- class(data$time)

results_row_ind <- 1L
count_error <- 0L
updated_sarima_fit <- sarima_fit
count_warning <- 0L

for(prediction_horizon in all_prediction_horizons) {
  for(prediction_time_ind in prediction_time_inds) {
    data_set_results$prediction_horizon[results_row_ind] <-
      prediction_horizon
    data_set_results$prediction_time[results_row_ind] <-
      data$time[prediction_time_ind]
    
    ## Get index of analysis time in data set
    ## (time from which we predict forward).
    analysis_time_ind <- prediction_time_ind - prediction_horizon
    
    ## Observed value at prediction time -- used in calculating log
    ## score, absolute error and DSS.
    observed_prediction_target <-
      data[prediction_time_ind, prediction_target_var]
    
    # Data used for reestimation and forecast.
    new_data <- logit_prediction_target[
      seq_len(max(0, analysis_time_ind))]
    
    # Reestimate model at each step.
    # If error or warning occurs, use the estimation result at previous successful step and update new data.
    updated_sarima_fit_try <- try(Arima(
      new_data,
      order = model_order[1:3],
      seasonal = list(order = model_order[4:6], period = model_order[7]),
      include.mean = "intercept" %in% names(sarima_fit$coef),
      include.drift = "drift" %in% names(sarima_fit$coef),
      lambda = sarima_fit$lambda)
    )      
    
    
    if(is(updated_sarima_fit_try,"try-error")){
      updated_sarima_fit <- Arima(new_data, model = updated_sarima_fit)
      
      count_error <- count_error + 1L
    } else {
      try_result <- tryCatch(summary(updated_sarima_fit_try),
                             warning = function(warning_message) {
                               return("warning")
                             },
                             error = function(error_message) {
                               return("error")
                             }
      )
      if(try_result[1] == "warning"){
        updated_sarima_fit <- Arima(new_data, model = updated_sarima_fit)
        count_warning <- count_warning + 1L
      } else if (try_result[1] == "error"){
        updated_sarima_fit <- Arima(new_data, model = updated_sarima_fit)
        count_error <- count_error + 1L
      }else   updated_sarima_fit <- updated_sarima_fit_try
    }
    
    predict_result <- 
      predict(updated_sarima_fit, n.ahead = prediction_horizon)
    
    data_set_results$mu[results_row_ind] <- predictive_logit_mean <- 
      as.numeric(predict_result$pred[prediction_horizon]) 
    data_set_results$sigma[results_row_ind] <- predictive_logit_sd <- 
      as.numeric(predict_result$se[prediction_horizon])
    
    ## Compute log score of distribution prediction using logit normal distribution.
    ## the original value (percentage/100) follows logit normal
    data_set_results$log_score[results_row_ind] <-
      dlogitnorm(observed_prediction_target/100,
                 mu = predictive_logit_mean,
                 sigma = predictive_logit_sd,
                 log = TRUE) - log(100)
    
    ## Compute point prediction
    data_set_results$pt_pred[results_row_ind] <- pt_predi <-
      momentsLogitnorm(mu = predictive_logit_mean,
                       sigma = predictive_logit_sd)["mean"] * 100
    
    ## Compute absolute error of point prediction
    data_set_results$AE[results_row_ind] <-
      abs(pt_predi - observed_prediction_target)
    
    # Compute variance of point prediction in original scale.
    data_set_results$var[results_row_ind] <- vari <- as.numeric(momentsLogitnorm(mu = predictive_logit_mean,
                                                                                 sigma = predictive_logit_sd)["var"]) *10000
    
    data_set_results$DS_score[results_row_ind] <- log(vari) +
      (observed_prediction_target - pt_predi)^2/vari
    
    ## Compute prediction interval bounds
    data_set_results[results_row_ind,
                     c("interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] <-
      qlogitnorm(c(0.025, 0.975, 0.25, 0.75),
                 mu = predictive_logit_mean,
                 sigma = predictive_logit_sd) * 100
    
    # For 1-week ahead forecast, quantiles at (1:99)/100 in each step are saved in quantile_matrix,
    # to make fan plot.
    if(prediction_horizon == 1){
      quantile_matrix[results_row_ind,] <-  qlogitnorm((1:99)/100,
                                                       mu = predictive_logit_mean,
                                                       sigma = predictive_logit_sd) * 100
      
    }
    ## Increment results row
    results_row_ind <- results_row_ind + 1L
    print(results_row_ind)
  } # prediction_time_ind
} # prediction_horizon

run_time <- proc.time() - ptm

SARIMA_Result <- list(count_error = count_error,
                      count_warning = count_warning,
                      result = data_set_results, 
                      last_fit = updated_sarima_fit,
                      first_fit = sarima_fit,
                      run_time = run_time,
                      quantile_matrix= quantile_matrix)

saveRDS(SARIMA_Result,file = "./Results/Forecast_ph1-4/SARIMAResults.rds")

####################################################################
