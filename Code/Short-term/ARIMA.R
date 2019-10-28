library(lubridate)
library(forecast)
library(data.table)
library(logitnorm)
library(here)

load(file = here("./Data/data_holidays.RData"))
data <- data.table(data)
logit_FUN <- function(x){
  qlogis(x/100)
}

logistic_FUN <- function(x){
  plogis(x) * 100
}
#-----------------------------------------------------------------------------
# remove missing values at the beginning.
missing <- is.na(data$weighted_ili)
firstnonmiss <- head(which(!missing), 1)
data <- data[firstnonmiss : dim(data)[1], ]

sarima_cov_total <- data[, .(sin_InPeriod1, 
                             cos_InPeriod1, 
                             sin_InPeriod2, 
                             cos_InPeriod2,
                             sin_InPeriod3, 
                             cos_InPeriod3,
                             sin_InPeriod4, 
                             cos_InPeriod4,
                             x,
                             y)]
sarima_cov_total <- as.matrix(sarima_cov_total)

data_train <- data[data$train == TRUE,]

logit_prediction_target_train <- logit_FUN(data_train[, weighted_ili])
sarima_cov <- sarima_cov_total[data$train == TRUE, , drop = FALSE]

arima_fit <-
  auto.arima(logit_prediction_target_train, xreg = as.matrix(sarima_cov), seasonal = FALSE, trace = TRUE)
# Save fit result for peak week prediction.
saveRDS(arima_fit, file = here("./Results/Peak/arima.rds"))

model_order <- arima_fit$arma[c(1, 6, 2, 3, 7, 4, 5)]

# check_fit <- Arima(
#   logit_prediction_target_train,
#   order = model_order[1:3],
#   seasonal = list(order = model_order[4:6], period = model_order[7]),
#   include.mean = "intercept" %in% names(arima_fit$coef),
#   include.drift = "drift" %in% names(arima_fit$coef),
#   lambda = arima_fit$lambda,
#   xreg = as.matrix(sarima_cov))
# prediction --------------------------------------------------------

ptm <- proc.time()
all_prediction_horizons <- 1:4
## Row indices in data corresponding to times at which we want to make a prediction
prediction_time_inds <- which(data$train == FALSE)

logit_prediction_target <- logit_FUN(data[, weighted_ili])

num_rows <- length(all_prediction_horizons) *
  length(prediction_time_inds)

quantile_matrix <- matrix(nrow = length(prediction_time_inds), ncol = 99)
data_set_results <- data.frame(data_set =  data[data$train == FALSE,],
                               model = "ARIMA",
                               prediction_horizon = rep(NA_integer_, num_rows),
                               prediction_time = rep(NA, num_rows),
                               log_score = rep(NA_real_, num_rows),
                               pt_pred = rep(NA_real_, num_rows),
                               AE = rep(NA_real_, num_rows),
                               interval_pred_lb_95 = rep(NA_real_,  num_rows),
                               interval_pred_ub_95 = rep(NA_real_, num_rows),
                               interval_pred_lb_50 = rep(NA_real_, num_rows),
                               interval_pred_ub_50 = rep(NA_real_, num_rows),
                               mu = rep(NA_real_, num_rows),
                               sigma = rep(NA_real_, num_rows),
                               var = rep(NA_real_, num_rows),
                               DS_score = rep(NA_real_, num_rows),
                               stringsAsFactors = FALSE, row.names = NULL)


class(data_set_results$prediction_time) <- class(data$time)

results_row_ind <- 1L
count_error <- 0L
updated_sarima_fit <- arima_fit

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
      data[prediction_time_ind, weighted_ili]
    
    new_data <- logit_prediction_target[
      seq_len(max(0, analysis_time_ind))]
    
    new_sarima_cov <- sarima_cov_total[
      seq_len(max(0, analysis_time_ind)), ]
    
    updated_sarima_fit_try <- try(Arima(
      new_data,
      order = model_order[1:3],
      seasonal = list(order = model_order[4:6], period = model_order[7]),
      include.mean = "intercept" %in% names(arima_fit$coef),
      include.drift = "drift" %in% names(arima_fit$coef),
      lambda = arima_fit$lambda,
      xreg = new_sarima_cov)
    )
    
    
    
    if(is(updated_sarima_fit_try,"try-error")){
      updated_sarima_fit_try <- Arima(new_data, model = updated_sarima_fit, xreg = new_sarima_cov)
      count_error <- count_error + 1L
    } else updated_sarima_fit <- updated_sarima_fit_try
    
    predict_result <- predict(updated_sarima_fit, n.ahead = prediction_horizon, 
                              newxreg = matrix(sarima_cov_total[prediction_time_ind, ], nrow = 1))
    
    data_set_results$mu[results_row_ind] <- predictive_logit_mean <- 
      as.numeric(predict_result$pred[prediction_horizon]) 
    data_set_results$sigma[results_row_ind] <- predictive_logit_sd <- 
      as.numeric(predict_result$se[prediction_horizon])
    
    ## Compute log score of distribution prediction
    data_set_results$log_score[results_row_ind] <-
      dlogitnorm(observed_prediction_target/100,
                 mu = predictive_logit_mean,
                 sigma = predictive_logit_sd,
                 log = TRUE) - log(100)
    
    
    ## Compute point prediction
    data_set_results$pt_pred[results_row_ind] <- pt_predi <-
      momentsLogitnorm(mu = predictive_logit_mean,
                       sigma = predictive_logit_sd)["mean"] * 100
    
    ## Compute DS score
    data_set_results$var[results_row_ind] <- vari <- as.numeric(momentsLogitnorm(mu = predictive_logit_mean,
                                                                                 sigma = predictive_logit_sd)["var"]) *10000
    
    data_set_results$DS_score[results_row_ind] <- log(vari) +
      (observed_prediction_target - pt_predi)^2/vari
    
    ## Compute absolute error of point prediction
    data_set_results$AE[results_row_ind] <-
      as.numeric(abs(pt_predi -
                       observed_prediction_target))
    
    ## Compute prediction interval bounds
    
    data_set_results[results_row_ind,
                     c("interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] <-
      qlogitnorm(c(0.025, 0.975, 0.25, 0.75),
                 mu = predictive_logit_mean,
                 sigma = predictive_logit_sd) * 100
    
    if(prediction_horizon == 1){
      quantile_matrix[results_row_ind,] <-  qlogitnorm((1:99)/100,
                                                       mu = predictive_logit_mean,
                                                       sigma = predictive_logit_sd) * 100
    }
    
    ## Increment results row
    results_row_ind <- results_row_ind + 1L
  } # prediction_time_ind
} # prediction_horizon

run_time <- proc.time() - ptm

ARIMA_sin_cos <- list(count_error = count_error,
                      result = data_set_results, 
                      last_fit = updated_sarima_fit, 
                      first_fit = arima_fit,
                      run_time = run_time,
                      quantile_matrix = quantile_matrix)

saveRDS(ARIMA_sin_cos, file = here("./Results/Forecast_ph1-4/ArimaResults.rds"))

