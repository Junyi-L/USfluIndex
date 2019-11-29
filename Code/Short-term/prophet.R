
# Prophet model after logit transform

## the design is based on the "ili_national" example in
## https://github.com/reichlab/article-disease-pred-with-kcde/blob/master/inst/code/prediction/sarima-prediction.R



library(prophet)
library(data.table)
library(lubridate)
library(logitnorm)
library(here)
set.seed(1)

load(file = here("./Data/data_holidays.RData"))
logit_FUN <- function(x){
  qlogis(x/100)
}

logistic_FUN <- function(x){
  plogis(x) * 100
}

data <- data.table(data)
holidaytab <- data[(x + y) > 0, ]
holidaytab[, holiday := ifelse(x == TRUE, "week22", 
                               "week23")]

prophet_change <- data.frame(holiday = holidaytab$holiday, ds = holidaytab$time)

data <- data.frame(data)
prediction_target_var <- "weighted_ili"
logit_prediction_target <- logit_FUN(data[, prediction_target_var])
time_target <- data[, "time"]


# prophetfit_control <- prophet(
#   yearly.seasonality = TRUE, weekly.seasonality = FALSE, daily.seasonality = FALSE,
#   holidays = prophet_change,
#   mcmc.samples = 0, # invokes rstan::optimizing (fast MAP estimation)
#   interval.width = 0.95, fit = FALSE)
# prophetfit <- fit.prophet(
#   m = prophetfit_control,
#   df = data.frame(ds = time_target, y = logit_prediction_target)
# )

##################################### 1:4 step ahead forecast
ptm <- proc.time()
all_prediction_horizons <- 1 : 4
results_row_ind <- 1L
prediction_time_inds <- which(data$train == FALSE)
num_rows <- length(all_prediction_horizons) *
  length(prediction_time_inds)
quantile_matrix <- matrix(nrow = length(prediction_time_inds),ncol = 99)

prophetfit_control <- prophet(
  yearly.seasonality = TRUE, weekly.seasonality = FALSE, daily.seasonality = FALSE,
  holidays = prophet_change,
  mcmc.samples = 0, # invokes rstan::optimizing (fast MAP estimation)
  interval.width = 0.95, fit = FALSE)

data_set_results <- data.frame(data_set = data[data$train == FALSE,],
                               model = "prophet",
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
                               stringsAsFactors = FALSE,row.names = NULL)


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
      data[prediction_time_ind, prediction_target_var]
    
    new_data <- data.frame(ds = time_target[seq_len(max(0, analysis_time_ind))],
                           y = logit_prediction_target[seq_len(max(0, analysis_time_ind))])
    
    prophetfit <- fit.prophet(
      m = prophetfit_control,
      df = new_data
    )
    
    future <- data.frame(ds = time_target[prediction_time_ind])
    predict_result <- predict(prophetfit,future)
    
    ## Get mean and variance of predictive distribution on log scale
    
    data_set_results$mu[results_row_ind] <- predictive_logit_mean <-
      as.numeric(predict_result$yhat) 
    data_set_results$sigma[results_row_ind] <- predictive_logit_sd <- 
      as.numeric(with(predict_result, ((yhat_upper-yhat_lower)/2)/qnorm(0.975)))
    
    
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
      abs(pt_predi - observed_prediction_target)
    
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

## save results
ProphetResults <- list(result = data_set_results,
                       quantile_matrix = quantile_matrix,
                       run_time = run_time, 
                       last_fit = prophetfit)

saveRDS(ProphetResults, file = here("./Results/Forecast_ph1-4/ProphetResults.rds"))
