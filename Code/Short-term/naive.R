# Prophet model after logit transform

## the design is based on the "ili_national" example in
## https://github.com/reichlab/article-disease-pred-with-kcde/blob/master/inst/code/prediction/sarima-prediction.R


library(logitnorm)
library(here)

load(file = here("./Data/data_holidays.RData"))
logit_FUN <- function(x){
  qlogis(x/100)
}

logistic_FUN <- function(x){
  plogis(x) * 100
}


ptm <- proc.time()
prediction_target_var <- "weighted_ili"
all_prediction_horizons <- 1 : 4
results_row_ind <- 1L
prediction_time_inds <- which(data$train == FALSE)
num_rows <- length(all_prediction_horizons) *
  length(prediction_time_inds)
quantile_matrix <- matrix(nrow = length(prediction_time_inds),ncol = 99)

data_set_results <- data.frame(data_set = data[data$train == FALSE,],
                               model = "naive",
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
                               stringsAsFactors = FALSE,
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
      data[prediction_time_ind, prediction_target_var]
    
    prediction_week <- data[prediction_time_ind, ]$week
    prediction_year <- data[prediction_time_ind, ]$year
    # if new_week == 53, past week data is not enough. use data on week 52.
    if(prediction_week == 53) {
      week_idx <- which(data$week %in% c(52,53) & data$year < prediction_year) 
    }else week_idx <- which(data$week == prediction_week & data$year < prediction_year)
    # data_set_results[results_row_ind,c("mu","sigma")] <- 
    #   MASS::fitdistr(as.numeric(na.omit(data[week_idx,prediction_target_var])), 
    #                  densfun = "lognormal")$estimate
    
    
    data_set_results[results_row_ind,c("mu","sigma")] <- 
      MASS::fitdistr(logit_FUN(as.numeric(na.omit(data[week_idx,prediction_target_var]))), 
                     densfun = "normal")$estimate
    ## Compute log score of distribution prediction
    
    data_set_results$log_score[results_row_ind] <-
      dlogitnorm(observed_prediction_target/100,
                 mu = data_set_results[results_row_ind,]$mu,
                 sigma = data_set_results[results_row_ind,]$sigma,
                 log = TRUE) - log(100)
    
    
    ## Compute point prediction
    data_set_results$pt_pred[results_row_ind] <- pt_predi <-
      momentsLogitnorm(mu = data_set_results[results_row_ind,]$mu,
                       sigma = data_set_results[results_row_ind,]$sigma)["mean"] * 100
    
    
    ## Compute absolute error of point prediction
    data_set_results$AE[results_row_ind] <-
      abs(pt_predi - observed_prediction_target)
    
    ## Compute DS score
    data_set_results$var[results_row_ind] <- vari <-  momentsLogitnorm(mu = data_set_results[results_row_ind,]$mu,
                                                                       sigma = data_set_results[results_row_ind,]$sigma)["var"]  *10000
    
    data_set_results$DS_score[results_row_ind] <- log(vari) + 
      (observed_prediction_target - pt_predi)^2/vari
    
    ## Compute prediction interval bounds
    data_set_results[results_row_ind,
                     c("interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] <-
      qlogitnorm(c(0.025, 0.975, 0.25, 0.75),
                 mu = data_set_results[results_row_ind,]$mu,
                 sigma = data_set_results[results_row_ind,]$sigma) * 100
    
    if(prediction_horizon == 1){
      quantile_matrix[results_row_ind,] <-  qlogitnorm((1:99)/100,
                                                       mu = data_set_results[results_row_ind,]$mu,
                                                       sigma = data_set_results[results_row_ind,]$sigma) * 100
      
    }
    ## Increment results row
    results_row_ind <- results_row_ind + 1L
  } # prediction_time_ind
} # prediction_horizon

run_time <- proc.time() - ptm

NaiveResults <- list(result = data_set_results,
                     quantile_matrix = quantile_matrix,
                     run_time = run_time)
saveRDS(NaiveResults, file = here("./Results/Forecast_ph1-4/NaiveResults.rds"))
