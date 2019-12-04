# KCDE short-term forecast
# Author: Evan Ray, modified by Junyi Lu.
# Date: 16 Aug 2017
# Code comes from 
## https://github.com/reichlab/article-disease-pred-with-kcde/blob/master/inst/code/prediction/kcde-prediction.R

# Modifications:
# - The loading data and saving file path are changed. 
# - Forecasts up to four-weeks ahead are made.
# - In order to calculate Dawid-Sabastiani score for each forecast, 
# 100000 samples are made to calculate the empirical variance.
# - For one-week ahead forecast, predictive quantiles are saved in quantile_matrix.

library(kcde)
library(pdtmvn)
library(plyr)
library(dplyr)
library(lubridate)
library(doMC)
library(here)


all_data_sets <- "ili_national"
all_prediction_horizons <- 1 : 4 # make predictions at every horizon from 1 to 4 weeks ahead

all_max_lags <- 1L # use incidence at times t^* and t^* - 1 to predict incidence after t^*
all_max_seasonal_lags <- 0L # not used
all_filtering_values <- FALSE # not used
all_differencing_values <- FALSE # not used
all_seasonality_values <- TRUE # specifications without and with periodic kernel
all_bw_parameterizations <-  "full" # specifications with diagonal and full bandwidth

num_cores <- 3L
registerDoMC(cores = num_cores)
results_path <- file.path(
  here( "./Results/Forecast_ph1-4/KCDEresults"))

load(file = here("./Data/data_holidays.RData"))


## Row indices in data corresponding to times at which we want to make a prediction
prediction_time_inds <- which(data$train == FALSE)


quantile_matrix <- matrix(nrow = length(prediction_time_inds),ncol = 99)
for(prediction_horizon in all_prediction_horizons){
  ptm <- proc.time()
  results_row_ind <- 1L
  
  ## Allocate data frame to store results for this prediction horizon
  num_rows <- length(all_max_lags) *
    length(all_max_seasonal_lags) *
    length(all_filtering_values) *
    length(all_differencing_values) *
    length(all_seasonality_values) *
    length(all_bw_parameterizations) *
    length(prediction_time_inds)
  
  ph_results <- data.frame(data_set = data[data$train == FALSE,],
                           prediction_horizon = rep(NA_integer_, num_rows),
                           max_lag = rep(NA_integer_, num_rows),
                           max_seasonal_lag = rep(NA_integer_, num_rows),
                           filtering = rep(NA, num_rows),
                           differencing = rep(NA, num_rows),
                           seasonality = rep(NA, num_rows),
                           bw_parameterization = rep(NA_character_, num_rows),
                           model = "kcde",
                           prediction_time = rep(NA, num_rows),
                           log_score = rep(NA_real_, num_rows),
                           pt_pred = rep(NA_real_, num_rows),
                           AE = rep(NA_real_, num_rows),
                           interval_pred_lb_95 = rep(NA_real_, num_rows),
                           interval_pred_ub_95 = rep(NA_real_, num_rows),
                           interval_pred_lb_50 = rep(NA_real_, num_rows),
                           interval_pred_ub_50 = rep(NA_real_, num_rows),
                           stringsAsFactors = FALSE,
                           vari = rep(NA_integer_, num_rows),
                           DS_score = rep(NA_integer_, num_rows)
  )
  class(ph_results$prediction_time) <- class(data$time)
  
  for(max_lag in all_max_lags) {
    for(max_seasonal_lag in all_max_seasonal_lags) {
      for(filtering in all_filtering_values) {
        for(differencing in all_differencing_values) {
          ## Set prediction target var
          
          if(differencing) {
            prediction_target_var <- "weighted_ili_ratio"
            orig_prediction_target_var <- "weighted_ili"
          } else {
            prediction_target_var <- "weighted_ili"
          } 
          
          
          for(seasonality in all_seasonality_values) {
            for(bw_parameterization in all_bw_parameterizations) {
              for(prediction_time_ind in prediction_time_inds) {
                ## Set values describing case in ph_results
                ph_results$prediction_horizon[results_row_ind] <-
                  prediction_horizon
                ph_results$max_lag[results_row_ind] <-
                  max_lag
                ph_results$max_seasonal_lag[results_row_ind] <-
                  max_seasonal_lag
                ph_results$filtering[results_row_ind] <-
                  filtering
                ph_results$differencing[results_row_ind] <-
                  differencing
                ph_results$seasonality[results_row_ind] <-
                  seasonality
                ph_results$bw_parameterization[results_row_ind] <-
                  bw_parameterization
                ph_results$prediction_time[results_row_ind] <-
                  data$time[prediction_time_ind]
                
                ## Load kcde_fit object.  Estimation was performed previously.
                case_descriptor <- paste0(
                  data_set,
                  "-prediction_horizon_", prediction_horizon,
                  "-max_lag_", max_lag,
                  "-max_seasonal_lag_", max_seasonal_lag,
                  "-filtering_", filtering,
                  "-differencing_", differencing,
                  "-seasonality_", seasonality,
                  "-bw_parameterization_", bw_parameterization
                )
                
                kcde_fit_file_path <- file.path(results_path,
                                                paste0("kcde_fit-", case_descriptor, ".rds"))
                kcde_fit <- readRDS(kcde_fit_file_path)
                
                
                ## fix rkernel_fn for pdtmvn-based kernel functions
                for(kernel_component_ind in seq(from = (as.logical(seasonality) + 1), to = length(kcde_fit$kcde_control$kernel_components))) {
                  kcde_fit$kcde_control$kernel_components[[kernel_component_ind]]$rkernel_fn <-
                    function(n,
                             conditioning_obs,
                             center,
                             bw,
                             bw_continuous,
                             conditional_bw_discrete,
                             conditional_center_discrete_offset_multiplier,
                             continuous_vars,
                             discrete_vars,
                             continuous_var_col_inds,
                             discrete_var_col_inds,
                             discrete_var_range_fns,
                             lower,
                             upper,
                             x_names,
                             ...) {
                      if(missing(conditioning_obs) || is.null(conditioning_obs)) {
                        log_conditioning_obs <- NULL
                      } else {
                        log_conditioning_obs <- log(conditioning_obs)
                      }
                      if(missing(bw_continuous)) {
                        bw_continuous <- NULL
                      }
                      if(missing(conditional_bw_discrete)) {
                        conditional_bw_discrete <- NULL
                      }
                      if(missing(conditional_center_discrete_offset_multiplier)) {
                        conditional_center_discrete_offset_multiplier <- NULL
                      }
                      
                      ## center parameter of pdtmvn_kernel is mean of log
                      ## mode of resulting log-normal distribution is
                      ## mode = exp(mu - bw %*% 1) (where 1 is a column vector of 1s)
                      ## therefore mu = log(mode) + bw %*% 1
                      reduced_x_names <- names(center)
                      inds_x_vars_in_orig_vars <- which(x_names %in% reduced_x_names)
                      x_names_for_call <- x_names[inds_x_vars_in_orig_vars]
                      
                      mean_offset <- apply(bw, 1, sum)[x_names %in% colnames(center)]
                      
                      return(exp(rpdtmvn_kernel(n = n,
                                                conditioning_obs = log_conditioning_obs,
                                                center = sweep(log(center)[, x_names_for_call, drop = FALSE], 2, mean_offset, `+`),
                                                bw = bw,
                                                bw_continuous = bw_continuous,
                                                conditional_bw_discrete = conditional_bw_discrete,
                                                conditional_center_discrete_offset_multiplier = conditional_center_discrete_offset_multiplier,
                                                continuous_vars = continuous_vars,
                                                discrete_vars = discrete_vars,
                                                continuous_var_col_inds = continuous_var_col_inds,
                                                discrete_var_col_inds = discrete_var_col_inds,
                                                discrete_var_range_fns = discrete_var_range_fns,
                                                lower = lower,
                                                upper = upper,
                                                x_names = x_names)[, reduced_x_names, drop = FALSE]))
                    }
                }
                
                
                ## Get index of analysis time in data set
                ## (time from which we predict forward)
                analysis_time_ind <- prediction_time_ind - prediction_horizon
                
                ## Compute log score
                observed_prediction_target <-
                  data[prediction_time_ind, prediction_target_var, drop = FALSE]
                colnames(observed_prediction_target) <-
                  paste0(prediction_target_var, "_horizon", prediction_horizon)
                ph_results$log_score[results_row_ind] <-
                  kcde_predict(
                    kcde_fit = kcde_fit,
                    prediction_data =
                      data[seq_len(analysis_time_ind), , drop = FALSE],
                    leading_rows_to_drop = 0L,
                    trailing_rows_to_drop = 0L,
                    additional_training_rows_to_drop = NULL,
                    prediction_type = "distribution",
                    prediction_test_lead_obs = observed_prediction_target,
                    log = TRUE
                  )
                
                if(differencing) {
                  ph_results$log_score[results_row_ind] <-
                    ph_results$log_score[results_row_ind] -
                    (abs(data[analysis_time_ind - 52, orig_prediction_target_var]))
                }
                
                ## Compute point prediction and interval predictions -- quantiles
                # ph_results[results_row_ind,
                #            c("pt_pred", "interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] <-
                #   kcde_predict(
                #     p = c(0.5, 0.025, 0.975, 0.25, 0.75),
                #     n = 100000,
                #     kcde_fit = kcde_fit,
                #     prediction_data =
                #       data[seq_len(analysis_time_ind), , drop = FALSE],
                #     leading_rows_to_drop = 0L,
                #     trailing_rows_to_drop = 0L,
                #     additional_training_rows_to_drop = NULL,
                #     prediction_type = "quantile"
                #   )
                
                # Simulate samples to calculate the empirical variance for Dawid-Sebastiani score.
                samples <-
                  kcde_predict(
                    n = 100000,
                    kcde_fit = kcde_fit,
                    prediction_data =
                      data[seq_len(analysis_time_ind), , drop = FALSE],
                    leading_rows_to_drop = 0L,
                    trailing_rows_to_drop = 0L,
                    additional_training_rows_to_drop = NULL,
                    prediction_type = "sample"
                  )
                
                if(prediction_horizon == 1){
                  quantile_matrix[results_row_ind,] <-  quantile(samples,probs = (1:99)/100)
                  
                }
                ph_results[results_row_ind,
                           c("pt_pred", "interval_pred_lb_95", "interval_pred_ub_95",
                             "interval_pred_lb_50", "interval_pred_ub_50")] <- 
                  quantile(samples,probs = c(0.5, 0.025, 0.975, 0.25, 0.75))
                
                ph_results[results_row_ind,"vari"] <- var(samples)
                
                # Compute Dawid-Sebastiani score
                ph_results[results_row_ind,"DS_score"] <- log(ph_results[results_row_ind,"vari"]) + 
                  ((observed_prediction_target - ph_results[results_row_ind,"pt_pred"])^2)/ph_results[results_row_ind,"vari"]
                
                if(differencing) {
                  ph_results[results_row_ind,
                             c("pt_pred", "interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] <-
                    ph_results[results_row_ind,
                               c("pt_pred", "interval_pred_lb_95", "interval_pred_ub_95", "interval_pred_lb_50", "interval_pred_ub_50")] * 
                    data[prediction_time_ind - 52, prediction_target_var]
                }
                
                ## Compute absolute error of point prediction
                ph_results$AE[results_row_ind] <-
                  as.numeric(abs(ph_results$pt_pred[results_row_ind] -
                                   observed_prediction_target))
                
                ## Increment results row
                results_row_ind <- results_row_ind + 1L
                print(results_row_ind)
              } # prediction_time_ind
            } # bw_parameterization
          } # seasonality
        } # differencing
      } # filtering
    } # max_seasonal_lag
  } # max_lag
  run_time <- proc.time() - ptm
  
  saveRDS(ph_results, file = file.path(
    here( "./Results/Forecast_ph1-4/KCDEresults"),
    paste0("kcde-predictions-ph_", prediction_horizon, ".rds")))
  saveRDS(run_time, file = file.path(
    here("./Results/Forecast_ph1-4/KCDEresults"),
    paste0("kcde-predictions-time-ph_", prediction_horizon, ".rds")))
  
} # prediction_horizon 

saveRDS(quantile_matrix,file = here("./Results/Forecast_ph1-4/KCDEresults/quantile_matrix.rds"))

