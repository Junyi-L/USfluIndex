library(lubridate)
library(plyr)
library(dplyr)
library(reshape)
library(forecast)
library(mvtnorm)
library(doMC)
library(here)

set.seed(1)

n_sims <- 10000

prediction_save_path <- file.path(here("./Results/Peak/"))
logit_FUN <- function(x){
  qlogis(x/100)
}

logistic_FUN <- function(x){
  plogis(x) * 100
}

sample_predictive_trajectories_arima <- function (object, h = ifelse(object$arma[5] > 1, 2 * object$arma[5], 
                                                                     10), level = c(80, 95), fan = FALSE, xreg = NULL, lambda = object$lambda, 
                                                  npaths = 5000, ...) 
{
  sim <- matrix(NA, nrow = npaths, ncol = h)
  
  for (i in 1:npaths) {
    sim[i, ] <- simulate(object,
                         nsim = h, bootstrap = TRUE, future = TRUE, xreg = xreg)
  }
  
  return(sim)
}

load(file = here("./Data/data_holidays.RData"))

data <- data.table(data)
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

analysis_seasons <- c("2014/2015","2015/2016", "2016/2017", "2017/2018")
first_analysis_time_season_week <- 10 # == week 40 of year
last_analysis_time_season_week <- 41 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead

prediction_target_var <- "weighted_ili"
data <- as.data.frame(data)
logit_prediction_target <- logit_FUN(data[, prediction_target_var])


ili_incidence_bins <- data.frame(
  lower = seq(from = 0, to = 13, by = 0.5),
  upper = c(seq(from = 0.5, to = 13, by = 0.5), Inf))

arima_fit <- readRDS(file = here("./Results/Peak/arima.rds"))
model_order <- arima_fit$arma[c(1, 6, 2, 3, 7, 4, 5)]

results <- cbind(
  expand.grid(analysis_seasons,
              seq(from = first_analysis_time_season_week, to = last_analysis_time_season_week - 1),
              stringsAsFactors = FALSE),
  matrix(NA,
         nrow = length(analysis_seasons) * (last_analysis_time_season_week - first_analysis_time_season_week),
         ncol = 3 * n_sims + 2)
) %>%
  `colnames<-`(c("analysis_time_season",
                 "analysis_time_season_week",
                 "peak_week_log_score",
                 "peak_height_log_score",
                 paste0("peak_week_", seq_len(n_sims)),
                 paste0("peak_height_", seq_len(n_sims)),
                 paste0("unbinned_peak_height_", seq_len(n_sims))))

ptm <- proc.time()

count_error <- 0
count_warning <- 0
updated_sarima_fit <- arima_fit
for(analysis_time_season in analysis_seasons) {
  ## get observed quantities, for computing log score
  observed_peak_height <- max(data[data$season == analysis_time_season, prediction_target_var])
  observed_peak_week_ind <- which((data$season == analysis_time_season) &
                                    (data[, prediction_target_var] == observed_peak_height))
  observed_peak_week <- data$season_week[observed_peak_week_ind]
  
  
  observed_peak_height <- which(
    ili_incidence_bins$lower <= observed_peak_height &
      ili_incidence_bins$upper > observed_peak_height)
  
  
  for(analysis_time_season_week in seq(from = first_analysis_time_season_week, to = last_analysis_time_season_week - 1)) {
    analysis_time_ind <- which(data$season == analysis_time_season &
                                 data$season_week == analysis_time_season_week)
    
    new_data <- logit_prediction_target[
      seq_len(max(0, analysis_time_ind))]
    
    new_cov <- sarima_cov_total[seq_len(max(0, analysis_time_ind)),]
    updated_sarima_fit_try <- try(Arima(
      new_data,
      order = model_order[1:3],
      seasonal = list(order = model_order[4:6], period = model_order[7]),
      xreg = new_cov)
    )      
    if(is(updated_sarima_fit_try, "try-error")){
      updated_sarima_fit <- Arima(new_data, model = updated_sarima_fit, xreg = new_cov)
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
        updated_sarima_fit <- Arima(new_data, model = updated_sarima_fit, xreg = new_cov)
        count_warning <- count_warning + 1L
      } else if (try_result[1] == "error"){
        updated_sarima_fit <- Arima(new_data, model = updated_sarima_fit, xreg = new_cov)
        count_error <- count_error + 1L
      }else   updated_sarima_fit <- updated_sarima_fit_try
    }
    
    ## simulate n_sims trajectories recursively from sarima
    max_prediction_horizon <-
      last_analysis_time_season_week + 1 -
      analysis_time_season_week
    sim_cov <- sarima_cov_total[(analysis_time_ind + 1) : (analysis_time_ind + max_prediction_horizon),]
    
    trajectory_samples <- sample_predictive_trajectories_arima(
      updated_sarima_fit,
      h = max_prediction_horizon,
      npaths = n_sims,
      xreg = sim_cov)
    
    trajectory_samples <- logistic_FUN(trajectory_samples)
    
    ## Augment trajectory samples with previous observed incidence values
    season_start_ind <- which(data$season == analysis_time_season &
                                data$season_week == 1)
    if(season_start_ind < analysis_time_ind) {
      trajectory_samples <- cbind(
        matrix(
          rep(data[seq(from = season_start_ind, to = analysis_time_ind), prediction_target_var], each = n_sims),
          nrow = n_sims
        ),
        trajectory_samples
      )
    }
    
    ## Get peak week and height at peak week for each simulated trajectory 
    results_save_row <- which(results$analysis_time_season == analysis_time_season &
                                results$analysis_time_season_week == analysis_time_season_week)
    
    peak_week_by_sim_ind <- apply(trajectory_samples, 1, which.max)
    results[results_save_row, paste0("peak_week_", seq_len(n_sims))] <-
      peak_week_by_sim_ind
    
    peak_week_height_by_sim_ind <- trajectory_samples[cbind(seq_len(n_sims), peak_week_by_sim_ind)]
    results[results_save_row, paste0("unbinned_peak_height_", seq_len(n_sims))] <-
      peak_week_height_by_sim_ind
    
    
    peak_week_height_by_sim_ind <- sapply(peak_week_height_by_sim_ind,
                                          function(height) {
                                            which(ili_incidence_bins$lower <= height &
                                                    ili_incidence_bins$upper > height)
                                          })
    
    results[results_save_row, paste0("peak_height_", seq_len(n_sims))] <-
      peak_week_height_by_sim_ind
    
    ## Get log scores
    results[results_save_row, "peak_week_log_score"] <- log(sum(peak_week_by_sim_ind == observed_peak_week)) - log(n_sims)
    results[results_save_row, "peak_height_log_score"] <- log(sum(peak_week_height_by_sim_ind == observed_peak_height)) - log(n_sims)
  }
}


run_time <- proc.time() - ptm
saveRDS(results,
        file.path(prediction_save_path,
                  paste0("peak-week-arima", ".rds")))

others <- list(run_time = run_time, count_error = count_error, count_warning = count_warning)
saveRDS(run_time,
        file.path(prediction_save_path,
                  paste0("peak-week-arima-others",  ".rds")))




