# Beta model with one lag, peak prediction

## the design is based on the "ili_national" example in
## https://github.com/reichlab/article-disease-pred-with-kcde/blob/master/inst/code/prediction/sarima-peak-prediction.R



library(betareg)
library(data.table)
library(lubridate)
library(plyr)
library(dplyr)
library(reshape)
library(mvtnorm)
library(doMC)
library(here)
source(file = here("./Code/Beta_forecast_p.R"))
load(file = here("./Data/data_holidays.RData"))
# The number of harmonics is chosen in Beta.R
data <- data.table(data)
n_sims <- 10000
prediction_save_path <- file.path(here("./Results/Peak/"))

lags <- 1
for (i in 1 : lags){
  data[, paste0("p", i) := logit_FUN(shift(weighted_ili,n = i, type = "lag"))]
}


analysis_seasons <- c("2014/2015","2015/2016", "2016/2017", "2017/2018")
first_analysis_time_season_week <- 10 # == week 40 of year
last_analysis_time_season_week <- 41 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead



ili_incidence_bins <- data.frame(
  lower = seq(from = 0, to = 13, by = 0.5),
  upper = c(seq(from = 0.5, to = 13, by = 0.5), Inf))

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


count_error <- 0
count_warning <- 0

ptm <- proc.time()

S_mean <- 3
S_Precision <- 4
lags <- 1

mean_covar <- paste0(c(rep("sin_InPeriod",S_mean), rep("cos_InPeriod",S_mean)), 1 : S_mean)
precisioin_covar <- paste0(c(rep("sin_InPeriod",S_Precision), rep("cos_InPeriod",S_Precision)), 1 : S_Precision)
mean_AR <- paste("weighted_ili_org ~", paste(paste0("p", 1 : lags), collapse = " + "))
mean_model <- paste(mean_AR, paste(mean_covar,collapse = " + "), sep = " + ")
precision_model <- paste(precisioin_covar, collapse = " + ")
full_model <- paste(mean_model, "|", precision_model)

max_S <- max(S_mean, S_Precision)
covar <- paste0(c(rep("sin_InPeriod",max_S), rep("cos_InPeriod",max_S)), 1 : max_S)

for(analysis_time_season in analysis_seasons) {
  ## get observed quantities, for computing log score
  observed_peak_height <- max(data[data$season == analysis_time_season,]$weighted_ili)
  observed_peak_week_ind <- which((data$season == analysis_time_season) &
                                    (data$weighted_ili == observed_peak_height))
  observed_peak_week <- data$season_week[observed_peak_week_ind]
  
  
  observed_peak_height <- which(
    ili_incidence_bins$lower <= observed_peak_height &
      ili_incidence_bins$upper > observed_peak_height)
  
  
  for(analysis_time_season_week in seq(from = first_analysis_time_season_week, to = last_analysis_time_season_week - 1)) {
    analysis_time_ind <- which(data$season == analysis_time_season &
                                 data$season_week == analysis_time_season_week)
    
    ## Update fit object with data up
    ## through analysis_time_ind
    new_data <- data[
      seq_len(max(0, analysis_time_ind)),]
    
    
    updated_Beta_Reg <- betareg(full_model,
                                data = new_data,
                                link = "logit")
    
    
    ## simulate n_sims trajectories 
    max_prediction_horizon <-
      last_analysis_time_season_week + 1 -
      analysis_time_season_week
    
    
    cov_subset <- data[analysis_time_ind + (1:max_prediction_horizon),
                       ..covar]
    
    trajectory_samples <- Beta_Forecast_p(object = updated_Beta_Reg,
                                          nsim = n_sims, 
                                          seed = 1,
                                          ph = max_prediction_horizon,
                                          p = lags,
                                          start_value = data[(analysis_time_ind - lags + 1) : analysis_time_ind, weighted_ili],
                                          coefs_subset = cov_subset)$pt
    
    
    ## Augment trajectory samples with previous observed incidence values
    season_start_ind <- which(data$season == analysis_time_season &
                                data$season_week == 1)
    if(season_start_ind < analysis_time_ind) {
      trajectory_samples <- cbind(
        matrix(
          rep(data[seq(from = season_start_ind, to = analysis_time_ind), ]$weighted_ili, each = n_sims),
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
saveRDS(results,
        file.path(prediction_save_path,
                  paste0("peak-week-beta", S_mean, "-", S_Precision,".rds")))

run_time <- proc.time() - ptm
others <- list(run_time = run_time,count_error = count_error, count_warning = count_warning)
saveRDS(run_time,
        file.path(prediction_save_path,
                  paste0("peak-week-beta",S_mean, "-", S_Precision,"-others.rds")))




