
library(lubridate)
library(plyr)
library(dplyr)
library(reshape)
library(mvtnorm)
library(doMC)
set.seed(1)

logit_FUN <- function(x){
  qlogis(x/100)
}

logistic_FUN <- function(x){
  plogis(x) * 100
}

load(file = "./data_holidays.RData")
prediction_target_var <- "weighted_ili"
prediction_target <- data[, prediction_target_var]
week_target <- data[, "week"]
n_sims <- 10000
prediction_save_path <- file.path("./Peak/")

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

week_fun <- function(base_data, new_week, nsim){
  # if new_week == 53, past week data is not enough. use data on week 52.
  if(new_week != 53){
    estimate <- MASS::fitdistr(logit_FUN(as.numeric(na.omit(base_data[base_data[,1] == new_week, 2]))), 
                               densfun = "normal")$estimate
  } else{
    estimate <- MASS::fitdistr(logit_FUN(as.numeric(na.omit(base_data[base_data[,1] %in% c(52, 53), 2]))), 
                               densfun = "normal")$estimate
  }
  sims <- rlogitnorm(nsim, mu = estimate[1], sigma = estimate[2]) * 100
  return(sims)
}

ptm <- proc.time()

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
    
    new_data <- data.frame(week = week_target[seq_len(max(0, analysis_time_ind))],
                           weighted_ili = prediction_target[seq_len(max(0, analysis_time_ind))])
    
    
    max_prediction_horizon <-
      last_analysis_time_season_week + 1 -
      analysis_time_season_week
    
    weeks_tab <- week_target[(analysis_time_ind + 1) : (analysis_time_ind + max_prediction_horizon)]
    trajectory_samples <- mapply(week_fun, new_week = weeks_tab, nsim = n_sims, MoreArgs = list(base_data = new_data))
    
    
    
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
run_time <- proc.time() - ptm

saveRDS(results,
        file.path(prediction_save_path,
                  paste0("peak-week-naive",".rds")))


others <- list(run_time = run_time,count_error = count_error, count_warning = count_warning)
saveRDS(run_time,
        file.path(prediction_save_path,
                  paste0("peak-week-naive","-others.rds")))


