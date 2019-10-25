library(here)
#------------------------------------------------------------
peak_arima <- readRDS(file = here("./Results/Peak/peak-week-arima.rds"))
peak_sarima <- readRDS(file = here("./Results/Peak/peak-week-sarima.rds"))
peak_Beta <- readRDS(file = here("./Results/Peak/peak-week-beta3-4.rds"))
peak_Beta_lags <- readRDS(file = here("./Results/Peak/peak-week-beta-lags.rds"))

# 
peak_prophet <- readRDS(file = here("./Results/Peak/peak-week-prophet.rds"))
peak_naive <- readRDS(file = here("./Results/Peak/peak-week-naive.rds"))
peak_KCDE <- readRDS(file = here("./Results/Peak/peak-week-KCDE.rds"))


result_table <- list(peak_Beta,
                     peak_Beta_lags,
                     peak_KCDE,
                     peak_arima,
                     peak_sarima,
                     peak_prophet,
                     peak_naive)


model_name <- c("Beta(1)",
                "Beta(4)",
                "KCDE",
                "ARIMA",
                "SARIMA",
                "Prophet",
                "Naive")


result_table0 <- list()
for (i in 1 : length(result_table)){
  result_table0[[i]] <- result_table[[i]][,c("analysis_time_season",
                                             "analysis_time_season_week",
                                             "peak_week_log_score",
                                             "peak_height_log_score")]
}

peak_week_log_score <- numeric()
peak_height_log_score <- numeric()
peak_week_log_score_BP <- numeric()
peak_height_log_score_BP <- numeric()

peak_week_log_score_max <- numeric()
peak_height_log_score_max <- numeric()
peak_week_log_score_BP_max <- numeric()
peak_height_log_score_BP_max <- numeric()

########################################################################################
# indicator for weeks before peak 
analysis_seasons <- c("2014/2015","2015/2016", "2016/2017", "2017/2018")
observed_peak_week <- numeric(length(analysis_seasons))
observed_peak_height <- numeric(length(analysis_seasons))

load(file = here("./Data/data_holidays.RData"))


first_analysis_time_season_week <- 10 # == week 40 of year
last_analysis_time_season_week <- 41 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead

prediction_target_var <- "weighted_ili"

## get observed quantities, for computing log score
for(i in 1 : length(analysis_seasons)) {
  analysis_time_season <- analysis_seasons[i]
  observed_peak_height[i] <- max(data[data$season == analysis_time_season, prediction_target_var])
  observed_peak_week_ind <- which((data$season == analysis_time_season) &
                                    (data[, prediction_target_var] == observed_peak_height[i]))
  observed_peak_week[i] <- data$season_week[observed_peak_week_ind]
}

###################################################################
for(i in 1: length(model_name)){
  target <- result_table0[[i]]
  peak_week_log_score[i] <- -mean(target$peak_week_log_score)
  peak_height_log_score[i] <- -mean(target$peak_height_log_score)
  peak_week_log_score_max[i] <- max(- target$peak_week_log_score)
  peak_height_log_score_max[i] <- max(- target$peak_height_log_score)
  
  sum_peak_week_log_score_BP <- 0
  sum_peak_height_log_score_BP <- 0
  peak_week_log_score_BP_max[i] <- 0
  peak_height_log_score_BP_max[i] <- 0
  for(j in 1 : length(analysis_seasons)){
    season <- analysis_seasons[j]
    # index for before peak predictions
    BP_idx <- target$analysis_time_season == season & (target$analysis_time_season_week < observed_peak_week[j])
    sum_peak_week_log_score_BP <- sum(target$peak_week_log_score[BP_idx]) + sum_peak_week_log_score_BP
    sum_peak_height_log_score_BP <- sum(target$peak_height_log_score[BP_idx]) + sum_peak_height_log_score_BP
    max_week_log_scorej <- max(-target$peak_week_log_score[BP_idx]) 
    max_height_log_scorej <- max(-target$peak_height_log_score[BP_idx]) 
    if(max_week_log_scorej > peak_week_log_score_BP_max[i]) peak_week_log_score_BP_max[i] <- max_week_log_scorej
    if(max_height_log_scorej > peak_height_log_score_BP_max[i]) peak_height_log_score_BP_max[i] <- max_height_log_scorej
  }
  # average log scores over all week before peak
  peak_week_log_score_BP[i] <- -sum_peak_week_log_score_BP/(sum(observed_peak_week) - length(observed_peak_week))
  peak_height_log_score_BP[i] <- -sum_peak_height_log_score_BP/(sum(observed_peak_week) - length(observed_peak_week))
}

# Equal Bin
# all possible peak weeks: 33 week each flu season, season week 10 - 42

model_name[length(model_name) + 1] <- "Equal bin"
peak_week_log_score[length(model_name)] <- -log(1/33)
peak_week_log_score_max[length(model_name)] <- -log(1/33)
peak_week_log_score_BP[length(model_name)] <- -log(1/33)
peak_week_log_score_BP_max[length(model_name)] <- -log(1/33)


peak_height_log_score[length(model_name)] <- -log(1/27)
peak_height_log_score_max[length(model_name)] <- -log(1/27)
peak_height_log_score_BP[length(model_name)] <- -log(1/27)
peak_height_log_score_BP_max[length(model_name)] <- -log(1/27)

Subset <- c("All weeks",
            rep(" ",(length(model_name) - 1)),
            "Before peak",
            rep(" ",(length(model_name) - 1)))


peak_height_log_score <- 
  paste0(formatC(round(peak_height_log_score, digits = 2), format='f', digits=2 ), "(",rank(peak_height_log_score), ")")

peak_height_log_score_BP <- 
  paste0(formatC(round(peak_height_log_score_BP, digits = 2), format='f', digits=2 ), "(",rank(peak_height_log_score_BP), ")")

peak_height_log_score_max <- 
  paste0(formatC(round(peak_height_log_score_max, digits = 2), format='f', digits=2 ), "(",rank(peak_height_log_score_max), ")")
peak_height_log_score_BP_max <- 
  paste0(formatC(round(peak_height_log_score_BP_max, digits = 2), format='f', digits=2 ), "(",rank(peak_height_log_score_BP_max), ")")

peak_week_log_score <- 
  paste0(formatC(round(peak_week_log_score, digits = 2), format='f', digits=2 ), "(",rank(peak_week_log_score), ")")
peak_week_log_score_BP <- 
  paste0(formatC(round(peak_week_log_score_BP, digits = 2), format='f', digits=2 ), "(",rank(peak_week_log_score_BP), ")")

peak_week_log_score_max <- 
  paste0(formatC(round(peak_week_log_score_max, digits = 2), format='f', digits=2 ), "(",rank(peak_week_log_score_max), ")")
peak_week_log_score_BP_max <- 
  paste0(formatC(round(peak_week_log_score_BP_max, digits = 2), format='f', digits=2 ), "(",rank(peak_week_log_score_BP_max), ")")

res_peak <- data.frame(Model = c(model_name,model_name),
                       Subset = Subset,
                       Peak.intensity__LS = c(peak_height_log_score,
                                              peak_height_log_score_BP),
                       Peak.intensity__maxLS = c(peak_height_log_score_max,
                                                 peak_height_log_score_BP_max),
                       Peak.timing__LS = c(peak_week_log_score,peak_week_log_score_BP),
                       Peak.timing__maxLS = c(peak_week_log_score_max,peak_week_log_score_BP_max))
colnames(res_peak)[3] <- "Peak intensity__LS"
colnames(res_peak)[4] <- "Peak intensity__maxLS"

colnames(res_peak)[5] <- "Peak timing__LS"
colnames(res_peak)[6] <- "Peak timing__maxLS"
saveRDS(res_peak, file = here("./Results/Peak/res_peak.rds"))

