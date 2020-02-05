
RNGkind(sample.kind="Rounding")
set.seed(20200202)

library(here)
library(surveillance)

peak_arima <- readRDS(file = here("./Results/Peak/peak-week-arima.rds"))
peak_sarima <- readRDS(file = here("./Results/Peak/peak-week-sarima.rds"))
peak_Beta <- readRDS(file = here("./Results/Peak/peak-week-beta3-4.rds"))
peak_Beta_lags <- readRDS(file = here("./Results/Peak/peak-week-beta-lags.rds"))

# peak timing
peak_prophet <- readRDS(file = here("./Results/Peak/peak-week-prophet.rds"))
peak_naive <- readRDS(file = here("./Results/Peak/peak-week-naive.rds"))
peak_KCDE <- readRDS(file = here("./Results/Peak/peak-week-KCDE.rds"))

Per_Beta_T <- permutationTest(peak_Beta$peak_week_log_score, 
                            peak_KCDE$peak_week_log_score, 
                            nPermutation = 9999)
Per_Beta4_T <- permutationTest(peak_Beta_lags$peak_week_log_score, 
                               peak_KCDE$peak_week_log_score, 
                            nPermutation = 9999)

Per_ARIMA_T <- permutationTest(peak_arima$peak_week_log_score, 
                               peak_KCDE$peak_week_log_score, 
                             nPermutation = 9999)

Per_SARIMA_T <- permutationTest(peak_sarima$peak_week_log_score, 
                                peak_KCDE$peak_week_log_score, 
                              nPermutation = 9999)

Per_Prophet_T <- permutationTest(peak_prophet$peak_week_log_score, 
                                 peak_KCDE$peak_week_log_score, 
                               nPermutation = 9999)

Per_Naive_T <- permutationTest(peak_naive$peak_week_log_score, 
                               peak_KCDE$peak_week_log_score, 
                             nPermutation = 9999)

Per_EB_T <- permutationTest(rep(-log(1/33), length(peak_KCDE$peak_week_log_score)), 
                               peak_KCDE$peak_week_log_score, 
                               nPermutation = 9999)

Per_list_T <-   as.numeric(c(Per_Beta_T$pVal.permut, 
                             Per_Beta4_T$pVal.permut, 
                            NA,
                            Per_ARIMA_T$pVal.permut,
                            Per_SARIMA_T$pVal.permut,
                            Per_Prophet_T$pVal.permut,
                            Per_Naive_T$pVal.permut,
                            Per_EB_T$pVal.permut))

#-----------------------------------------------------------------------------
# peak timing before peak
# indicator for weeks before peak 
analysis_seasons <- c("2014/2015","2015/2016", "2016/2017", "2017/2018")
observed_peak_week <- numeric(length(analysis_seasons))
observed_peak_height <- numeric(length(analysis_seasons))

load(file = here("./Data/data_holidays.RData"))


first_analysis_time_season_week <- 10 # == week 40 of year
last_analysis_time_season_week <- 41 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead

prediction_target_var <- "weighted_ili"

BP_idx <- numeric(dim(peak_arima)[1])
for(i in 1 : length(analysis_seasons)) {
  analysis_time_season <- analysis_seasons[i]
  observed_peak_height[i] <- max(data[data$season == analysis_time_season, prediction_target_var])
  observed_peak_week_ind <- which((data$season == analysis_time_season) &
                                    (data[, prediction_target_var] == observed_peak_height[i]))
  observed_peak_week[i] <- data$season_week[observed_peak_week_ind]
  BP_idx <- BP_idx + 
    (peak_arima$analysis_time_season == analysis_time_season) & 
    (peak_arima$analysis_time_season_week < observed_peak_week[i])
}

Per_Beta_T_BP <- permutationTest(peak_Beta$peak_week_log_score[BP_idx], 
                              peak_KCDE$peak_week_log_score[BP_idx], 
                              nPermutation = 9999)
Per_Beta4_T_BP <- permutationTest(peak_Beta_lags$peak_week_log_score[BP_idx], 
                               peak_KCDE$peak_week_log_score[BP_idx], 
                               nPermutation = 9999)

Per_ARIMA_T_BP <- permutationTest(peak_arima$peak_week_log_score[BP_idx], 
                               peak_KCDE$peak_week_log_score[BP_idx], 
                               nPermutation = 9999)

Per_SARIMA_T_BP <- permutationTest(peak_sarima$peak_week_log_score[BP_idx], 
                                peak_KCDE$peak_week_log_score[BP_idx], 
                                nPermutation = 9999)

Per_Prophet_T_BP <- permutationTest(peak_prophet$peak_week_log_score[BP_idx], 
                                 peak_KCDE$peak_week_log_score[BP_idx], 
                                 nPermutation = 9999)

Per_Naive_T_BP <- permutationTest(peak_naive$peak_week_log_score[BP_idx], 
                               peak_KCDE$peak_week_log_score[BP_idx], 
                               nPermutation = 9999)

Per_EB_T_BP <- permutationTest(rep(-log(1/33), length(peak_KCDE$peak_week_log_score[BP_idx])), 
                                  peak_KCDE$peak_week_log_score[BP_idx], 
                                  nPermutation = 9999)


Per_list_T_BP <-   as.numeric(c(Per_Beta_T_BP$pVal.permut, 
                             Per_Beta4_T_BP$pVal.permut, 
                             NA,
                             Per_ARIMA_T_BP$pVal.permut,
                             Per_SARIMA_T_BP$pVal.permut,
                             Per_Prophet_T_BP$pVal.permut,
                             Per_Naive_T_BP$pVal.permut,
                             Per_EB_T_BP$pVal.permut))

per_T <- c(Per_list_T, Per_list_T_BP)
###################################################################
# peak intensity
Per_Beta_I <- permutationTest(peak_Beta$peak_height_log_score, 
                              peak_KCDE$peak_height_log_score, 
                              nPermutation = 9999)
Per_Beta4_I <- permutationTest(peak_Beta_lags$peak_height_log_score, 
                               peak_KCDE$peak_height_log_score, 
                               nPermutation = 9999)

Per_ARIMA_I <- permutationTest(peak_arima$peak_height_log_score, 
                               peak_KCDE$peak_height_log_score, 
                               nPermutation = 9999)

Per_SARIMA_I <- permutationTest(peak_sarima$peak_height_log_score, 
                                peak_KCDE$peak_height_log_score, 
                                nPermutation = 9999)

Per_Prophet_I <- permutationTest(peak_prophet$peak_height_log_score, 
                                 peak_KCDE$peak_height_log_score, 
                                 nPermutation = 9999)

Per_Naive_I <- permutationTest(peak_naive$peak_height_log_score, 
                               peak_KCDE$peak_height_log_score, 
                               nPermutation = 9999)

Per_EB_I <- permutationTest(rep(-log(1/27), length(peak_KCDE$peak_height_log_score)), 
                            peak_KCDE$peak_height_log_score, 
                            nPermutation = 9999)

Per_list_I <-   as.numeric(c(Per_Beta_I$pVal.permut, 
                             Per_Beta4_I$pVal.permut, 
                             NA,
                             Per_ARIMA_I$pVal.permut,
                             Per_SARIMA_I$pVal.permut,
                             Per_Prophet_I$pVal.permut,
                             Per_Naive_I$pVal.permut,
                             Per_EB_I$pVal.permut))
#----------------------------------------------
# peak intensity before peak
Per_Beta_I_BP <- permutationTest(peak_Beta$peak_height_log_score[BP_idx], 
                                 peak_KCDE$peak_height_log_score[BP_idx], 
                                 nPermutation = 9999)
Per_Beta4_I_BP <- permutationTest(peak_Beta_lags$peak_height_log_score[BP_idx], 
                                  peak_KCDE$peak_height_log_score[BP_idx], 
                                  nPermutation = 9999)

Per_ARIMA_I_BP <- permutationTest(peak_arima$peak_height_log_score[BP_idx], 
                                  peak_KCDE$peak_height_log_score[BP_idx], 
                                  nPermutation = 9999)

Per_SARIMA_I_BP <- permutationTest(peak_sarima$peak_height_log_score[BP_idx], 
                                   peak_KCDE$peak_height_log_score[BP_idx], 
                                   nPermutation = 9999)

Per_Prophet_I_BP <- permutationTest(peak_prophet$peak_height_log_score[BP_idx], 
                                    peak_KCDE$peak_height_log_score[BP_idx], 
                                    nPermutation = 9999)

Per_Naive_I_BP <- permutationTest(peak_naive$peak_height_log_score[BP_idx], 
                                  peak_KCDE$peak_height_log_score[BP_idx], 
                                  nPermutation = 9999)

Per_EB_I_BP <- permutationTest(rep(-log(1/27), length(peak_KCDE$peak_height_log_score[BP_idx])), 
                               peak_KCDE$peak_height_log_score[BP_idx], 
                               nPermutation = 9999)


Per_list_I_BP <-   as.numeric(c(Per_Beta_I_BP$pVal.permut, 
                                Per_Beta4_I_BP$pVal.permut, 
                                NA,
                                Per_ARIMA_I_BP$pVal.permut,
                                Per_SARIMA_I_BP$pVal.permut,
                                Per_Prophet_I_BP$pVal.permut,
                                Per_Naive_I_BP$pVal.permut,
                                Per_EB_I_BP$pVal.permut))

per_I <- c(Per_list_I, Per_list_I_BP)

saveRDS(list(pvalue_I = per_I, pvalue_T = per_T), file = here("./Results/LT_pvalue.rds"))
