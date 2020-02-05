
library(here)

kcde_predictions_ph_1 <- readRDS(here("./Results/Forecast_ph1-4/KCDEresults/kcde-predictions-ph_1.rds"))
kcde_predictions_ph_2 <- readRDS(here("./Results/Forecast_ph1-4/KCDEresults/kcde-predictions-ph_2.rds"))
kcde_predictions_ph_3 <- readRDS(here("./Results/Forecast_ph1-4/KCDEresults/kcde-predictions-ph_3.rds"))
kcde_predictions_ph_4 <- readRDS(here("./Results/Forecast_ph1-4/KCDEresults/kcde-predictions-ph_4.rds"))
kcde_prediction <- rbind(kcde_predictions_ph_1,kcde_predictions_ph_2,kcde_predictions_ph_3,kcde_predictions_ph_4)


kcde_time_fit_ph_1 <-
  readRDS(here("./Results/Forecast_ph1-4/KCDEresults/fit_time-ili_national-prediction_horizon_1-max_lag_1-max_seasonal_lag_0-filtering_FALSE-differencing_FALSE-seasonality_TRUE-bw_parameterization_full.rds"))
kcde_time_fit_ph_2 <-
  readRDS(here("./Results/Forecast_ph1-4/KCDEresults/fit_time-ili_national-prediction_horizon_1-max_lag_1-max_seasonal_lag_0-filtering_FALSE-differencing_FALSE-seasonality_TRUE-bw_parameterization_full.rds"))
kcde_time_fit_ph_3 <-
  readRDS(here("./Results/Forecast_ph1-4/KCDEresults/fit_time-ili_national-prediction_horizon_1-max_lag_1-max_seasonal_lag_0-filtering_FALSE-differencing_FALSE-seasonality_TRUE-bw_parameterization_full.rds"))
kcde_time_fit_ph_4 <-
  readRDS(here("./Results/Forecast_ph1-4/KCDEresults/fit_time-ili_national-prediction_horizon_1-max_lag_1-max_seasonal_lag_0-filtering_FALSE-differencing_FALSE-seasonality_TRUE-bw_parameterization_full.rds"))
kcde_time_fit <- sum(kcde_time_fit_ph_1[3],kcde_time_fit_ph_2[3],kcde_time_fit_ph_3[3],kcde_time_fit_ph_4[4])

kcde_predictions_time_ph_1 <- readRDS(here("./Results/Forecast_ph1-4/KCDEresults/kcde-predictions-time-ph_1.rds"))
kcde_predictions_time_ph_2 <- readRDS(here("./Results/Forecast_ph1-4/KCDEresults/kcde-predictions-time-ph_2.rds"))
kcde_predictions_time_ph_3 <- readRDS(here("./Results/Forecast_ph1-4/KCDEresults/kcde-predictions-time-ph_3.rds"))
kcde_predictions_time_ph_4 <- readRDS(here("./Results/Forecast_ph1-4/KCDEresults/kcde-predictions-time-ph_4.rds"))
kcde_time_prediction <- sum(kcde_predictions_time_ph_1[3],kcde_predictions_time_ph_2[3],
                            kcde_predictions_time_ph_3[3],kcde_predictions_time_ph_4[3])

kcde_time <- kcde_time_fit + kcde_time_prediction

ARIMA <- readRDS(file = here("./Results/Forecast_ph1-4/ArimaResults.rds"))
SARIMA <- readRDS(file = here("./Results/Forecast_ph1-4/SARIMAResults.rds"))

Beta <- readRDS(file = here("./Results/Forecast_ph1-4/Beta.rds"))
Beta_lag4 <- readRDS(file = here("./Results/Forecast_ph1-4/Beta_lag.rds"))


ProphetResults <- readRDS(file = here("./Results/Forecast_ph1-4/ProphetResults.rds"))
NaiveResults <- readRDS(file = here("./Results/Forecast_ph1-4/NaiveResults.rds"))

Model_name <- c("Beta(1)",
                "Beta(4)",
                "KCDE",
                "ARIMA",
                "SARIMA",
                "Prophet",
                "Naive")
result_table <- list( Beta$result,
                      Beta_lag4$result,
                      kcde_prediction,
                      ARIMA$result,
                      SARIMA$result,
                      ProphetResults$result,
                      NaiveResults$result)
total_time <- c(Beta$run_time[3],
                Beta_lag4$run_time[3],
                kcde_time,
                ARIMA$run_time[3],
                SARIMA$run_time[3],
                ProphetResults$run_time[3],
                NaiveResults$run_time[3])

log_score_ph1 <- numeric()
log_score_ph2 <- numeric()
log_score_ph3 <- numeric()
log_score_ph4 <- numeric()
max_log_score <- numeric()
log_score <- numeric()
DSS <- numeric()
DSS_ph1 <- numeric()
DSS_ph2 <- numeric()
DSS_ph3 <- numeric()
DSS_ph4 <- numeric()
AE <- numeric()

log_score_ph1_sub <- numeric()
log_score_ph2_sub <- numeric()
log_score_ph3_sub <- numeric()
log_score_ph4_sub <- numeric()
max_log_score_sub <- numeric()
log_score_sub <- numeric()
DSS_sub <- numeric()
DSS_ph1_sub <- numeric()
DSS_ph2_sub <- numeric()
DSS_ph3_sub <- numeric()
DSS_ph4_sub <- numeric()
AE_sub <- numeric()

for(i in 1:length(result_table)){
  target <- result_table[[i]]
  log_score_ph1[i] <- - mean(target$log_score[target$prediction_horizon == 1])
  log_score_ph2[i] <- - mean(target$log_score[target$prediction_horizon == 2])
  log_score_ph3[i] <- - mean(target$log_score[target$prediction_horizon == 3])
  log_score_ph4[i] <- - mean(target$log_score[target$prediction_horizon == 4])
  log_score[i] <- - mean(target$log_score)
  max_log_score[i] <- -min(target$log_score)
  DSS_ph1[i] <- mean(target$DS_score[target$prediction_horizon == 1])
  DSS_ph2[i] <- mean(target$DS_score[target$prediction_horizon == 2])
  DSS_ph3[i] <- mean(target$DS_score[target$prediction_horizon == 3])
  DSS_ph4[i] <- mean(target$DS_score[target$prediction_horizon == 4])
  DSS[i] <- mean(target$DS_score)
  AE[i] <- mean(target$AE)

  log_score_ph1_sub[i] <- - mean(target$log_score[target$prediction_horizon == 1 & target$data_set.InSeason == TRUE])
  log_score_ph2_sub[i] <- - mean(target$log_score[target$prediction_horizon == 2 & target$data_set.InSeason == TRUE])
  log_score_ph3_sub[i] <- - mean(target$log_score[target$prediction_horizon == 3 & target$data_set.InSeason == TRUE])
  log_score_ph4_sub[i] <- - mean(target$log_score[target$prediction_horizon == 4 & target$data_set.InSeason == TRUE])
  log_score_sub[i] <- - mean(target$log_score[target$data_set.InSeason == TRUE])
  max_log_score_sub[i] <- -min(target$log_score[target$data_set.InSeason == TRUE])
  DSS_ph1_sub[i] <- mean(target$DS_score[target$prediction_horizon == 1 & target$data_set.InSeason == TRUE])
  DSS_ph2_sub[i] <- mean(target$DS_score[target$prediction_horizon == 2 & target$data_set.InSeason == TRUE ])
  DSS_ph3_sub[i] <- mean(target$DS_score[target$prediction_horizon == 3 & target$data_set.InSeason == TRUE ])
  DSS_ph4_sub[i] <- mean(target$DS_score[target$prediction_horizon == 4 & target$data_set.InSeason == TRUE])
  DSS_sub[i] <- mean(target$DS_score[target$data_set.InSeason == TRUE])
  AE_sub[i] <- mean(target$AE[target$data_set.InSeason == TRUE])
}

Subset <- c("All weeks",
            rep(" ",(length(Model_name) - 1)),
            "weeks 40--20",
            rep(" ",(length(Model_name) - 1)))

log_score_ph1 <-
  paste0(formatC(round(log_score_ph1, digits = 2), format='f', digits=2 ), " (",rank(log_score_ph1), ")")
log_score_ph1_sub <-
  paste0(formatC(round(log_score_ph1_sub, digits = 2), format='f', digits=2 ), " (",rank(log_score_ph1_sub), ")")
DSS_ph1 <-
  paste0(formatC(round(DSS_ph1, digits = 2), format='f', digits=2 ), " (",rank(DSS_ph1), ")")
DSS_ph1_sub <-
  paste0(formatC(round(DSS_ph1_sub, digits = 2), format='f', digits=2 ), " (",rank(DSS_ph1_sub), ")")

log_score_ph2 <-
  paste0(formatC(round(log_score_ph2, digits = 2), format='f', digits=2 ), " (",rank(log_score_ph2), ")")
log_score_ph2_sub <-
  paste0(formatC(round(log_score_ph2_sub, digits = 2), format='f', digits=2 ), " (",rank(log_score_ph2_sub), ")")
DSS_ph2 <-
  paste0(formatC(round(DSS_ph2, digits = 2), format='f', digits=2 ), " (",rank(DSS_ph2), ")")
DSS_ph2_sub <-
  paste0(formatC(round(DSS_ph2_sub, digits = 2), format='f', digits=2 ), " (",rank(DSS_ph2_sub), ")")

log_score_ph3 <-
  paste0(formatC(round(log_score_ph3, digits = 2), format='f', digits=2 ), " (",rank(log_score_ph3), ")")
log_score_ph3_sub <-
  paste0(formatC(round(log_score_ph3_sub, digits = 2), format='f', digits=2 ), " (",rank(log_score_ph3_sub), ")")
DSS_ph3 <-
  paste0(formatC(round(DSS_ph3, digits = 2), format='f', digits=2 ), " (",rank(DSS_ph3), ")")
DSS_ph3_sub <-
  paste0(formatC(round(DSS_ph3_sub, digits = 2), format='f', digits=2 ), " (",rank(DSS_ph3_sub), ")")

log_score_ph4 <-
  paste0(formatC(round(log_score_ph4, digits = 2), format='f', digits=2 ), " (",rank(log_score_ph4), ")")
log_score_ph4_sub <-
  paste0(formatC(round(log_score_ph4_sub, digits = 2), format='f', digits=2 ), " (",rank(log_score_ph4_sub), ")")
DSS_ph4 <-
  paste0(formatC(round(DSS_ph4, digits = 2), format='f', digits=2 ), " (",rank(DSS_ph4), ")")
DSS_ph4_sub <-
  paste0(formatC(round(DSS_ph4_sub, digits = 2), format='f', digits=2 ), " (",rank(DSS_ph4_sub), ")")



res_ph1_4 <- data.frame(Model = Model_name,
                        Subset = Subset,
                        ph1__LS = gsub("-","--",c(log_score_ph1, log_score_ph1_sub)),
                        ph1__DSS = gsub("-","--",c(DSS_ph1, DSS_ph1_sub)),
                        ph2__LS = gsub("-","--",c(log_score_ph2, log_score_ph2_sub)),
                        ph2__DSS = gsub("-","--",c(DSS_ph2, DSS_ph2_sub)),
                        ph3__LS = gsub("-","--",c(log_score_ph3, log_score_ph3_sub)),
                        ph3__DSS = gsub("-","--",c(DSS_ph3, DSS_ph3_sub)),
                        ph4__LS = gsub("-","--",c(log_score_ph4, log_score_ph4_sub)),
                        ph4__DSS = gsub("-","--",c(DSS_ph4, DSS_ph4_sub)))


saveRDS(res_ph1_4, file = here("./Results/Forecast_ph1-4/res_ph1_4.rds"))

npar <- c(19,20,28,16,3, 50,106)

log_score <-
  paste0(formatC(round(log_score, digits = 2), format='f', digits=2 ), " (",rank(log_score), ")")
log_score_sub <-
  paste0(formatC(round(log_score_sub, digits = 2), format='f', digits=2 ), " (",rank(log_score_sub), ")")

max_log_score <-
  paste0(formatC(round(max_log_score, digits = 2), format='f', digits=2 ), " (",rank(max_log_score), ")")
max_log_score_sub <-
  paste0(formatC(round(max_log_score_sub, digits = 2), format='f', digits=2 ), " (",rank(max_log_score_sub), ")")

DSS <-
  paste0(formatC(round(DSS, digits = 2), format='f', digits=2 ), " (",rank(DSS), ")")
DSS_sub <-
  paste0(formatC(round(DSS_sub, digits = 2), format='f', digits=2 ), " (",rank(DSS_sub), ")")

AE <-
  paste0(formatC(round(AE, digits = 2), format='f', digits=2 ), " (",rank(AE), ")")
AE_sub <-
  paste0(formatC(round(AE_sub, digits = 2), format='f', digits=2 ), " (",rank(AE_sub), ")")

total_time <-
  paste0(formatC(round(total_time/60, digits = 2), format='f', digits=2 ), " (",rank(total_time), ")")

npar <-
  paste0(formatC(round(npar, digits = 0), format='f', digits=0 ), " (",rank(npar), ")")

pvalue <- readRDS(file = here("./Results/ST_pvalue.rds"))

res_ph1_4_2 <- data.frame(Model = Model_name,
                          Subset = Subset,
                          LS = gsub("-","--",c(log_score, log_score_sub)),
                          pvalue = pvalue,
                          maxLS = gsub("-","--",c(max_log_score, max_log_score_sub)),
                          DSS = gsub("-","--",c(DSS, DSS_sub)),
                          AE = gsub("-","--",c(AE, AE_sub)),
                          Time = c(total_time, rep(NA, length(Model_name))),
                          npar = c(npar, rep(NA, length(Model_name))))

saveRDS(res_ph1_4_2, file = here("./Results/Forecast_ph1-4/res_ph1_4_2.rds"))
# ################# fanplots
# # General setting
# library("surveillance")
# library("HIDDA.forecasting")
# pal <- colorRampPalette(c("darkgreen", "gray93"))
# obs_data <- data[data$train == FALSE,]
# width <- 14
# height <- 7
# ylimi <- c(0,9)
#
# # KCDE
# kcde_quantile <- readRDS("./Forecast_ph1-4/KCDEresults/quantile_matrix.rds")
# means <- kcde_predictions_ph_1$pt_pred
# y <- obs_data$weighted_ili
# probs <- (1:99)/100
# quantiles <- kcde_quantile
#
# pdf(file = "./Forecast_ph1-4/kcde_fan.pdf",width = width,height = height)
# osaplot(
#   quantiles = quantiles, probs = 1:99/100,
#   observed = y, scores = cbind( logs = -kcde_predictions_ph_1$log_score,
#                                 DSS = kcde_predictions_ph_1$DS_score),
#   ylim = ylimi,xlab = "Week",ylab = "wILI", cex.lab=2,  cex.main=2, cex.sub=2, cex.axis=2,
#   fan.args = list(ln = c(0.1,0.9), rlab = NULL),
#   scores.args = list(ylim = c(-6, 3), cex.lab=2,  cex.main=2, cex.sub=2, cex.axis=2, las = TRUE),
#   legend.args = list(cex = 0.8, legend = c("LS","DSS")),
#   key.args = list(space = 1),
#   las = TRUE
# )
# abline(v = 2.2)
# dev.off()
#
# # ARIMA
# ARIMA_quantile <- ARIMA$quantile_matrix
# means <- ARIMA$result$pt_pred[ARIMA$result$prediction_horizon == 1]
# y <- obs_data$weighted_ili
# probs <- (1:99)/100
# quantiles <- ARIMA_quantile
#
# pdf(file = "./Forecast_ph1-4/ARIMA_fan.pdf",width = 14,height = 7)
# osaplot(
#   quantiles = quantiles, probs = 1:99/100,
#   observed = y, scores = cbind( logs = -ARIMA$result$log_score[ARIMA$result$prediction_horizon ==1],
#                                 DSS = ARIMA$result$DS_score[ARIMA$result$prediction_horizon ==1]),
#   ylim = ylimi,xlab = "Week",ylab = "wILI", cex.lab=2,  cex.main=2, cex.sub=2, cex.axis=2,
#   fan.args = list(ln = c(0.1,0.9), rlab = NULL),
#   scores.args = list(ylim = c(-6, 3), cex.lab=2,  cex.main=2, cex.sub=2, cex.axis=2, las = TRUE),
#   legend.args = list(cex = 0.8, legend = c("LS","DSS")),
#   key.args = list(space = 1),
#   las = TRUE
# )
#
# dev.off()
#
# # Beta
# Beta_quantile <- Beta$quantile_matrix
# means <- Beta$result$pt_pred[Beta$result$prediction_horizon == 1]
# y <- obs_data$weighted_ili
# probs <- (1:99)/100
# quantiles <- Beta_quantile
#
# pdf(file = "./Forecast_ph1-4/Beta_fan.pdf",width = 14,height = 7)
# par(mar =  c(5.1, 5.1, 4.1, 2.1))
# osaplot(
#   quantiles = quantiles, probs = 1:99/100,
#   observed = y, scores = cbind( logs = -Beta$result$log_score[Beta$result$prediction_horizon == 1],
#                                 DSS = Beta$result$DS_score[Beta$result$prediction_horizon == 1]),
#   xlab = "Week", ylim = ylimi,ylab = "wILI",cex.lab=2,  cex.main=2, cex.sub=2, cex.axis=2,
#   fan.args = list(ln = c(0.1,0.9), rlab = NULL),
#   scores.args = list(ylim = c(-6, 3), cex.lab=2,  cex.main=2, cex.sub=2, cex.axis=2, las = TRUE),
#   legend.args = list(cex = 0.8, legend = c("LS","DSS")),
#   key.args = list(space = 1),
#   las = TRUE
# )
#
# dev.off()
#
# # Beta
# Beta_quantile <- Beta_lag4$quantile_matrix
# means <- Beta_lag4$result$pt_pred[Beta_lag4$result$prediction_horizon == 1]
# y <- obs_data$weighted_ili
# probs <- (1:99)/100
# quantiles <- Beta_quantile
# ylimi <- c(0,9)
# pdf(file = "./Forecast_ph1-4/Beta_lag_fan.pdf",width = 14,height = 7)
# par(mar =  c(5.1, 5.1, 4.1, 2.1))
# osaplot(
#   quantiles = quantiles, probs = 1:99/100,
#   observed = y, scores = cbind( logs = -Beta$result$log_score[Beta$result$prediction_horizon == 1],
#                                 DSS = Beta$result$DS_score[Beta$result$prediction_horizon == 1]),
#   xlab = "Week", ylim = ylimi,ylab = "wILI",cex.lab=2,  cex.main=2, cex.sub=2, cex.axis=2,
#   fan.args = list(ln = c(0.1,0.9), rlab = NULL),
#   scores.args = list(ylim = c(-6, 3), cex.lab=2,  cex.main=2, cex.sub=2, cex.axis=2, las = TRUE),
#   legend.args = list(cex = 0.8, legend = c("LS","DSS")),
#   key.args = list(space = 1),
#   las = TRUE
# )
#
# dev.off()
#
