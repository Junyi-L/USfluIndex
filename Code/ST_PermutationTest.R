
RNGkind(sample.kind="Rounding")
library(xtable)
library(here)
library(surveillance)

kcde_predictions_ph_1 <- readRDS(here("./Results/Forecast_ph1-4/KCDEresults/kcde-predictions-ph_1.rds"))
kcde_predictions_ph_2 <- readRDS(here("./Results/Forecast_ph1-4/KCDEresults/kcde-predictions-ph_2.rds"))
kcde_predictions_ph_3 <- readRDS(here("./Results/Forecast_ph1-4/KCDEresults/kcde-predictions-ph_3.rds"))
kcde_predictions_ph_4 <- readRDS(here("./Results/Forecast_ph1-4/KCDEresults/kcde-predictions-ph_4.rds"))
kcde_prediction <- rbind(kcde_predictions_ph_1,kcde_predictions_ph_2,kcde_predictions_ph_3,kcde_predictions_ph_4)


ARIMA <- readRDS(file = here("./Results/Forecast_ph1-4/ArimaResults.rds"))
SARIMA <- readRDS(file = here("./Results/Forecast_ph1-4/SARIMAResults.rds"))

Beta <- readRDS(file = here("./Results/Forecast_ph1-4/Beta.rds"))
Beta_lag4 <- readRDS(file = here("./Results/Forecast_ph1-4/Beta_lag.rds"))


ProphetResults <- readRDS(file = here("./Results/Forecast_ph1-4/ProphetResults.rds"))
NaiveResults <- readRDS(file = here("./Results/Forecast_ph1-4/NaiveResults.rds"))




Per_Beta <- permutationTest(Beta$result$log_score, 
                            Beta_lag4$result$log_score, 
                            nPermutation = 9999)
Per_KCDE <- permutationTest(kcde_prediction$log_score, 
                            Beta_lag4$result$log_score, 
                            nPermutation = 9999)

Per_ARIMA <- permutationTest(ARIMA$result$log_score, 
                             Beta_lag4$result$log_score, 
                             nPermutation = 9999)

Per_SARIMA <- permutationTest(SARIMA$result$log_score, 
                              Beta_lag4$result$log_score, 
                              nPermutation = 9999)

Per_Prophet <- permutationTest(ProphetResults$result$log_score, 
                               Beta_lag4$result$log_score, 
                               nPermutation = 9999)

Per_Naive <- permutationTest(NaiveResults$result$log_score, 
                             Beta_lag4$result$log_score, 
                             nPermutation = 9999)

Per_list1 <-   as.numeric(c(Per_Beta$pVal.permut, 
                           NA, 
                           Per_KCDE$pVal.permut,
                           Per_ARIMA$pVal.permut,
                           Per_SARIMA$pVal.permut,
                           Per_Prophet$pVal.permut,
                           Per_Naive$pVal.permut))
#--------------------------------------------------

Per_sub_Beta <- permutationTest(Beta$result$log_score[Beta$result$data_set.InSeason == TRUE], 
                                Beta_lag4$result$log_score[Beta_lag4$result$data_set.InSeason == TRUE], 
                                nPermutation = 9999)
Per_sub_KCDE <- permutationTest(kcde_prediction$log_score[kcde_prediction$data_set.InSeason == TRUE], 
                                Beta_lag4$result$log_score[Beta_lag4$result$data_set.InSeason == TRUE], 
                                nPermutation = 9999)

Per_sub_ARIMA <- permutationTest(ARIMA$result$log_score[ARIMA$result$data_set.InSeason == TRUE], 
                                 Beta_lag4$result$log_score[Beta_lag4$result$data_set.InSeason == TRUE], 
                                 nPermutation = 9999)

Per_sub_SARIMA <- permutationTest(SARIMA$result$log_score[SARIMA$result$data_set.InSeason == TRUE], 
                                  Beta_lag4$result$log_score[Beta_lag4$result$data_set.InSeason == TRUE], 
                                  nPermutation = 9999)

Per_sub_Prophet <- permutationTest(ProphetResults$result$log_score[ProphetResults$result$data_set.InSeason == TRUE], 
                                   Beta_lag4$result$log_score[Beta_lag4$result$data_set.InSeason == TRUE], 
                                   nPermutation = 9999)

Per_sub_Naive <- permutationTest(NaiveResults$result$log_score[NaiveResults$result$data_set.InSeason == TRUE], 
                                 Beta_lag4$result$log_score[Beta_lag4$result$data_set.InSeason == TRUE], 
                                 nPermutation = 9999)
Per_list2 <-   as.numeric(c(Per_sub_Beta$pVal.permut, 
                           NA, 
                           Per_sub_KCDE$pVal.permut,
                           Per_sub_ARIMA$pVal.permut,
                           Per_sub_SARIMA$pVal.permut,
                           Per_sub_Prophet$pVal.permut,
                           Per_sub_Naive$pVal.permut))
pvalue <- c(Per_list1,Per_list2)
saveRDS(pvalue, file = here("./Results/ST_pvalue.rds"))

