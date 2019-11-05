# USfluIndex

## Description
This repository contains code and data to reproduce the model comparison result in the paper "Forecasting Flu Activity in the United States: Benchmarking an Endemic-Epidemic Beta Model". 

Seven models are compared in this paper,

* Beta model with 1 lag, using [**betareg**](https://CRAN.R-project.org/package=betareg)`::betareg()`,

* Beta model with more lags, using [**betareg**](https://CRAN.R-project.org/package=betareg)`::betareg()`,

* Harmonic regression with ARIMA errors (ARIMA), using [**forecast**](https://CRAN.R-project.org/package=forecast)`::auto.arima()`,

* SARIMA, using [**forecast**](https://CRAN.R-project.org/package=forecast)`::auto.arima()`,

* Prophet, using [**prophet**](https://CRAN.R-project.org/package=prophet)`::prophet()`,

* KCDE, using [**kcde**](https://github.com/reichlab/kcde)`::kcde()` and code in [**article-disease-pred-with-kcde**](https://github.com/reichlab/article-disease-pred-with-kcde),

* A naive approach.

Two kinds of prediction targets are scored,

* Short-term forecast for prediction horizon 1-4,

* Peak prediction, including peak intensity and peak timing prediction. 

## This repository is organized as follows:

* Data/ contains code to load data and the loaded data.

    * Data/load_data.R loads wILI data through the R package [**cdcfluview**](https://cran.r-project.org/web/packages/cdcfluview/index.html). Details about this wILI data please see http://www.cdc.gov/flu/weekly/. Note that when ILINet members provide revisions or backfill reports for past weeks, the wILI data will be updated accordingly. This means that data will have samll difference when loading at different time points. The data used for comparison in this paper are loaded on March 11 2019. The loaded wILI data are saved in Data/usflu.RData.

    * Data/load_data_holiday.R adds necessary columns in loaded wILI data, including holidays and sin, cos terms and saves the data in Data/data_holidays.RData.

* Code/ contains code used to estimate and make forecasts.

    * Code/Short-term/ contains code to estimate and make short-term forecast.
    
        * Code/Short-term/Beta.R selects a Beta model with 1 lag based on AICc, estimates parameters and perform short-term forecasts, and scores of forecast are calculated are saved in Results/Forecast_ph1-4/Beta.rds.
        * Code/Short-term/Beta_lags.R selects a Beta model with more lags based on AICc, estimates parameters and performs short-term forecasts, and scores of forecast are calculated are saved in Results/Forecast_ph1-4/Beta_lag4.rds.
        * Code/Short-term/ARIMA_model_choice.R selects an ARIMA model using training data based on AICc. Code/Short-term/ARIMA.R performs estimation and short-term forecasts, and scores of forecast are calculated are saved in Results/Forecast_ph1-4/ArimaResults.rds.
        * Code/Short-term/SARIMA.R selects a SARIMA model and performs estimation and short-term forecasts, and saves scores of forecast in Results/Forecast_ph1-4/SARIMAResults.rds.
        * Code/Short-term/naive.R performs estimation and short-term forecasts of the naive approach, and saves scores of forecast in Results/Forecast_ph1-4/NaiveResults.rds.
        * Code/Short-term/prophet.R performs estimation and short-term forecasts of the Prophet model, and saves scores of forecast in Results/Forecast_ph1-4/ProphetResults.rds.
        * Code/Short-term/KCDE_step_estimation.R estimates parameters for the KCDE model, and saves them in Results/Forecast_ph1-4/KCDEResults. Code/Short-term/KCDE_step_prediction.R performs short-term forecasts, and saves scores of forecast in Results/Forecast_ph1-4/KCDEResults/.
        
    * Code/Peak/ contains code to predict the peak timing and the peak intensity.
    
        * Code/Peak/Beta_peak.R predicts the peak timing and intensity using the selected Beta model with 1 lag and saves the prediction and scores in Results/Peak/peak-week-beta3-4.rds. 
        * Code/Peak/Beta_peak_lags.R predicts the peak timing and intensity using the selected Beta model with more lags and saves the prediction and scores in Results/Peak/peak-week-beta-lags.rds.
        * Code/Peak/ARIMA_peak.R predicts the peak timing and intensity using the selected ARIMA model and saves the prediction and scores in Results/Peak/peak-week-arima.rds. 
        * Code/Peak/SARIMA_peak.R predicts the peak timing and intensity using the selected SARIMA model and saves the prediction and scores in Results/Peak/peak-week-sarima.rds. 
        * Code/Peak/Prophet_peak.R predicts the peak timing and intensity using the Prophet model and saves the prediction and scores in Results/Peak/peak-week-prophet.rds.
        * Code/Peak/KCDE_copula_estimation.R estimates the corresponding copulas of the KCDE model and saves the estimated parameters in Code/USfluIndex/Results/Peak/copula-estimation-results/. Code/Peak/Prophet_peak.R predicts the peak timing and intensity using the fitted copulas and estimated KCDE, and saves the prediction and scores in Results/Peak/peak-week-KCDE.rds.
        * Code/Peak/naive_peak.R predicts the peak timing and intensity using the naive approach and saves the prediction and scores in Results/Peak/peak-week-naive.rds.
  
    * Code/Beta_forecast_p.R defines a function `Beta_Forecast_p` to simulate trajectories from a given Beta model for a given prediction horizon. This function is used for both short-term forecast and peak prediction of the Beta model. A function `cal_log_score` is defined here, which calculates the log score of short-term forecast, based on the averaged predictive distribution of each trajectory. 
    
    * Code/summary_ph1-4.R summarizes the short-term forecast and saves the summary tables in Results/Forecast_ph1-4/res_ph1_4.rds and Results/Forecast_ph1-4/res_ph1_4_2.rds.
    
    * Code/summary_peak.R summarizes the peak prediction and saves the summary table in Results/Peak/res_peak.rds.
    
    * Code/plots_paper.R generates some plots for the paper and saves them in Results/plots/.
    
* tex/ contains tex/cdcflu_Forecast.Rnw and related files for the article.

        

        

        

 