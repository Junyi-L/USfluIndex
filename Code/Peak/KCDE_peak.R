library(lubridate)
library(plyr)
library(dplyr)
library(reshape)
library(kcde)
library(copula)
library(mvtnorm)
library(doMC)
library(here)
set.seed(1)
## Function borrowed from copula package and modified
makePosDef <- function (mat, delta = 0.001) 
{
  while(min(eigen(mat)$values) < 10^{-6}) {
    decomp <- eigen(mat)
    Lambda <- decomp$values
    Lambda[Lambda < 0] <- delta
    Gamma <- decomp$vectors
    newmat <- Gamma %*% diag(Lambda) %*% t(Gamma)
    D <- 1/sqrt(diag(newmat))
    mat <- diag(D) %*% newmat %*% diag(D)
  }
  return(mat)
}


## Run over cases defined by combinations of
## data set, max lag, max seasonal lag, filtering, differencing, seasonality, bw_parameterization


num_cores <- 4L
registerDoMC(cores = num_cores)


data_set <- "ili_national"
all_prediction_horizons <- as.character(seq_len(52))
all_max_lags <- as.character(c(1L)) # use incidence at times t^* and t^* - 1 to predict incidence after t^*
all_max_seasonal_lags <- as.character(0L) # not used
all_filtering_values <- c("FALSE") # not used
all_differencing_values <- "FALSE" # not used
all_seasonality_values <- "TRUE" # specifications without and with periodic kernel
all_bw_parameterizations <-  "full" # specifications with diagonal and full bandwidth
all_sim_n <- "NA" # not used for applications
all_sim_families <- "NA" # not used for applications
all_sim_run_inds <- 1L # not used for applications

incidence_bins <- data.frame(
  lower = seq(from = 0, to = 13, by = 0.5),
  upper = c(seq(from = 0.5, to = 13, by = 0.5), Inf))



case_definitions <-
  expand.grid(
    data_set,
    all_max_lags,
    all_max_seasonal_lags,
    all_filtering_values,
    all_differencing_values,
    all_seasonality_values,
    all_bw_parameterizations,
    stringsAsFactors = FALSE) %>%
  `colnames<-`(c("data_set",
                 "max_lag",
                 "max_seasonal_lag",
                 "filtering",
                 "differencing",
                 "seasonality",
                 "bw_parameterization"))


junk <- foreach(case_row_ind = seq_len(nrow(case_definitions)),
                .packages = c("kcde", "plyr", "dplyr", "lubridate", "reshape", "copula", "mvtnorm"),
                .combine = "rbind") %dopar% {
                  ptm <- proc.time()
                  data_set <- case_definitions$data_set[case_row_ind]
                  max_lag <- case_definitions$max_lag[case_row_ind]
                  max_seasonal_lag <- case_definitions$max_seasonal_lag[case_row_ind]
                  filtering <- case_definitions$filtering[case_row_ind]
                  differencing <- case_definitions$differencing[case_row_ind]
                  seasonality <- case_definitions$seasonality[case_row_ind]
                  bw_parameterization <- case_definitions$bw_parameterization[case_row_ind]
                  
                  n_sims <- 10000
                  
                  
                  copula_save_path <- file.path(here("./Results/Peak/copula-estimation-results"))
                  estimation_save_path <- file.path(here("./Results/Forecast_ph1-4/KCDEresults"))
                  prediction_save_path <- file.path(here("./Results/Peak"))
                  
                  case_descriptor <- paste0(
                    data_set,
                    "-max_lag_", max_lag,
                    "-max_seasonal_lag_", max_seasonal_lag,
                    "-filtering_", filtering,
                    "-differencing_", differencing,
                    "-seasonality_", seasonality,
                    "-bw_parameterization_", bw_parameterization
                  )
                  file_name <- paste0("kcde-copula-fits-",
                                      case_descriptor,
                                      ".rds")
                  copula_fits <- readRDS(file = file.path(copula_save_path, file_name))
                  analysis_time_season_week_by_copula_fit <- unlist(lapply(copula_fits,
                                                                           function(copula_fit) {copula_fit$analysis_time_season_week}))
                  
                  
                  kcde_fits_by_prediction_horizon <- lapply(seq_len(32),
                                                            function(prediction_horizon) {
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
                                                              
                                                              kcde_fit_file_path <- file.path(estimation_save_path,
                                                                                              paste0("kcde_fit-", case_descriptor, ".rds"))
                                                              kcde_fit <- readRDS(kcde_fit_file_path)
                                                              
                                                              ## fix rkernel_fn for pdtmvn-based kernel functions
                                                              ## I supplied a buggy version of this in the call to the estimation routine
                                                              ## that did not ensure that variables were supplied in a consistent order.
                                                              ## This did not affect estimation as rkernel_fn is not called there
                                                              ## But it does need to be fixed here.
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
                                                              
                                                              return(kcde_fit)
                                                            })
                  
                  
                  
                  load(file = here("./Data/data_holidays.RData"))
                  
                  orig_prediction_target_var <- "weighted_ili"
                  prediction_target_var <- "weighted_ili"
                  
                  analysis_seasons <- c("2014/2015", "2015/2016", "2016/2017", "2017/2018")
                  first_analysis_time_season_week <- 10 # == week 40 of year
                  last_analysis_time_season_week <- 41 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
                  
                  
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
                  
                  ## generate peak week timing and height estimates
                  for(analysis_time_season in analysis_seasons) {
                    ## get observed quantities, for computing log score
                    observed_peak_height <- max(data[data$season == analysis_time_season, orig_prediction_target_var])
                    observed_peak_week_ind <- which((data$season == analysis_time_season) &
                                                      (data[, orig_prediction_target_var] == observed_peak_height))
                    observed_peak_week <- data$season_week[observed_peak_week_ind]
                    
                    observed_peak_height <- which(
                      incidence_bins$lower <= observed_peak_height &
                        incidence_bins$upper > observed_peak_height)
                    
                    for(analysis_time_season_week in seq(from = first_analysis_time_season_week, to = last_analysis_time_season_week - 1)) {
                      ### simulate from copula that ties the marginal predictive distributions together
                      
                      ## get the right copula for analysis_time_season_week
                      predictive_copula_ind <- which(analysis_time_season_week_by_copula_fit == analysis_time_season_week)
                      copula_fit <- copula_fits[[predictive_copula_ind]]$copula_fit
                      predictive_copula <- copula_fit@copula
                      
                      ## simulate n_sims sequences from copula
                      max_prediction_horizon <-
                        last_analysis_time_season_week + 1 -
                        analysis_time_season_week
                      sim_sequences <- rCopula(n_sims, predictive_copula)
                      
                      ## get quantiles from marginal predictive distributions corresponding to
                      ## values simulated from copula
                      analysis_time_ind <- which(data$season == analysis_time_season &
                                                   data$season_week == analysis_time_season_week)
                      trajectory_samples <- matrix(NA, nrow = n_sims, ncol = max_prediction_horizon)
                      for(prediction_horizon in seq_len(max_prediction_horizon)) {
                        trajectory_samples[, prediction_horizon] <-
                          kcde_predict(
                            p = sim_sequences[, prediction_horizon],
                            n = 100000, 
                            kcde_fit = kcde_fits_by_prediction_horizon[[prediction_horizon]],
                            prediction_data =
                              data[seq_len(analysis_time_ind), , drop = FALSE],
                            leading_rows_to_drop = 0L,
                            trailing_rows_to_drop = 0L,
                            additional_training_rows_to_drop = NULL,
                            prediction_type = "quantile",
                            log = TRUE
                          )
                      }
                      
                      ## Augment trajectory samples with previous observed incidence values
                      season_start_ind <- which(data$season == analysis_time_season &
                                                  data$season_week == 1)
                      if(season_start_ind < analysis_time_ind) {
                        trajectory_samples <- cbind(
                          matrix(
                            rep(data[seq(from = season_start_ind, to = analysis_time_ind), orig_prediction_target_var], each = n_sims),
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
                                                              which(incidence_bins$lower <= height &
                                                                      incidence_bins$upper > height)
                                                            })
                      results[results_save_row, paste0("peak_height_", seq_len(n_sims))] <-
                        peak_week_height_by_sim_ind
                      
                      ## Get log scores
                      results[results_save_row, "peak_week_log_score"] <- log(sum(peak_week_by_sim_ind == observed_peak_week)) - log(n_sims)
                      results[results_save_row, "peak_height_log_score"] <- log(sum(peak_week_height_by_sim_ind == observed_peak_height)) - log(n_sims)
                    }
                  }
                  run_time <- proc.time() - ptm
                  case_descriptor <- paste0(
                    data_set,
                    "-max_lag_", max_lag,
                    "-max_seasonal_lag_", max_seasonal_lag,
                    "-filtering_", filtering,
                    "-differencing_", differencing,
                    "-seasonality_", seasonality,
                    "-bw_parameterization_", bw_parameterization
                  )
                  saveRDS(results,
                          file.path(prediction_save_path,
                                    paste0("peak-week-KCDE.rds")))
                  
                  saveRDS(run_time,
                          file.path(prediction_save_path,
                                    paste0("peak-week-KCDE-others.rds")))
                }

