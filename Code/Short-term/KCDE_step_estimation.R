# KCDE short-term estimation
# Author: Evan Ray, modified by Junyi Lu.
# Date: 16 Aug 2017
# Code comes from 
# https://github.com/reichlab/article-disease-pred-with-kcde/blob/master/inst/code/estimation/kcde-estimation-step.R

# Modification:
# - The loading data and saving file path are changed. 

library(lubridate)
library(plyr)
library(dplyr)
library(reshape)
library(kcde)
library(doMC)
library(here)

## Initialization rule -- "cov" or "scott" for covariance or scott's rule, respectively
bw_init_rule <- "scott"
data_set <- "ili_national"
max_lag <- 1L
max_seasonal_lag <- 0L
filtering <- FALSE
differencing <- FALSE
seasonality <- TRUE
bw_parameterization <- "full"
save_path <- here( "./Results/Forecast_ph1-4/KCDEresults/")

load(file = here("./Data/data_holidays.RData"))

data <- data[data$train == TRUE,]

prediction_target_var <- "weighted_ili"
continuous_var_names <- c(
  paste0(c("weighted_ili", "filtered_weighted_ili"), "_horizon", rep(1:52, each=2)),
  paste0(c("weighted_ili", "filtered_weighted_ili"), "_lag", rep(seq(from = 0, to = max_lag + 52 * max_seasonal_lag), each=2))
)
discrete_var_names <- NULL
predictive_vars <- c("weighted_ili")
time_var <- "time"

kernel_fn <- log_pdtmvn_mode_centered_kernel
rkernel_fn <- rlog_pdtmvn_mode_centered_kernel
initialize_kernel_params_fn <- initialize_params_log_pdtmvn_kernel
get_theta_optim_bounds_fn <- get_theta_optim_bounds_log_pdtmvn_kernel
vectorize_kernel_params_fn <- vectorize_params_log_pdtmvn_kernel
update_theta_from_vectorized_theta_est_fn <- update_theta_from_vectorized_theta_est_log_pdtmvn_kernel

variable_selection_method <- "all_included"
crossval_buffer <- ymd("2014-01-01") - ymd("2013-01-01")



for (prediction_horizon in 1:32) {
  ### Assemble control parameters for KCDE estimation process
  
  ## List describing kernel components -- depends on the values of
  ## prediction_horizon, filtering, seasonality, and bw_parameterization
  kernel_components <- list()
  
  ## sample size for initialize_kernel_params_args
  if(identical(bw_init_rule, "cov")) {
    init_sample_size <- 1L
  } else {
    init_sample_size <- nrow(data)
  }
  
  ## If requested, periodic kernel component capturing seasonality
  if(seasonality) {
    kernel_components <- c(kernel_components,
                           list(list(
                             vars_and_offsets = data.frame(var_name = "time_index",
                                                           offset_value = 0L,
                                                           offset_type = "lag",
                                                           combined_name = "time_index_lag0",
                                                           stringsAsFactors = FALSE),
                             kernel_fn = periodic_kernel,
                             theta_fixed = list(period = pi / 365.2425), # 365.2425 is the mean number of days in a year
                             theta_est = list("bw"),
                             initialize_kernel_params_fn = initialize_params_periodic_kernel,
                             initialize_kernel_params_args = list(
                               sample_size = init_sample_size
                             ),
                             get_theta_optim_bounds_fn = get_theta_optim_bounds_periodic_kernel,
                             get_theta_optim_bounds_args = NULL,
                             vectorize_kernel_params_fn = vectorize_params_periodic_kernel,
                             vectorize_kernel_params_args = NULL,
                             update_theta_from_vectorized_theta_est_fn =
                               update_theta_from_vectorized_theta_est_periodic_kernel,
                             update_theta_from_vectorized_theta_est_args = NULL
                           )))
  }
  
  ## Kernel components for observed values of incidence
  ## First step is setup: create list of data frames specifying groups of
  ## variables and offsets included in each kernel component
  lag_values <- NULL
  for(seasonal_lag in seq(from = 0, to = max_seasonal_lag)) {
    lag_values <- c(lag_values,
                    seq(from = 0, to = max_lag) + 52 * seasonal_lag)
  }
  
  if(identical(bw_parameterization, "diagonal")) {
    ## Separate kernel components for each prediction target variable and
    ## predictive variable
    
    vars_and_offsets_groups <- list()
    
    ## Group of variable names and offsets for prediction target
    new_vars_and_offsets_group <- data.frame(
      var_name = prediction_target_var,
      offset_value = prediction_horizon,
      offset_type = "horizon",
      stringsAsFactors = FALSE
    )
    new_vars_and_offsets_group$combined_name <- paste0(
      new_vars_and_offsets_group$var_name,
      "_",
      new_vars_and_offsets_group$offset_type,
      new_vars_and_offsets_group$offset_value
    )
    vars_and_offsets_groups <- c(vars_and_offsets_groups,
                                 list(new_vars_and_offsets_group))
    
    ## Groups of variable names and offsets for lagged predictive variables
    
    for(lag_value in lag_values) {
      for(predictive_var in predictive_vars) {
        if(filtering) {
          ## If requested, group for lagged filtered observed incidence
          new_vars_and_offsets_group <- data.frame(
            var_name = paste0("filtered_", predictive_var),
            offset_value = lag_value,
            offset_type = "lag",
            stringsAsFactors = FALSE
          )
          new_vars_and_offsets_group$combined_name <- paste0(
            new_vars_and_offsets_group$var_name,
            "_",
            new_vars_and_offsets_group$offset_type,
            new_vars_and_offsets_group$offset_value
          )
          vars_and_offsets_groups <- c(vars_and_offsets_groups,
                                       list(new_vars_and_offsets_group))
        } else {
          ## Else, group for lagged "raw"/unfiltered observed incidence
          new_vars_and_offsets_group <- data.frame(
            var_name = predictive_var,
            offset_value = lag_value,
            offset_type = "lag",
            stringsAsFactors = FALSE
          )
          new_vars_and_offsets_group$combined_name <- paste0(
            new_vars_and_offsets_group$var_name,
            "_",
            new_vars_and_offsets_group$offset_type,
            new_vars_and_offsets_group$offset_value
          )
          vars_and_offsets_groups <- c(vars_and_offsets_groups,
                                       list(new_vars_and_offsets_group))
        }
      }
    }
  } else if(identical(bw_parameterization, "full")) {
    ## One kernel component for prediction target variable and all predictive
    ## variables
    
    ## Prediction target variable
    new_vars_and_offsets_group <- data.frame(
      var_name = prediction_target_var,
      offset_value = prediction_horizon,
      offset_type = "horizon",
      stringsAsFactors = FALSE
    )
    
    ## Lagged prediction target == predictive variables
    for(lag_value in lag_values) {
      for(predictive_var in predictive_vars) {
        if(filtering) {
          ## If requested, lagged filtered incidence
          new_vars_and_offsets_group <- rbind(
            new_vars_and_offsets_group,
            data.frame(
              var_name = paste0("filtered_", predictive_var),
              offset_value = lag_value,
              offset_type = "lag",
              stringsAsFactors = FALSE
            )
          )
        } else {
          ## Else, lagged "raw"/unfiltered observed incidence
          new_vars_and_offsets_group <- rbind(
            new_vars_and_offsets_group,
            data.frame(
              var_name = predictive_var,
              offset_value = lag_value,
              offset_type = "lag",
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
    
    ## Add combined_name column and put in a list for further processing below
    new_vars_and_offsets_group$combined_name <- paste0(
      new_vars_and_offsets_group$var_name,
      "_",
      new_vars_and_offsets_group$offset_type,
      new_vars_and_offsets_group$offset_value
    )
    vars_and_offsets_groups <- list(new_vars_and_offsets_group)
  } else {
    stop("Invalid bandwidth parameterization")
  }
  
  ## Second step is to actually append the kernel component descriptions to the
  ## kernel_components list
  
  if(data_set %in% c("ili_national", "dengue_sj")) {
    #' Compute whether each element of x is equal to the log of an integer,
    #' up to a specified tolerance level.
    #' 
    #' @param x numeric
    #' @param tolerance numeric tolerance for comparison of integer values
    #' 
    #' @return logical vector of same length as x; entry i is TRUE if
    #'     x[i] is within tol of as.integer(x[i])
    equals_log_integer <- function(x, tolerance = .Machine$double.eps ^ 0.5) {
      return(sapply(x, function(x_i) {
        return(isTRUE(all.equal(x_i, log(as.integer(exp(x_i))))))
      }))
    }
    
    #' Compute log(exp(x) - 0.5)
    #' Used as default "a" function
    #' 
    #' @param x numeric
    #' 
    #' @return floor(x) - 1
    log_exp_x_minus_0.5 <- function(x) {
      temp <- exp(x) - 0.5
      temp[temp < 0] <- 0
      return(log(temp))
    }
    
    #' Compute log(exp(x) + 0.5)
    #' Used as default "a" function
    #' 
    #' @param x numeric
    #' 
    #' @return floor(x) - 1
    log_exp_x_plus_0.5 <- function(x) {
      return(log(exp(x) + 0.5))
    }
    
    #' Compute log(round(exp(x))) in such a way that the rounding function always rounds up _.5
    #' 
    #' @param x numeric
    #' 
    #' @return floor(x) - 1
    log_round_up_0.5_exp <- function(x) {
      exp_x <- exp(x)
      
      inds_ceil <- exp_x - floor(exp_x) >= 0.5
      
      exp_x[inds_ceil] <- ceiling(exp_x[inds_ceil])
      exp_x[!inds_ceil] <- floor(exp_x[!inds_ceil])
      
      return(log(exp_x))
    }
    
    #' Compute log(round(exp(x))) in such a way that the rounding function
    #' always rounds up or down to an integer + 0.5, and
    #' an integer always gets rounded up.
    #' 
    #' @param x numeric
    #' 
    #' @return floor(x) - 1
    log_round_to_integer_plus_0.5_exp <- function(x) {
      exp_x <- exp(x) + 0.5
      
      inds_ceil <- exp_x - floor(exp_x) >= 0.5
      
      exp_x[inds_ceil] <- ceiling(exp_x[inds_ceil])
      exp_x[!inds_ceil] <- floor(exp_x[!inds_ceil])
      
      return(log(exp_x - 0.5))
    }
    
    
    a_fn <- log_exp_x_minus_0.5
    b_fn <- log_exp_x_plus_0.5
    #    discretizer_fn <- log_round_up_0.5_exp
    discretizer_fn <- log_round_to_integer_plus_0.5_exp
    in_range_fn <- function(x, tolerance = 0.5 * .Machine$double.eps^0.5) {
      return(sapply(x, function(x_i) {
        return(
          isTRUE(all.equal(
            x_i,
            #                    log_round_up_0.5_exp(x_i),
            log_round_to_integer_plus_0.5_exp(x_i),
            tolerance = tolerance
          ))
        )
      }))
    }
  } else {
    a_fn <- x_minus_0.25
    b_fn <- x_plus_0.25
    in_range_fn <- equals_half_integer
    discretizer_fn <- round_to_half_integer
  }
  
  kernel_components <- c(kernel_components,
                         lapply(vars_and_offsets_groups, function(vars_and_offsets) {
                           lower_trunc_bds <- rep(-Inf, nrow(vars_and_offsets))
                           names(lower_trunc_bds) <- vars_and_offsets$combined_name
                           upper_trunc_bds <- rep(Inf, nrow(vars_and_offsets))
                           names(upper_trunc_bds) <- vars_and_offsets$combined_name
                           
                           
                           discrete_var_range_fns <- NULL
                           
                           
                           return(list(
                             vars_and_offsets = vars_and_offsets,
                             kernel_fn = kernel_fn,
                             rkernel_fn = rkernel_fn,
                             theta_fixed = list(
                               parameterization = "bw-chol-decomp",
                               continuous_vars = vars_and_offsets$combined_name[
                                 vars_and_offsets$combined_name %in% continuous_var_names],
                               discrete_vars = vars_and_offsets$combined_name[
                                 vars_and_offsets$combined_name %in% discrete_var_names],
                               discrete_var_range_fns = discrete_var_range_fns,
                               lower = lower_trunc_bds,
                               upper = upper_trunc_bds,
                               validate_in_support = FALSE
                             ),
                             theta_est = list("bw"),
                             initialize_kernel_params_fn = initialize_kernel_params_fn,
                             initialize_kernel_params_args = list(
                               sample_size = init_sample_size
                             ),
                             get_theta_optim_bounds_fn = get_theta_optim_bounds_fn,
                             get_theta_optim_bounds_args = NULL,
                             vectorize_kernel_params_fn = vectorize_kernel_params_fn,
                             vectorize_kernel_params_args = NULL,
                             update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_fn,
                             update_theta_from_vectorized_theta_est_args = NULL
                           ))
                         })
  )
  
  
  ## Set up filter_control to do no filtering
  filter_control <- NULL
  
  ## Assemble kernel_components and filter_control created above,
  ## along with other parameters controling KCDE definition and estimation
  kcde_control <- create_kcde_control(X_names = "time_index",
                                      y_names = prediction_target_var,
                                      time_name = time_var,
                                      prediction_horizons = prediction_horizon,
                                      filter_control = filter_control,
                                      kernel_components = kernel_components,
                                      crossval_buffer = crossval_buffer,
                                      loss_fn = neg_log_score_loss,
                                      loss_fn_prediction_args = list(
                                        prediction_type = "distribution",
                                        log = TRUE),
                                      loss_args = NULL,
                                      par_cores = 4L,
                                      variable_selection_method = variable_selection_method,
                                      prediction_inds_not_included = NULL)
  
  # add prediction_inds_not_included = NULL/11.01.2019
  
  ### Do estimation
  ## Read in output from an earlier run if it exists.
  ## We started some runs that cut off by the cluster either because of run time
  ## limits or some sort of cluster I/O issue.
  ## Read in output to get initial values for parameters that are the
  ## best that were realized in that earlier run.
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
  
 
  init_theta_vector <- NULL
  init_phi_vector <- NULL
  library("doMC")
  registerDoMC(cores = kcde_control$par_cores)
  
  ## Get the KCDE fit
  fit_time <- system.time({
    kcde_fit <- kcde(data = data,
                     kcde_control = kcde_control,
                     init_theta_vector = init_theta_vector,
                     init_phi_vector = init_phi_vector)
  })
  
  
  ### Save results
  saveRDS(kcde_fit,
          file = file.path(save_path,
                           paste0("kcde_fit-",
                                  case_descriptor,
                                  ".rds")
          )
  )
  
  saveRDS(fit_time,
          file = file.path(save_path,
                           paste0("fit_time-",
                                  case_descriptor,
                                  ".rds")
          )
  )
}
