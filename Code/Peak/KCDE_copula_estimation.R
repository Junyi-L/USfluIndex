library(lubridate)
library(plyr)
library(dplyr)
library(reshape)
library(kcde)
library(copula)
library(mvtnorm)
library(here)
## There are some problems with the implementation of copula estimation in the
## copula package, mainly centering around numerical issues and handling
## negative values.  I've made some minor changes to fix these problems.
## I haven't made these a part of the copula package, but have just put the
## revised functions here.

fitCopula.ml <- function (copula, u, start, lower, upper, method, optim.control, 
                          estimate.variance, hideWarnings, bound.eps = .Machine$double.eps^0.5)
{
  if (any(u < 0) || any(u > 1)) 
    stop("'u' must be in [0,1] -- probably rather use pobs(.)")
  stopifnot(is.numeric(d <- ncol(u)), d >= 2)
  if (copula@dimension != d) 
    stop("The dimension of the data and copula do not match")
  if (is.null(start)) 
    start <- fitCopStart(copula, u)
  if (any(is.na(start))) 
    stop("'start' contains NA values")
  
  ## NEW -- make start give a p.d. matrix.
  copula@parameters <- start
  while(min(eigen(getSigma(copula))$values) < 10^{-6}) {
    sigma <- getSigma(copula)
    sigma <- makePosDef(sigma)
    start <- sigma[1, seq(from = 2, to = ncol(sigma))]
    copula@parameters <- start
  }
  
  q <- length(copula@parameters)
  if (q != length(start)) 
    stop(gettextf("The lengths of 'start' (= %d) and copula@parameters (=%d) differ", 
                  length(start), q), domain = NA)
  control <- c(optim.control, fnscale = -1)
  control <- control[!vapply(control, is.null, NA)]
  if (!is.null(optim.control[[1]])) 
    control <- c(control, optim.control)
  meth.has.bounds <- method %in% c("Brent", "L-BFGS-B")
  if (is.null(lower)) 
    lower <- if (meth.has.bounds) 
      copula@param.lowbnd + bound.eps
  else -Inf
  if (is.null(upper)) 
    upper <- if (meth.has.bounds) 
      copula@param.upbnd - bound.eps
  else Inf
  (if (hideWarnings) 
    suppressWarnings
    else identity)(fit <- optim(start, loglikCopula, lower = lower, 
                                upper = upper, method = method, copula = copula, x = u, 
                                control = control))
  copula@parameters[1:q] <- fit$par
  loglik <- fit$val
  has.conv <- fit[["convergence"]] == 0
  if (is.na(estimate.variance)) 
    estimate.variance <- has.conv
  if (!has.conv) 
    warning("possible convergence problem: optim() gave code=", 
            fit$convergence)
  varNA <- matrix(NA_real_, q, q)
  var.est <- if (estimate.variance) {
    fit.last <- optim(copula@parameters, loglikCopula, lower = lower, 
                      upper = upper, method = method, copula = copula, 
                      x = u, control = c(control, maxit = 0), hessian = TRUE)
    vcov <- tryCatch(solve(-fit.last$hessian), error = function(e) e)
    if (is(vcov, "error")) {
      warning("Hessian matrix not invertible: ", vcov$message)
      varNA
    }
    else vcov
  }
  else varNA
  new("fitCopula", estimate = fit$par, var.est = var.est, method = "maximum likelihood", 
      loglik = loglik, fitting.stats = c(list(method = method), 
                                         fit[c("convergence", "counts", "message")], control), 
      nsample = nrow(u), copula = copula)
}

fitCopula <- function (copula, data, method = c("mpl", "ml", "itau", "irho"), 
                       start = NULL, lower = NULL, upper = NULL, optim.method = "BFGS", 
                       optim.control = list(maxit = 1000), estimate.variance = NA, 
                       hideWarnings = TRUE) 
{
  if (!is.matrix(data)) {
    warning("coercing 'data' to a matrix.")
    data <- as.matrix(data)
    stopifnot(is.matrix(data))
  }
  switch(match.arg(method), ml = fitCopula.ml(copula, data, 
                                              start = start, lower = lower, upper = upper, method = optim.method, 
                                              optim.control = optim.control, estimate.variance = estimate.variance, 
                                              hideWarnings = hideWarnings), mpl = fitCopula.mpl(copula, 
                                                                                                data, start = start, lower = lower, upper = upper, optim.method = optim.method, 
                                                                                                optim.control = optim.control, estimate.variance = estimate.variance, 
                                                                                                hideWarnings = hideWarnings), itau = fitCopula.itau(copula, 
                                                                                                                                                    data, estimate.variance = estimate.variance), irho = fitCopula.irho(copula, 
                                                                                                                                                                                                                        data, estimate.variance = estimate.variance))
}

loglikCopula <- function (param, x, copula, hideWarnings) 
{
  stopifnot(length(copula@parameters) == length(param))
  if (!missing(hideWarnings)) 
    warning("'hideWarnings' is deprecated and has no effect anymore")
  copula@parameters <- param
  if (chkParamBounds(copula)) { 
    #        cat(sum(dCopula(x, copula, log = TRUE, checkPar = FALSE)))
    #        cat("\n")
    retval <- max(-10^6, sum(dCopula(x, copula, log = TRUE, checkPar = FALSE)))
  }
  else retval <- -10^6
  #    cat(retval)
  #    cat("\n")
  return(retval)
}

fitCopStart <- function (copula, data, default = copula@parameters) 
{
  if (hasMethod("iTau", clc <- class(copula))) {
    ccl <- getClass(clc)
    .par.df <- has.par.df(copula, ccl)
    start <- fitCopula.itau(if (.par.df) 
      as.df.fixed(copula, ccl)
      else copula, data, estimate.variance = FALSE, warn.df = FALSE)@estimate
    if (.par.df) 
      start <- c(start, copula@df)
    if (!is.finite(loglikCopula(start, data, copula))) 
      default
    else start
  }
  else default
}

fitCopula.itau <- function (copula, x, estimate.variance, warn.df = TRUE) 
{
  ccl <- getClass(class(copula))
  isEll <- extends(ccl, "ellipCopula")
  if (has.par.df(copula, ccl, isEll)) {
    if (warn.df) 
      warning("\"itau\" fitting ==> copula coerced to 'df.fixed=TRUE'")
    copula <- as.df.fixed(copula, classDef = ccl)
  }
  q <- length(copula@parameters)
  tau <- cor(x, method = "kendall")
  itau <- fitKendall(copula, tau)
  itau <- P2p(itau)
  if(ncol(x) > 2) {
    X <- getXmat(copula)
    estimate <- as.vector(if (isEll && copula@dispstr == "ar1") 
      exp(lm.fit(X, y = log(itau))$coefficients)
      else lm.fit(X, y = itau)$coefficients)
  } else {
    estimate <- itau
  }
  copula@parameters <- estimate
  var.est <- if (is.na(estimate.variance) || estimate.variance) 
    varKendall(copula, x)/nrow(x)
  else matrix(NA, q, q)
  new("fitCopula", estimate = estimate, var.est = var.est, 
      method = "inversion of Kendall's tau", loglik = NA_real_, 
      fitting.stats = list(convergence = NA_integer_), nsample = nrow(x), 
      copula = copula)
}

has.par.df <- function (cop, classDef = getClass(class(cop)), isEllip = extends(classDef, 
                                                                                "ellipCopula")) 
{
  ((isEllip && extends(classDef, "tCopula")) || extends(classDef, 
                                                        "tevCopula")) && !cop@df.fixed
}

fitKendall <- function (cop, tau) 
{
  stopifnot(is.numeric(p <- ncol(tau)), p == nrow(tau))
  sigma <- matrix(1, p, p)
  for (j in 1:(p - 1)) for (i in (j + 1):p) {
    sigma[i, j] <- iTau(cop, tau[i, j])
    sigma[j, i] <- sigma[i, j]
  }
  if (is(cop, "ellipCopula")) 
    makePosDef(sigma, delta = 0.001)
  else sigma
}

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

getXmat <- function (copula) 
{
  p <- copula@dimension
  pp <- p * (p - 1)/2
  if (!is(copula, "ellipCopula")) 
    matrix(1, nrow = pp, ncol = 1)
  else {
    switch(copula@dispstr, ex = matrix(1, nrow = pp, ncol = 1), 
           un = diag(pp), toep = , ar1 = {
             dgidx <- outer(1:p, 1:p, "-")
             dgidx <- P2p(dgidx)
             if (copula@dispstr == "toep") model.matrix(~factor(dgidx) - 
                                                          1) else {
                                                            cbind(dgidx, deparse.level = 0L)
                                                          }
           }, stop("Not implemented yet for the dispersion structure ", 
                   copula@dispstr))
  }
}

P2p <- function (P) {
  P[lower.tri(P)]
}

chkParamBounds <- function (copula) 
{
  d <- length(param <- copula@parameters)
  m <- length(upper <- copula@param.upbnd)
  if (d != m) 
    return(FALSE)
  m <- length(lower <- copula@param.lowbnd)
  if (d != m) 
    return(FALSE)
  !(any(is.na(param) | param > upper | param < lower))
}

dnormalCopula <- function(u, copula, log=FALSE, ...) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  #    sigma <- makePosDef(sigma)
  if(!is.matrix(u)) u <- matrix(u, ncol = dim)
  r <- numeric(nrow(u)) # i.e. 0  by default (i.e. "outside")
  ok <- !apply(u, 1, function(x) any(is.na(x)))
  x <- qnorm(u[ok, , drop=FALSE])
  ## work in log-scale [less over-/under-flow, then (maybe) transform:
  r[ok] <- dmvnorm(x, sigma = sigma, log=TRUE) - rowSums(dnorm(x, log=TRUE))
  ## now happens in dCopula(): -- dnormalCopula() not called directly by user
  ## if(any(out <- !is.na(u) & (u <= 0 | u >= 1)))
  ##   val[apply(out, 1, any)] <- -Inf
  if(log) r else exp(r)
}

setMethod("dCopula", signature("matrix", "normalCopula"), dnormalCopula)
setMethod("dCopula", signature("numeric", "normalCopula"),dnormalCopula)

data_set <- "ili_national"


save_path <- file.path(here("./Results/Peak/copula-estimation-results"))
estimation_results_path <- file.path(here("./Results/Forecast_ph1-4/KCDEresults"))


all_max_lags <- as.character(c(1L)) # use incidence at times t^* and t^* - 1 to predict incidence after t^*
all_max_seasonal_lags <- as.character(0L) # not used
all_filtering_values <- c("FALSE") # not used
all_differencing_values <- "FALSE" # not used
all_seasonality_values <- "TRUE" # specifications without and with periodic kernel
all_bw_parameterizations <- "full" # specifications with diagonal and full bandwidth
all_sim_n <- "NA" # not used for applications
all_sim_families <- "NA" # not used for applications
all_sim_run_inds <- 1L # not used for applications

for(sim_family in all_sim_families) {
  for(sim_n in all_sim_n) {
    for(sim_run_ind in all_sim_run_inds) {
      for(max_lag in all_max_lags) {
        for(max_seasonal_lag in all_max_seasonal_lags) {
          for(filtering in all_filtering_values) {
            for(differencing in all_differencing_values) {
              for(seasonality in all_seasonality_values) {
                for(bw_parameterization in all_bw_parameterizations) {
                  
                  ptm <- proc.time()
                  copula_fits_case_descriptor <- paste0(
                    data_set,
                    "-max_lag_", max_lag,
                    "-max_seasonal_lag_", max_seasonal_lag,
                    "-filtering_", filtering,
                    "-differencing_", differencing,
                    "-seasonality_", seasonality,
                    "-bw_parameterization_", bw_parameterization
                  )
                  copula_fits_file_name <- paste0("kcde-copula-fits-",
                                                  copula_fits_case_descriptor,
                                                  ".rds")
                  
                  if(identical(data_set, "ili_national")) {
                    load(file = here("./Data/data_holidays.RData"))
                    data <- data[data$train == TRUE, ]
 
                    prediction_target_var <- "weighted_ili"
                    train_seasons <- paste0(seq(from = 1998, to = 2013), "/", seq(from = 1999, to = 2014))
                    
                    season_length <- 33L
                    first_analysis_time_season_week <- 10 # == week 40 of year
                    last_analysis_time_season_week <- 41 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
                  } 
                  
                  copula_fits <- vector("list", last_analysis_time_season_week - first_analysis_time_season_week + 1)
                  copula_fits_ind <- 1L
                  
                  cat(paste0("seasonality = ", seasonality, "\n"))
                  cat(paste0("bw_parameterization = ", bw_parameterization, "\n"))
                  
                  ## We obtain a separate copula estimate for each combination of
                  ## data set and number of weeks left in season
                  ## omit the last week since there is only a prediction horizon of 1 then.
                  for(analysis_time_season_week in seq(from = first_analysis_time_season_week, to = last_analysis_time_season_week - 1)) {
                    cat(paste0("analysis_time_season_week = ", analysis_time_season_week, "\n"))
                    
                    ## Assemble matrix of probability-integral-transformed observed incidence trajectories in each season
                    ## The estimated KCDE density is used for the integral transform.
                    ## Rows correspond to seasons, columns to prediction horizons
                    max_prediction_horizon <-
                      first_analysis_time_season_week + season_length - 1 -
                      analysis_time_season_week
                    
                    pit_incidence_trajectories <- t(sapply(train_seasons,
                                                           function(train_season) {
                                                             ## Calculate conditional probability that incidence <= observed
                                                             ## incidence for each prediction horizon separately
                                                             analysis_time_ind <- which(data$season == train_season & data$season_week == analysis_time_season_week)
                                                             
                                                             sapply(seq_len(max_prediction_horizon),
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
                                                                      
                                                                      kcde_fit_file_path <- file.path(estimation_results_path,
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
                                                                      
                                                                      ## Get the probability integral transforms of the observed incidence trajectories.
                                                                      tryCatch(
                                                                        kcde_predict(
                                                                          q = data[analysis_time_ind + prediction_horizon, prediction_target_var, drop = FALSE],
                                                                          n = 10000,
                                                                          kcde_fit = kcde_fit,
                                                                          prediction_data =
                                                                            data[seq_len(analysis_time_ind), , drop = FALSE],
                                                                          leading_rows_to_drop = 0L,
                                                                          trailing_rows_to_drop = 0L,
                                                                          additional_training_rows_to_drop = NULL,
                                                                          prediction_type = "prob",
                                                                          log = TRUE
                                                                        ),
                                                                        error = function(e) {
                                                                          if(identical(e$message, "0 rows in cross validation data after filtering, computing offsets and dropping NA rows")) {
                                                                            ## Data at analysis time did not include sufficient lags without NA values to use
                                                                            ## Returning NA allows us to continue without including this analysis time
                                                                            return(NA)
                                                                          } else {
                                                                            stop(e$message)
                                                                          }
                                                                        }
                                                                      )
                                                                    }
                                                             )
                                                           }
                    ))
                    
                    ## Drop NA rows
                    ## These may occur in early seasons if lag > 0 was used
                    na_rows <- which(apply(pit_incidence_trajectories, 1, function(pit_row) {any(is.na(pit_row))}))
                    if(length(na_rows) > 0) {
                      pit_incidence_trajectories <- pit_incidence_trajectories[-na_rows, ]
                    }
                    ## The evaluation of predictive distribution above is not exact, and may
                    ## give values of 0 or 1.  We require values between 0 and 1 for copula
                    pit_incidence_trajectories[pit_incidence_trajectories == 0] <- 10^-6
                    pit_incidence_trajectories[pit_incidence_trajectories == 1] <- 1 - 10^-6
                    
                    ## Fit copula
                    copula_fit <- fitCopula(
                      copula = normalCopula(param = rep(0.3, max_prediction_horizon - 1),
                                            dim = max_prediction_horizon,
                                            dispstr = "toep"),
                      data = pit_incidence_trajectories,
                      estimate.variance = TRUE,
                      method = "ml",
                      optim.method = "L-BFGS-B",
                      start = NULL)
                    
                    copula_fits[[copula_fits_ind]] <- list(
                      analysis_time_season_week = analysis_time_season_week,
                      pit_incidence_trajectories = pit_incidence_trajectories,
                      copula_fit = copula_fit
                    )
                    
                    copula_fits_ind <- copula_fits_ind + 1L
                  }
                  run_time <- proc.time() - ptm
                  saveRDS(copula_fits,
                          file = file.path(save_path, copula_fits_file_name))
                  saveRDS(run_time,
                          file = file.path(save_path, "kcde-copula-fits-others.rds"))
                } # bw_parameterization
              } # seasonality
            } # differencing
          } # filtering
        } # max_seasonal_lag
      } # max_lag
    } # sim_run_ind
  } # sim_n
} # sim_family


