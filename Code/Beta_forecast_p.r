library(betareg)

logit_FUN <- function(x){
  qlogis(x/100)
}

# simulation of beta model with higher lags.
Beta_Forecast_p <- function(object,
                            nsim = 1,
                            seed = 1,
                            ph = 1,
                            p = 1,
                            start_value,
                            coefs_subset){
  set.seed(seed)
  
  pt <- shape1 <- shape2 <- mean <- precision <-
    as.data.frame(matrix(NA_real_, nrow = nsim, ncol = ph,
                         dimnames = list(NULL, paste("ph", 1 : ph))))
  
  starti <- as.data.frame(matrix(sapply(start_value, rep, each = nsim), nrow = nsim))
  for(i in 1 : ph){
    vari <- sapply(starti, logit_FUN)
    coefi <- data.frame(vari,
                        coefs_subset[i, ])
    names(coefi)[1 : p] <- paste0("p", p : 1)
    mean[, i] <- meani <- predict(object,
                                  coefi,
                                  type = "response")
    precision[, i] <- precisioni <- predict(object,
                                            coefi,
                                            type = "precision")
    shape1[, i] <- shape1i <- meani * precisioni
    shape2[, i] <- shape2i <- precisioni - shape1i
    if(p > 1){
      starti[, 1 : (p - 1)] <- starti[, 2 : p]
      starti[, p] <-  100 * rbeta(n = nsim,
                                  shape1 = shape1i,
                                  shape2 = shape2i)
      pt[, i] <- starti[, p]
      
    }else{
      starti <- 100 * rbeta(n = nsim,
                            shape1 = shape1i,
                            shape2 = shape2i)
      pt[, i] <- starti
    }
  }
  
  
  res <- list(pt = pt,
              mean = mean,
              precision = precision,
              shape1 = shape1,
              shape2 = shape2)
  return(res)
}

# obs is percentage
cal_log_score <- function(sim_object,
                          obs,
                          ph){
  shape1 <- sim_object$shape1[, ph]
  shape2 <- sim_object$shape2[, ph]
  dens <- mapply(dbeta,
                 x = obs/100,
                 shape1 = shape1,
                 shape2 = shape2)
  log(mean(dens)) - log(100)
}

