#' Resample using permutations
#'
#' This function reruns the desired wintime package method on a given number of permutations. This resampling method is recommended
#' for the Expected win time (EWT) and Restricted mean survival in favor of treatment (RMT) methods.
#'
#' @param type A string value indicating the wintime package method that will run with resampling.
#' @param rmst_restriction The RMT cutoff value (days).
#' @param model A string value indicating the model used on observed data ('markov' or 'km').
#' @param n The total number of trial participants.
#' @param m The number of events in the hierarchy.
#' @param Time A m x n matrix of event times (days). Rows should represent events and columns should represent participants. Event rows should be
#' in increasing order of clinical severity.
#' @param Delta A m x n matrix of event indicators. Rows should represent events and columns should represent participants. Event rows should
#' be in increasing order of clinical severity.
#' @param trt A numeric vector of treatment arm indicators (1 for treatment, 0 for control).
#' @param cov A n x p matrix of covariate values, where p is the number of covariates. Rows should represent participants and columns should
#' represent covariate values.
#' @param z_ewtr The Z-statistic of EWTR.
#' @param z_comp The Z-statistic of the composite event approach.
#' @param resample_num The number of desired permutations.
#' @param seed The seed used for random number generation.
#' @return A vector of length resample_num containing the treatment effect estimates (for type='max' these are z-statistics) for each permutation.

# -------------------------------
# Permutations
# -------------------------------
perm <- function(type,rmst_restriction,model,n,m,Time,Delta,trt,cov,z_ewtr,z_comp,resample_num,seed) {
  # Initialize vector to hold permuted data values
  y <- rep(0,times = resample_num)
  time <- Time[m:1, ]
  delta <- Delta[m:1, ]
  iperm <- 1
  while (iperm <= resample_num) {
    # Generate permutations
    perm_indices <- sample(1:n, size = length(trt), replace = FALSE)
    trt_perm <- trt[perm_indices]
    Time_perm <- Time[, perm_indices, drop = FALSE]
    Delta_perm <- Delta[, perm_indices, drop = FALSE]
    time_perm <- time[, perm_indices, drop = FALSE]
    delta_perm <- delta[, perm_indices, drop = FALSE]
    cov_perm <- cov[perm_indices,, drop = FALSE]

    # Set parameters for model fitting calls
    n0 <- sum(trt_perm == 0)
    n1 <- sum(trt_perm == 1)
    trt_perm[1:n0] <- 0
    trt_perm[(n0+1):n] <- 1
    m <- nrow(Time_perm)
    tau <- max(Time_perm[m, ])

    # Type function calls on permuted data
    z <- NULL
    markov_ind <- FALSE
    if (!is.null(model)) {
      # If a Markov model is specified, call the function to fit a Markov model
      if (model == "markov") {
        z <- markov(n0,n1,m,Time_perm,Delta_perm)
        markov_ind <- TRUE
      }
      # If no model is specified, use the default (KM) model
      else {
        z <- km(n0,n1,m,Time_perm,Delta_perm)
      }
    }

    if (!is.null(z)) {
      # Parameters for type function calls
      dist_state0_perm <- z[[1]]
      dist_state1_perm <- z[[2]]
      untimes0_perm <- z[[3]]
      untimes1_perm <- z[[4]]
      nuntimes0_perm <- z[[5]]
      nuntimes1_perm <- z[[6]]
      max_follow0_perm <- z[[7]]
      max_follow1_perm <- z[[8]]
    }

    # Type function calls
    type <- tolower(type)
    if (type == "ewt") {
      y[iperm] <- EWT(m,dist_state0_perm,dist_state1_perm,untimes0_perm,untimes1_perm,nuntimes0_perm,nuntimes1_perm)
    }
    else if (type == "ewtr") {
      y[iperm] <- EWTR(n,m,nuntimes0_perm,max_follow0_perm,untimes0_perm,Time_perm,Delta_perm,dist_state0_perm,markov_ind,cov_perm,trt_perm)[[1]]
    }
    else if (type == "rmt") {
      y[iperm] <- RMT(m,rmst_restriction,dist_state0_perm,dist_state1_perm,untimes0_perm,untimes1_perm,nuntimes0_perm,nuntimes1_perm)
    }
    else if (type == "wtr") {
      y[iperm] <- WTR(n,m,tau,trt_perm,Time_perm,Delta_perm)[[1]]
    }
    else if (type == "rwtr") {
      y[iperm] <- RWTR(n,m,tau,trt_perm,Time_perm,Delta_perm)[[1]]
    }
    else if (type == "pwt") {
      y[iperm] <- PWT(n,n0,n1,Time_perm,Delta_perm,trt_perm,m,tau)
    }
    else if (type == "max") {
      # Permutation Z statistic
      z_ewtr_perm <- EWTR(n,m,nuntimes0_perm,max_follow0_perm,untimes0_perm,Time_perm,Delta_perm,dist_state0_perm,markov_ind,cov_perm,trt_perm)[[3]]
      z_comp_perm <- COMP(n,Time_perm,Delta_perm,cov_perm,trt_perm)[[1]]
      y[iperm] <- max((z_ewtr_perm - z_ewtr),(z_comp_perm - z_comp))
    }
    else {
      stop(paste("Invalid type:", type, "- Please specify one of 'ewt', 'ewtr', 'rmt', 'wtr', 'rwtr', 'pwt', 'max'"))
    }
    iperm <- iperm + 1
  }
  return(y)
}
