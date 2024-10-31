#' Resample using bootstraps
#'
#' This function reruns the desired wintime package method on a given number of bootstrap samples. This resampling method is recommended
#' for all pairwise wintime methods including Win time ratio (WTR), Restricted win time ratio (RWTR), and Pairwise win time (PWT). This
#' function is also recommended for the EWTR_composite max test (MAX).
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
#' @param resample_num The number of desired bootstraps.
#' @param seed The seed used for random number generation.
#' @return A vector of length resample_num containing the calculated treatment effect estimates (for type='max' these are z-statistics) for each bootstrap.

# ---------------------------
# Bootstrap
# ---------------------------
bootstrap <- function(type,rmst_restriction,model,n,m,Time,Delta,trt,cov,z_ewtr,z_comp,resample_num,seed) {
  # Initialize vector to hold bootstrap values
  y <- rep(0,times = resample_num)

  time <- Time[m:1, ]
  delta <- Delta[m:1, ]
  iboot <- 1
  while (iboot <= resample_num) {
    # Generate bootstrap samples
    bootstrap_indices <- sample(1:n, size = n, replace = TRUE)

    trt_boot <- trt[bootstrap_indices]
    Time_boot <- Time[ ,bootstrap_indices,drop = FALSE]
    Delta_boot <- Delta[ ,bootstrap_indices,drop = FALSE]
    time_boot <- time[ ,bootstrap_indices,drop = FALSE]
    delta_boot <- delta[ ,bootstrap_indices,drop = FALSE]
#    cat("dim(cov)=","\n")
#    print(dim(cov))
    cov_boot <- cov[bootstrap_indices, ,drop = FALSE]
    # cat("cov_boot =", "\n")
    # print(cov_boot)

    indices <- order(trt_boot)

    # Reorder bootstrap samples (control before treatment)
    trt_boot <- trt_boot[indices]
    Time_boot <- Time_boot[ ,indices]
    Delta_boot <- Delta_boot[ ,indices]
    time_boot <- time_boot[ ,indices]
    delta_boot <- delta_boot[ ,indices]
    # cat("dim cov_boot =", dim(cov_boot), "\n")
    cov_boot <- cov_boot[indices, ]

    # Set parameters for model fitting calls
    n0 <- sum(trt_boot == 0)
    n1 <- sum(trt_boot == 1)
    m <- nrow(Time_boot)
    tau <- max(Time_boot[m, ])

    # Type function calls on bootstrap data
    z <- NULL
    markov_ind <- FALSE

    if (!is.null(model)) {
      # If a Markov model is specified, call the function to fit a Markov model
      if (model == "markov") {
        z <- markov(n0,n1,m,Time_boot,Delta_boot)
        markov_ind <- TRUE
      }
      # If no model is specified, use the default (KM) model
      else {
        z <- km(n0,n1,m,Time_boot,Delta_boot)
      }
    }

    if (!is.null(z)) {
      # Parameters for type function calls
      dist_state0_boot <- z[[1]]
      dist_state1_boot <- z[[2]]
      untimes0_boot <- z[[3]]
      untimes1_boot <- z[[4]]
      nuntimes0_boot <- z[[5]]
      nuntimes1_boot <- z[[6]]
      max_follow0_boot <- z[[7]]
      max_follow1_boot <- z[[8]]
    }

    # Type function calls
    if (type == "ewt") {
      y[iboot] <- EWT(m,dist_state0_boot,dist_state1_boot,untimes0_boot,untimes1_boot,nuntimes0_boot,nuntimes1_boot)
    }
    else if (type == "ewtr") {
      y[iboot] <- EWTR(n,m,nuntimes0_boot,max_follow0_boot,untimes0_boot,Time_boot,Delta_boot,dist_state0_boot,markov_ind,cov_boot,trt_boot)[[1]]
    }
    else if (type == "rmt") {
      y[iboot] <- RMT(m,rmst_restriction,dist_state0_boot,dist_state1_boot,untimes0_boot,untimes1_boot,nuntimes0_boot,nuntimes1_boot)
    }
    else if (type == "wtr") {
      y[iboot] <- WTR(n,m,tau,trt_boot,Time_boot,Delta_boot)[[1]]
      # cat("wtr =", "\n")
      # print(y[iboot])
    }
    else if (type == "rwtr") {
      y[iboot] <- RWTR(n,m,tau,trt_boot,Time_boot,Delta_boot)[[1]]
    }
    else if (type == "pwt") {
      y[iboot] <- PWT(n,n0,n1,m,Time_boot,Delta_boot,trt_boot,tau)
    }
    else if (type == "max") {
      # Observed Z statistic
      # ewtr_time <- EWTR(n,m,nuntimes0,max_follow0,untimes0,Time,Delta,dist_state0,markov_ind,cov)[[1]]
      # ewtr_time_var <- EWTR(n,m,nuntimes0,max_follow0,untimes0,Time,Delta,dist_state0,markov_ind,cov)[[2]]
      # z_ewtr <- ewtr_time/sqrt(ewtr_time_var)

      # Bootstrap Z statistic
      # ewtr_time_boot <- EWTR(n,m,nuntimes0_boot,max_follow0_boot,untimes0_boot,Time_boot,Delta_boot,dist_state0_boot,markov_ind,cov_boot)[[1]]
      # ewtr_time_var_boot <- EWTR(n,m,nuntimes0_boot,max_follow0_boot,untimes0_boot,Time_boot,Delta_boot,dist_state0_boot,markov_ind,cov_boot)[[2]]
      z_ewtr_boot <- EWTR(n,m,nuntimes0_boot,max_follow0_boot,untimes0_boot,Time_boot,Delta_boot,dist_state0_boot,markov_ind,cov_boot,trt_boot)[[3]]
      z_comp_boot <- COMP(n,Time_boot,Delta_boot,cov_boot,trt_boot)[[1]]
      y[iboot] <- max((z_ewtr_boot - z_ewtr),(z_comp_boot - z_comp))
    }
    else {
      stop(paste("Invalid type:", type, "- Please specify one of 'ewt', 'ewtr', 'rmt', 'wtr', 'rwtr', 'pwt', 'max'"))
    }
    iboot <- iboot + 1
  }
  return(y)
}
