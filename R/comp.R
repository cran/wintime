#' Run composite analysis
#'
#' This function fits a Cox Model to time-to-event data and calculates the z statistic. In the wintime package, this function is used for the
#' EWTR-composite max test (MAX) method.
#'
#' @param n The total number of trial participants.
#' @param Time A m x n matrix of event times (days), where m is the number of events in the hierarchy. Rows should represent events and columns
#' should represent participants. Event rows should be in increasing order of clinical severity.
#' @param Delta A m x n matrix of event indicators, where m is the number of events in the hierarchy. Rows should represent events and
#' columns should represent participants. Event rows should be in increasing order of clinical severity.
#' @param cov A n x p matrix of covariate values, where p is the number of covariates. Rows should represent participants and columns
#' should represent covariate values.
#' @param trt A vector of length n containing treatment arm indicators (1 for treatment, 0 for control).
#' @return A list containing: The z-statistic of the treatment effect from the Cox Model fit, the treatment effect estimate, the variance of the treatment effect estimate, the p-value for treatment effect.

# --------------------------------
# Composite Analysis
# --------------------------------
COMP <- function(n,Time,Delta,cov,trt) {
  TIME_COMP <- numeric(n)
  DELTA_COMP <- numeric(n)
  m=length(Time[,1])

  for (i in 1:n) {
    TIME_COMP[i] <- Time[m,i]
    TIME_COMP[i] <- min(TIME_COMP[i],Time[Delta[,i]==1,i])
#    TIME_COMP[i] <- min(Time[ ,i])
    DELTA_COMP[i] <- max(Delta[ ,i])
  }

  if (!is.null(cov)) {
    fit=survival::coxph(survival::Surv(TIME_COMP,DELTA_COMP)~trt+cov)
  }
  else {
    fit <- survival::coxph(survival::Surv(TIME_COMP,DELTA_COMP)~trt)
  }
  z_COMP <- -1*coef(fit)[1]/sqrt(vcov(fit)[1,1])
  est_COMP <- exp(coef(fit)[1])
  var_COMP=vcov(fit)[1,1]
  p_COMP <- 2*(1-pnorm(abs(z_COMP),mean=0,sd=1))
  return(list(z_COMP,est_COMP,var_COMP,p_COMP))
}
