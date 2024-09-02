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
#' @return The z-statistic of the treatment effect from the Cox Model fit.

# --------------------------------
# Composite Analysis
# --------------------------------
COMP <- function(n,Time,Delta,cov,trt) {
  TIME_COMP <- numeric(n)
  DELTA_COMP <- numeric(n)
  for (i in 1:n) {
    TIME_COMP[i] <- min(Time[ ,i])
    DELTA_COMP[i] <- max(Delta[ ,i])
  }

  if (!is.null(cov)) {
    fit=survival::coxph(survival::Surv(TIME_COMP,DELTA_COMP)~trt+cov)
  }
  else {
    fit <- survival::coxph(survival::Surv(TIME_COMP,DELTA_COMP)~trt)
  }
  z_COMP <- -1*coef(fit)[1]/sqrt(vcov(fit)[1,1])
  return(z_COMP)
}
