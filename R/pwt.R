#' Pairwise win time
#'
#' This function calculates the sum of each pair's win time difference divided by the total number of pairs.
#'
#' @param n The total number of trial participants.
#' @param n0 The number of control arm patients.
#' @param n1 The number of treatment arm patients.
#' @param m The number of events in the hierarchy.
#' @param Time A m x n matrix of event time (days). Rows should represent events and columns should represent participants. Event rows should
#' be in increasing order of clinical severity.
#' @param Delta A m x n matrix of event indicators. Rows should represent events and columns should represent participants. Event rows
#' should be in increasing order of clinical severity.
#' @param tg A numeric vector containing treatment arm indicators (1 for treatment, 0 for control).
#' @param tau The maximum follow up time (days).
#' @return The pairwise win time.


# -----------------------------------------
# Pairwise win time
# -----------------------------------------
PWT <- function(n,n0,n1,m,Time,Delta,tg,tau) {
  total <- 0
  npairs <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      if (tg[i] == 1 && tg[j] == 0) {
        # Set max follow-up time for the control person to their event m time
        follow0 <- Time[m,j]
        # If event m is observed in person j
        if (Delta[m,j] == 1)  {
          follow0 <- tau
        }
        # Set max follow-up time for the treatment person to their event m time
        follow1 <- Time[m,i]

        # If event m is observed in person i
        if (Delta[m,i] == 1)  {
          follow1 <- tau
        }
        follow <- min(follow0,follow1)

        # Temp variables for times and deltas
        time0 <- Time[1:m,j]
        time1 <- Time[1:m,i]
        delta0 <- Delta[1:m,j]
        delta1 <- Delta[1:m,i]

        order <- setEventTimes(m,delta0,delta1,time0,time1,follow)
        total <- total + getWintimeIntegral(m,order,time0,time1,delta0,delta1)
        npairs <- npairs + 1
      }
    }
  }
  return(total/npairs)
}
