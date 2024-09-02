#' Win time ratio
#'
#' This function calculates the ratio of losses to wins on treatment. It iterates through all pairs of treatment and control patients and uses
#' their win time difference as the deciding factor of a win or loss.
#'
#' @param n The total number of trial participants.
#' @param m The number of events in the hierarchy.
#' @param tau The maximum follow up time (days).
#' @param tg A numeric vector containing treatment arm indicators (1 for treatment, 0 for control).
#' @param Time A m x n matrix of event times (days). Rows should represent events and columns should represent participants. Event rows should be
#' in increasing order of clinical severity.
#' @param Delta A m x n matrix of event indicators. Rows should represent events and columns should represent participants. Event rows should
#' be in increasing order of clinical severity.
#' @return A list containing: The ratio of losses to wins on treatment, the total number of wins, and the total number of losses.

# ------------------------------------
# Win time ratio
# ------------------------------------
WTR <- function(n,m,tau,tg,Time,Delta) {
  win <- matrix(0,nrow=n,ncol=n)
  loss <- matrix(0,nrow=n,ncol=n)

  wintot <- 0
  losstot <- 0

  total <- 0
  npairs <- 0

  # Start main loop
  for (i in 1:n) {
    for (j in 1:n) {
      if ((tg[i] == 1) && (tg[j] == 0))  {
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

        # Sorted array of event times
        order <- setEventTimes(m,delta0,delta1,time0,time1,follow)

        # Calculate integral
        integral <- getWintimeIntegral(m,order,time0,time1,delta0,delta1)

        # Update win and loss matrices
        if (integral > 0.0) {
          win[i, j] <- 1
        }
        if (integral < 0.0) {
          loss[i, j] <- 1
        }
        npairs <- npairs + 1
        total <- total + integral

        # Update win and loss totals
        if (integral > 0.0) {
          wintot <- wintot + 1
        }
        if (integral < 0.0) {
          losstot <- losstot + 1
        }
      }
    }
  }
  return(list(losstot/wintot,wintot,losstot))
}
