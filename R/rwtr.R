#' Restricted win time ratio
#'
#' This function calculates the ratio of losses to wins on treatment. It iterates through all pairs of treatment and control patients and uses
#' their time-to-death (or most severe clinical event) to determine a win or loss. If death is inconclusive, then a winner is determined based
#' on wintime.
#'
#' @param n The total number of trial participants.
#' @param m The number of events in the hierarchy.
#' @param tau The maximum follow up time (days).
#' @param tg A numeric vector containing treatment arm indicators (1 for treatment, 0 for control).
#' @param Time A m x n matrix of event times (days), where m is the number of events in the hierarchy, and n is the total number of trial participants.
#' Rows should represent events and columns should represent participants. Event rows should be in increasing order of clinical severity.
#' @param Delta A m x n matrix of event indicators, where m is the number of events in the hierarchy, and n is the total number of trial participants.
#' Rows should represent events and columns should represent participants. Event rows should be in increasing order of clinical severity.
#' @return A list containing: The ratio of losses to wins on treatment, the total number of wins, and the total number of losses.

# ---------------------------------------------------
# Restricted win time ratio
# ---------------------------------------------------
RWTR <- function(n,m,tau,tg,Time,Delta) {
  rwin <- matrix(0,nrow=n,ncol=n)
  rloss <- matrix(0,nrow=n,ncol=n)

  rwintot <- 0
  rlosstot <- 0

  for (i in 1:n) {
    for (j in 1:n) {
      if (tg[i] == 1 && tg[j] == 0) {
        # If event m favors the treatment person
        if ((Delta[m,j] == 1 && Time[m,i] > Time[m,j]) ||
            (Delta[m,j] == 1 && Time[m,i] == Time[m,j] && Delta[m,i] == 0)) {
          # Win for the treatment person
          rwin[i, j] <- 1
          rwintot <- rwintot + 1
        }

        # If event m favors the control person
        else if ((Delta[m,i] == 1 && Time[m,j] > Time[m,i]) ||
                 (Delta[m,i] == 1 && Time[m,j] == Time[m,i] && Delta[m,j] == 0)) {
          # Loss for the treatment person
          rloss[i, j] <- 1
          rlosstot <- rlosstot + 1
        }
        else {
          # Pick a winner based on win time

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
          temp <- getWintimeIntegral(m,order,time0,time1,delta0,delta1)
          integral <- temp[[1]]

          if (integral > 0)
          {
            rwin[i, j] <- 1
            rwintot <- rwintot + 1
          }
          if(integral < 0)
          {
            rloss[i, j] <- 1
            rlosstot <- rlosstot + 1
          }
        }
      }
    }
  }
  return(list(rlosstot/rwintot,rwintot,rlosstot))
}
