#' Time Restricted Pairwise win time
#'
#' This function calculates the sum of each pair's win time difference (truncated at the user-specified time_restriction (days)) divided by the total number of pairs.
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
#' @param time_restriction The time restriction (days) for calculation.
#' @return The pairwise win time, and the components of the pairwise win time.


# -----------------------------------------
# Time Restricted Pairwise win time
# -----------------------------------------
RPWT <- function(n,n0,n1,m,Time,Delta,tg,tau,time_restriction) {
  total <- 0
  npairs <- 0
  count <- 0
  icount <- 0
  WTD <- rep(0,times=n0*n1)
  dim(WTD) <- c(n1,n0)
  tgroup <- c(1:n)[tg==1]
  cgroup <- c(1:n)[tg==0]
  components <- rep(0,m)
  for (i in tgroup) {
    icount <- icount+1
    jcount <- 0
    for (j in cgroup) {
      jcount <- jcount+1
#  for (i in 1:n) {
#    for (j in 1:n) {
#      if (tg[i] == 1 && tg[j] == 0) {
        count <- count+1
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

        # if (icount==1 & jcount==2) {
        #   cat('icount=',icount,'\n')
        #   cat('time1=',time1,'\n')
        #   cat('delta1=',delta1,'\n')
        #   cat('jcount=',jcount,'\n')
        #   cat('time0=',time0,'\n')
        #   cat('delta0=',delta0,'\n')
        # }

        order <- setEventTimes(m,delta0,delta1,time0,time1,follow)

        # if (icount==1 & jcount==2) {
        #   cat('order=',order,'\n')
        # }

        temp <- getWintimeIntegral_rest(m,order,time0,time1,delta0,delta1,time_restriction)
        WTD[icount,jcount] <- temp[[1]]
        comp <- temp[[2]]
        total <- total + WTD[icount,jcount]
        components <- components + comp
        npairs <- npairs + 1
#      }
    }
  }
#  cat('final total=',total,'\n')
#  cat('WTD[1,2]=','\n')
#  print(WTD[1,2])
#  cat('npairs=',npairs,'\n')

  return(list(total/npairs,components/npairs))
}
