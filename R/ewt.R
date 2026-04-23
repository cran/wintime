#' Expected win time
#'
#' Calculates the state space probabilities using a Kaplan-Meier model (recommended) or a Markov model. This function uses these probabilities
#' to compare both arms and calculate the expected win time of the treatment arm.
#'
#' @param m The number of events in the hierarchy.
#' @param dist_state0 A matrix of control arm state probabilities (returned from wintime::km() or wintime::markov()).
#' @param dist_state1 A matrix of treatment arm state probabilities (returned from wintime::km() or wintime::markov()).
#' @param unique_event_times0 A vector of unique control arm event times (days) (returned from wintime::km() or wintime::markov()).
#' @param unique_event_times1 A vector of unique treatment arm event times (days) (returned from wintime::km() or wintime::markov()).
#' @param nunique_event_times0 The number of unique control arm event times (returned from wintime::km() or wintime::markov()).
#' @param nunique_event_times1 The number of unique treatment arm event times (returned from wintime::km() or wintime::markov()).
#' @return A list of the expected win time of the treatment arm, the components of the treatment effect, and the maximum time used in comparisons.

# ------------------------
# Expected win time
# ------------------------
EWT <- function(m,dist_state0,dist_state1,unique_event_times0,unique_event_times1,nunique_event_times0,nunique_event_times1) {
  # cat("Start EWT", "\n")
  components <- rep(0,m)
  unique_event_times=unique_event_times0
  nunique_event_times=length(unique_event_times)
  new_dist_state1 <- matrix(data=0,nrow=m+1,ncol=nunique_event_times)
  max_time <- 0

  max_follow0=unique_event_times0[nunique_event_times0]
  max_follow1=unique_event_times1[nunique_event_times1]

  # cat('From ewt: nunique_event_times0=',nunique_event_times0,'\n')
  # cat('From ewt: unique_event_times0=',unique_event_times0,'\n')
  # cat('From ewt: control group max event time=',max_follow0,'\n')
  # cat('From ewt: nunique_event_times1=',nunique_event_times1,'\n')
  # cat('From ewt: unique_event_times1=',unique_event_times1,'\n')
  # cat('From ewt: trt group max event time=',max_follow1,'\n')
  # cat('------------------------------------------------------------','\n')

  # Loop 1
  i <- 1
  while(i <= nunique_event_times) {
    x <- length(unique_event_times1[unique_event_times1 <= unique_event_times[i]])
    for (event_num in 1:(m+1)) {
      if (x > 0) {
        new_dist_state1[event_num,i] <- dist_state1[event_num,x]
      }
      else {
        if (event_num==1) {new_dist_state1[event_num,i] <- 1}
      }
    }
    i <- i + 1
  }

  # Loop 2

  # ADD EVENT TIMES FOR TREATED SUBJECTS TO LIST FROM CONTROL SUBJECTS
  #
  tdist_state0 <- matrix(data = 0, nrow = m+1, ncol = ncol(dist_state0))
  #cat("dim tdist_state0 =", dim(tdist_state0), "\n")
  tnew_dist_state1 <- matrix(data = 0, nrow = m+1, ncol = ncol(new_dist_state1))
  i=1
  #    while (i <= 1) {
  while (i <= nunique_event_times1) {
    #cat("loop 2 i =", i, "\n")
    #if (!(unique_event_times1[i] %in% unique_event_times) & unique_event_times1[i] <= unique_event_times[nunique_event_times]) {
    if (!(unique_event_times1[i] %in% unique_event_times)) {
      #        cat('Enter ADD EVENT time to List','\n')
      jstop=length(unique_event_times[unique_event_times < unique_event_times1[i]])
      jstar=jstop+1
      if (jstop==0) {jstop=1}

      tdist_state0 <- matrix(data = 0, nrow = m+1, ncol = ncol(tdist_state0)+1)
      tnew_dist_state1 <- matrix(data = 0, nrow = m+1, ncol = ncol(tnew_dist_state1)+1)
      # cat("dim tdist_state0 =", dim(tdist_state0), "\n")



      # Loop through each row and append 0 to the end
      for (k in 1:(m+1)) {
        current_dist_state0 <- c(dist_state0[k,], 0)
        tdist_state0[k, ] <- current_dist_state0

        current_dist_state1 <- c(new_dist_state1[k,], 0)
        tnew_dist_state1[k,] <- current_dist_state1
      }

      # Assign back to original matrices
      dist_state0 <- tdist_state0
      new_dist_state1 <- tnew_dist_state1

      j <- nunique_event_times + 1
      while (j > jstop) {
        for (event_num in 1:(m+1)) {
          dist_state0[event_num,j] <- dist_state0[event_num,j-1]
          new_dist_state1[event_num,j] <- new_dist_state1[event_num,j-1]
        }
        j <- j - 1
      }

      for (event_num in 1:(m+1)) {
        new_dist_state1[event_num,jstar] <- dist_state1[event_num,i]
      }

      unique_event_times=c(unique_event_times,unique_event_times1[i])
      perm=order(unique_event_times)
      unique_event_times=unique_event_times[perm]
      nunique_event_times=length(unique_event_times)
    }
    i=i+1
  }

  # cat('-------------------------------------','\n')
  # cat('Before restriction unique _event_times=','\n')
  # print(unique_event_times)
  # cat('-------------------------------------','\n')
  # cat('dist_state0[k,] for being in state k-1=','\n')
  # print(dist_state0)
  # cat('-------------------------','\n')
  # cat('new_dist_state1[k,] for being in state k-1=','\n')
  # print(new_dist_state1)
  # cat('-------------------------','\n')

  # RESTRICT LIST TO SMALLER OF MAX FOLLOW-UPS FOR ARMS
  # IF ONE ARM ENTERS 100% TERMINAL STATE, THEN CONTINUE TO LARGER FOLLOW-UP TIME
  #nunique_event_times=length(unique_event_times[unique_event_times <= unique_event_times1[nunique_event_times1]])
  if (max_follow0 <= max_follow1) {
    max_time=max_follow0
    #cat('max_time=',max_time,'\n')
    j=length(unique_event_times[unique_event_times <= max_time])
    #cat('j=',j,'\n')
    #cat('dist_state0[m+1,j]=',dist_state0[m+1,j],'\n')
    if (dist_state0[m+1,j]==1) {max_time=max_follow1}
    #cat('max_time=',max_time,'\n')
  } else {
    max_time=max_follow1
    j=length(unique_event_times[unique_event_times <= max_time])
    if (new_dist_state1[m+1,j]==1) {max_time=max_follow0}
  }
  #max_time=min(max_follow0,max_follow1)
  nunique_event_times=length(unique_event_times[unique_event_times <= max_time])

  # cat('-------------------------------------','\n')
  # cat('After restriction, # of times used in calculating EWT=',nunique_event_times,'\n')
  # cat('After restriction, Largest time used in calculating EWT=',unique_event_times[nunique_event_times],'\n')
  # cat('-------------------------------------','\n')

  # CALCULATE NET WINTIME
  # COMPARES DISTRIBUTIONS ACROSS TIMES IN COMBINED LIST

  ewt_time=0
  j=1

  # Loop 3
  while (j < nunique_event_times) {
    # Add wintime
    for (event_num in 1:m) {
      ewt_time <- ewt_time + new_dist_state1[event_num,j] * sum(dist_state0[((event_num+1):(m+1)),j]) * (unique_event_times[j+1]-unique_event_times[j])
      for (k in (event_num+1):(m+1)) {
        components[k-1] <- components[k-1] + new_dist_state1[event_num,j] * dist_state0[k,j] * (unique_event_times[j+1]-unique_event_times[j])
      }
    }

    # Subtract Wintime
    for (event_num in 2:(m+1)) {
      ewt_time <- ewt_time - new_dist_state1[event_num,j] * sum(dist_state0[(1:(event_num-1)),j]) * (unique_event_times[j+1]-unique_event_times[j])
      components[event_num-1] <- components[event_num-1] - new_dist_state1[event_num,j] * sum(dist_state0[(1:(event_num-1)),j]) * (unique_event_times[j+1]-unique_event_times[j])
    }
    j=j+1
  }

  # cat('-------------------------------------','\n')
  # cat('-------------------------------------','\n')
  # cat('Final ewt_time=',ewt_time,'\n')
  # cat('-------------------------------------','\n')

  return(list(ewt_time,components,max_time))
}
