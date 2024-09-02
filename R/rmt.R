#' Restricted mean survival in favor of treatment
#'
#' Calculates the state space probabilities using a Kaplan-Meier model (recommended) or a Markov model. This function uses these probabilities
#' to compare both arms and calculate the expected win time of the treatment arm up to a given time point.
#'
#' @param m The number of events in the hierarchy.
#' @param rmst_restriction The cutoff time point (days) for the calculation.
#' @param dist_state0 A matrix of control arm state probabilities (returned from wintime::km() or wintime::markov()).
#' @param dist_state1 A matrix of treatment arm state probabilities (returned from wintime::km() or wintime::markov()).
#' @param unique_event_times0 A vector of unique control arm event times (days) (returned from wintime::km() or wintime::markov()).
#' @param unique_event_times1 A vector of unique treatment arm event times (days) (returned from wintime::km() or wintime::markov()).
#' @param nunique_event_times0 The number of unique control arm event times (returned from wintime::km() or wintime::markov()).
#' @param nunique_event_times1 The number of unique treatment arm event times (returned from wintime::km() or wintime::markov()).
#' @return The restricted mean survival in favor of the treatment arm.

# ---------------------------------------------------
# Restricted mean survival in favor of treatment
# ---------------------------------------------------
RMT <- function(m,rmst_restriction,dist_state0,dist_state1,unique_event_times0,unique_event_times1,nunique_event_times0,nunique_event_times1) {
  #  cat("Starting RMT", "\n")
  unique_event_times=unique_event_times0
  nunique_event_times=length(unique_event_times)
  new_dist_state1 <- matrix(data=0,nrow=m+1,ncol=nunique_event_times)

  # Loop 1
  i <- 1
  while(i <= nunique_event_times) {
    x <- length(unique_event_times1[unique_event_times1 <= unique_event_times[i]])
    for (event_num in 1:(m+1)) {
      if (x > 0) {
        new_dist_state1[event_num,i] <- dist_state1[event_num,x]
      }
      else {
        new_dist_state1[event_num,i] <- 1
      }
    }
    i <- i + 1
  }

  # Loop 2

  # ADD EVENT TIMES FOR TREATED SUBJECTS TO LIST FROM CONTROL SUBJECTS
  #
  i=1
  #    while (i <= 1) {
  while (i <= nunique_event_times1) {
    #cat("loop 2 i =", i, "\n")
    if (!(unique_event_times1[i] %in% unique_event_times) & unique_event_times1[i] <= unique_event_times[nunique_event_times]) {
      #        cat('Enter ADD EVENT time to List','\n')
      jstop=length(unique_event_times[unique_event_times < unique_event_times1[i]])
      jstar=jstop+1
      if (jstop==0) {jstop=1}

      tdist_state0 <- matrix(data = 0, nrow = m+1, ncol = ncol(dist_state0)+1)
      tnew_dist_state1 <- matrix(data = 0, nrow = m+1, ncol = ncol(new_dist_state1)+1)
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

  # RESTRICT LIST TO SMALLER OF MAX FOLLOW-UPS FOR ARMS
  nunique_event_times=length(unique_event_times[unique_event_times <= unique_event_times1[nunique_event_times1]])

  rmst_time <- 0
  j <- 1

  # Loop 3
  while (j < nunique_event_times) {
    if (unique_event_times[j] <= rmst_restriction) {
      # Add RMST
      for (event_num in 1:m) {
        rmst_time <- rmst_time + new_dist_state1[event_num,j] * sum(dist_state0[((event_num+1):(m+1)),j]) * (unique_event_times[j+1]-unique_event_times[j])
      }

      # Subtract RMST
      for (event_num in 2:(m+1)) {
        rmst_time <- rmst_time - new_dist_state1[event_num,j] * sum(dist_state0[(1:(event_num-1)),j]) * (unique_event_times[j+1]-unique_event_times[j])
      }

    }
    j <- j + 1
  }
  return(rmst_time)
}
