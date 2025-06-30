#' Restricted mean survival in favor of treatment
#'
#' Calculates the state space probabilities using a Kaplan-Meier model (recommended) or a Markov model. This function uses these probabilities
#' to compare both arms and calculate the expected win time of the treatment arm up to a given time point.
#'
#' @param m The number of events in the hierarchy.
#' @param time_restriction The cutoff time point (days) for the calculation.
#' @param dist_state0 A matrix of control arm state probabilities (returned from wintime::km() or wintime::markov()).
#' @param dist_state1 A matrix of treatment arm state probabilities (returned from wintime::km() or wintime::markov()).
#' @param unique_event_times0 A vector of unique control arm event times (days) (returned from wintime::km() or wintime::markov()).
#' @param unique_event_times1 A vector of unique treatment arm event times (days) (returned from wintime::km() or wintime::markov()).
#' @param nunique_event_times0 The number of unique control arm event times (returned from wintime::km() or wintime::markov()).
#' @param nunique_event_times1 The number of unique treatment arm event times (returned from wintime::km() or wintime::markov()).
#' @return A list containing: The restricted mean survival in favor of the treatment arm, the components of the treatment effect.


# ---------------------------------------------------
# Restricted mean survival in favor of treatment
# ---------------------------------------------------
RMT <- function(m,time_restriction,dist_state0,dist_state1,unique_event_times0,unique_event_times1,nunique_event_times0,nunique_event_times1) {
#  cat("Starting RMT", "\n")
#  cat("time_restriction=",time_restriction,"\n")
  components <- rep(0,m)
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

#  cat("nunique_event_times=",nunique_event_times,"\n")
#  cat("unique_event_times=","\n")
#  print(unique_event_times)
#  cat("dist_state0=","\n")
#  print(dist_state0)
#  cat("new_dist_state1=","\n")
#  print(new_dist_state1)


  # RESTRICT LIST TO SMALLER OF MAX FOLLOW-UPS FOR ARMS
  nunique_event_times=length(unique_event_times[unique_event_times <= unique_event_times1[nunique_event_times1]])

#  cat(" updated nunique_event_times=",nunique_event_times,"\n")

  rmst_time <- 0
  j <- 1

  # Loop 3
  while (j < nunique_event_times) {
    if (unique_event_times[j+1] <= time_restriction) {
#      cat("j=",j,"\n")
      # Add RMST
      for (event_num in 1:m) {
        rmst_time <- rmst_time + new_dist_state1[event_num,j] * sum(dist_state0[((event_num+1):(m+1)),j]) * (unique_event_times[j+1]-unique_event_times[j])
        for (k in (event_num+1):(m+1)) {
          components[k-1] <- components[k-1] + new_dist_state1[event_num,j] * dist_state0[k,j] * (unique_event_times[j+1]-unique_event_times[j])
        }
#        if (j==2) {
#          cat("event_num for treated for win=",event_num,"\n")
#          cat("new_dist_state1[event_num,j]=",new_dist_state1[event_num,j],"\n")
#          cat("sum(dist_state0[((event_num+1):(m+1)),j])=",sum(dist_state0[((event_num+1):(m+1)),j]),"\n")
#          cat("updated RMT=",rmst_time,"\n")
#        }
      }

      # Subtract RMST
      for (event_num in 2:(m+1)) {
        rmst_time <- rmst_time - new_dist_state1[event_num,j] * sum(dist_state0[(1:(event_num-1)),j]) * (unique_event_times[j+1]-unique_event_times[j])
        components[event_num-1] <- components[event_num-1] - new_dist_state1[event_num,j] * sum(dist_state0[(1:(event_num-1)),j]) * (unique_event_times[j+1]-unique_event_times[j])
#        if (j==2) {
#          cat("event_num for treated for loss=",event_num,"\n")
#          cat("updated RMT=",rmst_time,"\n")
#        }
      }
#      cat("updated RMT=",rmst_time,"\n")

    } else {
      if (unique_event_times[j] <= time_restriction) {
#        cat("j=",j,"\n")
        # Add RMST
        for (event_num in 1:m) {
          rmst_time <- rmst_time + new_dist_state1[event_num,j] * sum(dist_state0[((event_num+1):(m+1)),j]) * (time_restriction-unique_event_times[j])
         for (k in (event_num+1):(m+1)) {
            components[k-1] <- components[k-1] + new_dist_state1[event_num,j] * dist_state0[k,j] * (time_restriction-unique_event_times[j])
          }
#           if (j==2) {
#            cat("event_num for treated for win=",event_num,"\n")
#            cat("new_dist_state1[event_num,j]=",new_dist_state1[event_num,j],"\n")
#            cat("sum(dist_state0[((event_num+1):(m+1)),j])=",sum(dist_state0[((event_num+1):(m+1)),j]),"\n")
#            cat("updated RMT=",rmst_time,"\n")
#          }
        }

        # Subtract RMST
        for (event_num in 2:(m+1)) {
          rmst_time <- rmst_time - new_dist_state1[event_num,j] * sum(dist_state0[(1:(event_num-1)),j]) * (time_restriction-unique_event_times[j])
          components[event_num-1] <- components[event_num-1] - new_dist_state1[event_num,j] * sum(dist_state0[(1:(event_num-1)),j]) * (time_restriction-unique_event_times[j])
          #        if (j==2) {
          #          cat("event_num for treated for loss=",event_num,"\n")
          #          cat("updated RMT=",rmst_time,"\n")
          #        }
        }
#        cat("updated RMT=",rmst_time,"\n")
      }
    }
    j <- j + 1
  }
#  cat("final RMT=",rmst_time,"\n")
#  cat("final components of RMT=",components,"\n")
  return(list(rmst_time,components))
}
