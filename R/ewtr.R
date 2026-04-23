#' Expected win time against reference
#'
#' Calculates the control group state space probabilities using a Markov model (recommended) or a Kaplan-Meier model. This function uses these
#' probabilities to compare each participant's clinical state to a distribution of control group states.
#'
#' @param n The total number of trial participants.
#' @param m The number of events in the hierarchy.
#' @param nunique The number of unique control group event times (returned from wintime::markov() or wintime::km()).
#' @param maxfollow The max control group follow up time (days) (returned from wintime::markov() or wintime::km()).
#' @param untimes A vector containing unique control group event times (days) (returned from wintime::markov() or wintime::km()).
#' @param Time A m x n matrix of event times (days). Rows should represent events and columns should represent participants. Rows should be
#' in increasing order of clinical severity.
#' @param Delta A m x n matrix of event indicators Rows should represent events and columns should represent participants. Rows should be
#' in increasing order of clinical severity.
#' @param dist A matrix of control group state probabilities (returned from wintime::markov() or wintime::km()).
#' @param markov_ind An indicator of the model type used (1 for Markov, 0 for Kaplan-Meier).
#' @param cov A n x p matrix of covariate values, where p is the number of covariates.
#' @param trt A vector of length n containing treatment arm indicators (1 for treatment, 0 for control).
#' @return A list containing: The estimated treatment effect from the linear regression model, the variance, the Z-statistic,
#' the components of the treatment effect, and the variance of the components.

# ----------------------------------------
# Expected win time against reference
# ----------------------------------------
EWTR <- function(n,m,nunique,maxfollow,untimes,Time,Delta,dist,markov_ind,cov,trt) {
  time <- Time[m:1, ]
  delta <- Delta[m:1, ]
  components <- rep(NA,m)
  components_var <- rep(NA,m)

  # cat('From ewtr: # unique control event times=',nunique,'\n')
  # cat('From ewtr: unique control times=',untimes,'\n')
  # cat('From ewtr: maximum event time for control arm=',maxfollow,'\n')
  # cat('From ewtr: control arm dist=','\n')
  # print(dist)

  # Check for extending comparisons when control arm is 100% in terminal state
  extend=0
  if (dist[m+1,nunique]==1) {extend=1}

  # cat('-------------------------------------------','\n')
  # cat('extend=',extend,'\n')
  # cat('-------------------------------------------','\n')

  # cat('Before addtime, control arm dist=','\n')
  # print(tdist)
  # cat('Before addtimes, # of tuntimes=',tnunique,'\n')
  # cat('Before addtimes, tuntimes=','\n')
  # print(tuntimes)

  # Initialize temporary variables
  tdist <- matrix(0,nrow=m+1,ncol=nunique+m)
  ewtr <- rep(0,times=n)
  ewtr_components <- matrix(0,nrow=m,ncol=n)

  # Start main loop
  for (i in 1:n) {
     # cat('-------------------------------------------','\n')
     # cat('-------------------------------------------','\n')
     # cat('i=',i,'\n')

    # Create temp variables
    tuntimes <- rep(0,times=nunique+m)
    tnunique <- nunique
    for (j in 1:nunique) {
      tuntimes[j] <- untimes[j]
    }

    # Copy control distributions
    for (event in 1:(m+1)) {
      for (t in 1:min(nunique,ncol(dist))) {
        tdist[event,t] <- dist[event,t]
      }
    }

     # if (i==20) {
     #   cat('Before addtime, control arm dist=','\n')
     #   print(tdist)
     #   cat('Before addtimes, # of tuntimes=',tnunique,'\n')
     #   cat('Before addtimes, tuntimes=','\n')
     #   print(tuntimes)
     # }

    # Create addtime matrix
    addtime <- matrix(data=NA,nrow=m,ncol=n)
    # addtime1=(!(dataset1$time1 %in% unique_event_times0) & dataset1$time1 <= unique_event_times0[nunique_event_times0])
    # addtime2=(!(dataset1$time2 %in% unique_event_times0) & dataset1$delta2==1 & dataset1$time2 <= unique_event_times0[nunique_event_times0])
    # addtime3=(!(dataset1$time3 %in% unique_event_times0) & dataset1$delta3==1 & dataset1$time3 <= unique_event_times0[nunique_event_times0])

    #addtime[1, ] <- (!(time[1, ] %in% untimes) & time[1, ] <= untimes[nunique])
    addtime[1, ] <- (!(time[1, ] %in% tuntimes) & (time[1,] <= tuntimes[tnunique] | extend==1))
    for (k in 2:m) {
      #addtime[k, ] <- (!(time[k, ] %in% untimes) & delta[k, ] == 1 & time[k, ] <= tuntimes[tnunique])
      addtime[k, ] <- (!(time[k, ] %in% tuntimes) & delta[k, ] == 1 & (time[k,] <= tuntimes[tnunique] | extend==1))
    }
    # cat("-----------------------------------------------", "\n")
    # cat("dim addtime =", dim(addtime), "\n")
    # cat("addtime =","\n")
    # print(addtime)
    # Add events

    # cat('addtime (rows in reverse order of event types, columns are subjects)=','\n')
    # print(addtime)

    for (event in 1:m) {
#      if (i==20) {cat('event=',event,'\n')}
      if (addtime[event,i] == TRUE) {
        jstop <- 0
        for (j in 1:tnunique) {
          if (tuntimes[j] < time[event,i]) {
            jstop <- jstop + 1
          }
          else {
            break
          }
        }
        if(jstop == 0) {
          jstop <- 1
        }
 #       if (i==20) {cat('Addtime==TRUE, jstop=',jstop,'\n')}
        #tdist[,tnunique] <- 0
        tnunique <- tnunique + 1
        # if (i==18) {cat('tnunique=',tnunique,'\n')}
        if (jstop==1) {
          tdist[,tnunique] <- c(1,rep(0,times=m))
        } else {
          tdist[,tnunique] <- tdist[,jstop]
        }
        tuntimes[tnunique] <- time[event,i]
        j <- tnunique

#        if (i==20) {cat('tuntimes[tnunique]=',tuntimes[tnunique],'\n')}
#        if (i==20) {cat('j=',j,'\n')}
        #while(j > jstop) {
        #  tdist[,j] <- tdist[,(j-1)]
        #  if (j > jstop + 1) {
            # swap indices j and (j-1)
        #    temp <- tuntimes[j]
        #    tuntimes[j] <- tuntimes[j-1]
        #    tuntimes[j-1] <- temp
        #  }
        #  j <- j - 1
        #}
 #        if (i==20) {cat('tuntimes[j-1]=',tuntimes[j-1],'\n')}
        while(tuntimes[j] < tuntimes[j-1]) {
 #--------------------------------------------------
 # swap indices j and (j-1)
 #--------------------------------------------------
 #         if (i==20) {cat('swap indicies j & j-1=',j,j-1,'\n')}
          temp <- tuntimes[j]
          tuntimes[j] <- tuntimes[j-1]
          tuntimes[j-1] <- temp
          vtemp <- tdist[,j]
          tdist[,j] <- tdist[,j-1]
          tdist[,j-1] <- vtemp
          j <- j - 1
          if (j==jstop) {break}
#          if (i==20) {cat('new index j=',j,'\n')}
        }

      }
    }
    tuntimes=tuntimes[1:tnunique]

    # if (i==20) {
    #   cat('After addtime, control arm dist=','\n')
    #   print(tdist)
    #   cat('After addtimes, # of tuntimes=',tnunique,'\n')
    #   cat('After addtimes, tuntimes=','\n')
    #   print(tuntimes)
    # }

    # cat(' # nunique event times after addition of current subject=',tnunique,'\n')
    # cat('unique event times after addition of current subject=',tuntimes,'\n')
    # cat('control arm dist=','\n')
    # print(tdist)


    # Set jmax
    jmax <- 0
    for (j in 1:tnunique) {
      if (tuntimes[j] < time[1,i]) {
        jmax <- jmax + 1
      }
    }
    #if (markov_ind == FALSE) {
    #  jmax <- min(maxfollow,jmax)
    #}

    if (delta[1,i] == 1) {
      jmax <- tnunique - 1
    }

    # cat('jmax=',jmax,'\n')
    # cat('---------------------------------------','\n')

    # Set state
    if (jmax != 0) {
      for (j in 1:jmax) {

        #cat('j=',j,'\n')

        state <- 0
        for (current_state in m:1) {
          temp_time_index <- m - current_state + 1
          if (time[temp_time_index,i] <= tuntimes[j] && delta[temp_time_index,i] == 1) {
            state <- current_state
            break
          }
        }

        #cat('state=',state,'\n')

        # Calculate ewtr

        # Calculate wins
        for (state_num in 0:(m-1)) {
          if (state_num == state) {
            # Add probabilities from higher states
            for (k in (state_num+1):m) {
              ewtr[i] <- ewtr[i] + tdist[k+1,j] * (tuntimes[j+1] - tuntimes[j])
              ewtr_components[k,i] <- ewtr_components[k,i] + tdist[k+1,j] * (tuntimes[j+1]-tuntimes[j])
            }
            break
          }
        }

        # Calculate Losses
        for (current in 1:m) {
          if (current == state) {
            # Subtract probabilities from lower states
            for (k in 1:current) {
              ewtr[i] <- ewtr[i] - tdist[k,j] * (tuntimes[j+1] - tuntimes[j])
              ewtr_components[current,i] <- ewtr_components[current,i] - tdist[k,j] * (tuntimes[j+1]-tuntimes[j])
            }
            break
          }
        }
        #cat('ewtr[',i,']=',ewtr[i],'\n')
      }
    }
  }

  # Get treatment estimate and variance for Z statistic
  fit_comp <- vector("list",m)
  if (!is.null(cov)) {
    fite=lm(ewtr~trt+cov)
    for (k in 1:m) {
      outcome <- ewtr_components[k,]
      dim(outcome) <- c(n)
      fit_comp[[k]] <- lm(outcome~trt+cov)
    }
  }
  else {
    fite <- lm(ewtr~trt)
    for (k in 1:m) {
      outcome <- ewtr_components[k,]
      dim(outcome) <- c(n)
      fit_comp[[k]] <- lm(outcome~trt)
    }
  }
  ewtr_time=coef(fite)[2]
  ewtr_time_var=vcov(fite)[2,2]
  z_ewtr <- ewtr_time/sqrt(ewtr_time_var)
  for (k in 1:m) {
    components[k] <- coef(fit_comp[[k]])[2]
    components_var[k] <- vcov(fit_comp[[k]])[2,2]
  }

  return(list(ewtr_time,ewtr_time_var,z_ewtr,components,components_var))
}
