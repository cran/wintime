#' Expected win time against trial population
#'
#' Calculates the combined arm state space probabilities using a Markov model or a Kaplan-Meier model (recommended). This function uses these
#' probabilities to compare each participant's clinical state to a distribution of combined arm states.
#'
#' @param n The total number of trial participants.
#' @param m The number of events in the hierarchy.
#' @param nunique The number of unique combined arm event times (returned from wintime::markov() or wintime::km()).
#' @param maxfollow The max combined arm follow up time (days) (returned from wintime::markov() or wintime::km()).
#' @param untimes A vector containing unique combined arm event times (days) (returned from wintime::markov() or wintime::km()).
#' @param Time A m x n matrix of event times (days). Rows should represent events and columns should represent participants. Rows should be
#' in increasing order of clinical severity.
#' @param Delta A m x n matrix of event indicators Rows should represent events and columns should represent participants. Rows should be
#' in increasing order of clinical severity.
#' @param dist A matrix of combined arm state probabilities (returned from wintime::markov() or wintime::km()).
#' @param markov_ind An indicator of the model type used (1 for Markov, 0 for Kaplan-Meier).
#' @param cov A n x p matrix of covariate values, where p is the number of covariates.
#' @param trt A vector of length n containing treatment arm indicators (1 for treatment, 0 for control).
#' @param time_restriction The time restriction (days) for calculation.
#' @return A list containing: The estimated treatment effect from the linear regression model, the variance, the Z-statistic, the components of the treatment effect, and the variance of the components.

# -------------------------------------------
# Expected win time against trial population
# -------------------------------------------
REWTP <- function(n,m,nunique,maxfollow,untimes,Time,Delta,dist,markov_ind,cov,trt,time_restriction) {
  time <- Time[m:1, ]
  delta <- Delta[m:1, ]
  components <- rep(NA,m)
  components_var <- rep(NA,m)

  # cat('----------------------------------------------','\n')
  # cat('----------------------------------------------','\n')
  # cat('nunique=',nunique,'\n')
  # cat('untimes=',untimes,'\n')
  # cat('time[,3]=',time[,3],'\n')
  # cat('delta[,3]=',delta[,3],'\n')
  # cat('----------------------------------------------','\n')

  # Initialize temporary variables
  tuntimes <- numeric(nunique+3)
  tdist <- matrix(0,nrow=m+1,ncol=nunique+m)
  rewtp <- rep(0,n)
  rewtp_components <- matrix(0,nrow=m,ncol=n)

  # Start main loop
  for (i in 1:n) {
    # Create temp variables
    tnunique <- nunique
    for (j in 1:nunique) {
      tuntimes[j] <- untimes[j]
    }

    # Copy combined arm distributions
    for (event in 1:(m+1)) {
      for (t in 1:min(nunique,ncol(dist))) {
        tdist[event,t] <- dist[event,t]
      }
    }

    # Create addtime matrix
    addtime <- matrix(data=NA,nrow=m,ncol=n)
    # addtime1=(!(dataset1$time1 %in% unique_event_times0) & dataset1$time1 <= unique_event_times0[nunique_event_times0])
    # addtime2=(!(dataset1$time2 %in% unique_event_times0) & dataset1$delta2==1 & dataset1$time2 <= unique_event_times0[nunique_event_times0])
    # addtime3=(!(dataset1$time3 %in% unique_event_times0) & dataset1$delta3==1 & dataset1$time3 <= unique_event_times0[nunique_event_times0])
    addtime[1, ] <- (!(time[1, ] %in% untimes) & time[1, ] <= untimes[nunique])
    for (k in 2:m) {
      addtime[k, ] <- (!(time[k, ] %in% untimes) & delta[k, ] == 1 & time[k, ] <= untimes[nunique])
    }
    # cat("-----------------------------------------------", "\n")
    # cat("dim addtime =", dim(addtime), "\n")
    # cat("addtime =","\n")
    # print(addtime)
    # Add events
    for (event in 1:m) {
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
        tdist[,tnunique] <- 0
        tnunique <- tnunique + 1
        tuntimes[tnunique] <- time[event,i]
        j <- tnunique

        while(j > jstop) {
          tdist[,j] <- tdist[,(j-1)]
          if (j > jstop + 1) {
            # swap indices j and (j-1)
            temp <- tuntimes[j]
            tuntimes[j] <- tuntimes[j-1]
            tuntimes[j-1] <- temp
          }
          j <- j - 1
        }
      }
    }

    # Set jmax
    jmax <- 0
    for (j in 1:tnunique) {
      if (tuntimes[j] < time[1,i]) {
        jmax <- jmax + 1
      }
    }
    if (markov_ind == FALSE) {
      jmax <- min(maxfollow,jmax)
    }
    if (delta[1,i] == 1 | jmax >= tnunique) {
      jmax <- tnunique - 1
    }

     # if (i==182) {
     #   cat('----------------------------------------------','\n')
     #   cat('----------------------------------------------','\n')
     #   cat('i=',i,'\n')
     #   cat('jmax=',jmax,'\n')
     #   cat('time_restriction=',time_restriction,'\n')
     #   cat('tnunique=',tnunique,'\n')
     #   cat('tuntimes=',tuntimes,'\n')
     # }

    if (jmax != 0) {
      for (j in 1:jmax) {
        if (tuntimes[j]>=time_restriction) {break}
        if (tuntimes[j+1]>time_restriction) {
          time_inc=time_restriction-tuntimes[j]
        } else {
          time_inc=tuntimes[j+1]-tuntimes[j]
        }
        # Set state
        state <- 0
        for (current_state in m:1) {
          temp_time_index <- m - current_state + 1
          if (time[temp_time_index,i] <= tuntimes[j] && delta[temp_time_index,i] == 1) {
            state <- current_state
            break
          }
        }
         # if (i==182) {
         #   cat('----------------------------------------------','\n')
         #   cat('j=',j,'\n')
         #   cat('state=',state,'\n')
         #   cat('tuntimes[j]=',tuntimes[j],'\n')
         #   cat('tuntimes[j+1]=',tuntimes[j+1],'\n')
         #   cat('time_inc=',time_inc,'\n')
         # }

        # Calculate rewtp

        # Calculate wins
        for (state_num in 0:(m-1)) {
          if (state_num == state) {
            # Add probabilities from higher states
            for (k in (state_num+1):m) {
              rewtp[i] <- rewtp[i] + tdist[k+1,j] * (time_inc)
              rewtp_components[k,i] <- rewtp_components[k,i] + tdist[k+1,j] * (time_inc)
            }
            break
          }
        }

        # Calculate Losses
        for (current in 1:m) {
          if (current == state) {
            # Subtract probabilities from lower states
            for (k in 1:current) {
              rewtp[i] <- rewtp[i] - tdist[k,j] * (time_inc)
              rewtp_components[current,i] <- rewtp_components[current,i] - tdist[k,j] * (time_inc)
            }
            break
          }
        }
        # if (i==182) {
        #   cat('rewtp_components[,i]=',rewtp_components[,i],'\n')
        # }
      }
    }
  }
   # cat('----------------------------------------------------','\n')
   # cat('final rewtp_component 1=','\n')
   # print(rewtp_components[1,])
   # cat('----------------------------------------------------','\n')

  # Get treatment estimate and variance for Z statistic
  fit_comp <- vector("list",m)
  if (!is.null(cov)) {
    fite=lm(rewtp~trt+cov)
    for (k in 1:m) {
      outcome <- rewtp_components[k,]
      dim(outcome) <- c(n)
      fit_comp[[k]] <- lm(outcome~trt+cov)
    }
  }
  else {
    fite <- lm(rewtp~trt)
    for (k in 1:m) {
      outcome <- rewtp_components[k,]
      dim(outcome) <- c(n)
      fit_comp[[k]] <- lm(outcome~trt)
    }
  }
  rewtp_time=coef(fite)[2]
  rewtp_time_var=vcov(fite)[2,2]
  z_rewtp <- rewtp_time/sqrt(rewtp_time_var)
  for (k in 1:m) {
    components[k] <- coef(fit_comp[[k]])[2]
    components_var[k] <- vcov(fit_comp[[k]])[2,2]
  }

  return(list(rewtp_time,rewtp_time_var,z_rewtp,components,components_var))
}
