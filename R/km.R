#' Fit a Kaplan-Meier model
#'
#' This function fits Kaplan-Meier models to calculate the state probabilities for each arm. In the wintime package, the returned state probability
#' distributions are used in all non-pairwise methods. The Kaplan-Meier model is recommended for the Expected win time (EWT) method and the Restricted mean survival
#' in favor of treatment (RMT) method.
#'
#' @param n0 The number of participants in the control arm.
#' @param n1 The number of participants in the treatment arm.
#' @param m The number of events in the hierarchy.
#' @param Time A `m x (n0 + n1)` matrix of event times (days). Rows should represent events and columns should represent participants. Event rows should be
#' in increasing order of clinical severity.
#' @param Delta A `m x (n0 + n1)` matrix of event indicators. Rows should represent events and columns should represent participants. Event rows should be
#' in increasing order of clinical severity.
#' @return A list containing: a matrix of control arm state probabilities, a matrix of treatment arm state probabilities,
#' a vector of unique control arm event times (days), a vector of unique treatment arm event times (days), the number of unique control
#' arm event times, the number of unique treatment arm event times, the control arm max follow time (days), the treatment arm
#' max follow time (days), a matrix of combined arm state probabilities, a vector of unique combined arm event times (days),
#' the number of unique combined arm event times, the combined arm max follow time (days), a (m x number unique combined arm event times)
#' matrix of combined arm km survival probabilities, matrix of trt arm km survival probabilities, matrix of control arm km survival probabilities.
#' @examples
#' # -----------------------------
#' # Example inputs
#' # -----------------------------
#'
#' # Event time vectors
#' TIME_1 <- c(256,44,29,186,29,80,11,380,102,33)
#' TIME_2 <- c(128,44,95,186,69,66,153,380,117,33)
#' TIME_3 <- c(435,44,95,186,69,270,1063,380,117,33)
#'
#' # Event time matrix
#' Time <- rbind(TIME_1, TIME_2, TIME_3)
#'
#' # Event indicator vectors
#' DELTA_1 <- c(1,0,1,0,1,1,1,0,1,0)
#' DELTA_2 <- c(1,0,0,0,0,1,1,0,0,0)
#' DELTA_3 <- c(0,0,0,0,0,0,0,0,0,0)
#'
#' # Event indicator matrix
#' Delta <- rbind(DELTA_1, DELTA_2, DELTA_3)
#'
#' # Treatment arm indicator vector
#' trt <- c(1,1,1,1,1,0,0,0,0,0)
#'
#' # Number of control arm patients
#' n0 <- sum(trt == 0)
#'
#' # Number of treatment arm patients
#' n1 <- sum(trt == 1)
#'
#' # Number of events in the hierarchy
#' m <- nrow(Time)
#'
#' # ------------------------
#' # km Examples
#' # ------------------------
#'
#' z <- km(n0, n1, m, Time, Delta)
#' print(z)
#'
#' @import survival
#' @export

# -------------------------
# Kaplan-Meier model
# -------------------------
km <- function(n0,n1,m,Time,Delta) {
  # -------------------------------
  # Validate inputs
  # -------------------------------
  if (!is.numeric(Time) || !is.numeric(Delta) || !is.numeric(n0) || !is.numeric(n1) || !is.numeric(m)) {
    stop("Input data must be numeric.")
  }
  if (length(Time) == 0 || length(Delta) == 0) {
    stop("Input matrices cannot be empty.")
  }
  if (nrow(Time) != nrow(Delta) || ncol(Time) != ncol(Delta)) {
    stop("Input matrices must have compatible dimensions.")
  }
  if (nrow(Time) != m || nrow(Delta) != m || ncol(Time) != (n0+n1) || ncol(Delta) != (n0+n1)) {
    stop("Rows of input matrices must correspond to events. Columns must correspond to participants. Ensure that the number of events
         is equal to the number of rows and that the number of participants is equal to the columns.")
  }

  n <- n0 + n1
  time <- Time[m:1, ]
  delta <- Delta[m:1, ]

  # Set KM time/delta matrices
  result <- setKM(n,m,time,delta)
  time_km <- result[[1]]
  delta_km <- result[[2]]

  #---------------------------------------------------------------------
  # CONTROL GROUP
  #
  #---------------------------------------------------------------------

  # Create KM survival curve matrices
  Time0 <- matrix(0,nrow=m,ncol=n)
  Surv0 <- matrix(0,nrow=m,ncol=n)
  conkm <- matrix(0,nrow=m,ncol=m*n)

  for (k in 1:m) {
    # temp <- survfit(Surv(time_km[k,dataset1$trt==0],delta_km[k,dataset1$trt==0])~1)$time
    # temp2 <- survfit(Surv(time_km[k,dataset1$trt==0],delta_km[k,dataset1$trt==0])~1)$surv
    temp <- survival::survfit(survival::Surv(time_km[k,1:n0],delta_km[k,1:n0])~1)$time
    temp2 <- survival::survfit(survival::Surv(time_km[k,1:n0],delta_km[k,1:n0])~1)$surv
    Time0[k,1:length(temp)] <- temp
    Surv0[k,1:length(temp2)] <- temp2
  }

  nkm <- numeric(m)
  for (k in 1:m) {
    nkm[k] <- sum(Time0[k,] != 0)
  }

  dist_state0 <- matrix(0,nrow=(m+1),ncol=sum(nkm))
  unique_event_times0 <- rep(0,times=sum(nkm))

  # Get max follow time
  control_follow_times <- rep(0,times=m)
  for (i in 1:m) {
    control_follow_times[i] <- Time0[i,nkm[i]]
  }
  max_follow0 <- min(control_follow_times)

  time <- 0
  count <- rep(0,m)
  max_times <- rep(0,m)
  nunique_event_times0 <- 0

  next_count <- rep(0,m)

  while (time < max_follow0) {
    for (i in 1:m) {
      next_count[i] <- count[i] + 1
    }
    for (i in 1:m) {
      if (next_count[i] > nkm[i]) {
        next_count[i] <- nkm[i]
        max_times[i] <- 1
      }
    }

    time <- max(Time0)
    for (k in 1:m) {
      if (max_times[k] == 0) {
        if (Time0[k,next_count[k]] < time) {
          time <- Time0[k,next_count[k]]
        }
      }
    }

    nunique_event_times0=nunique_event_times0+1
    unique_event_times0[nunique_event_times0]=time

    #  Update conkm[]
    #
    for (k in 1:m) {
      if (Time0[k,next_count[k]] == time) {
        conkm[k,nunique_event_times0] <- Surv0[k,next_count[k]]
      }
      else {
        if (nunique_event_times0==1) {
          conkm[k,nunique_event_times0] <- 1.0
        }
        else {
          conkm[k,nunique_event_times0] <- conkm[k,nunique_event_times0-1]
        }
      }
    }

    for (k in 1:m) {
      if (Time0[k,next_count[k]] <= time) {
        count[k] <- next_count[k]
      }
    }

    # Set probability of being in state 0
    if (count[m] > 0) {
      dist_state0[1,nunique_event_times0] <- Surv0[m,count[m]]
    }
    else {
      dist_state0[1,nunique_event_times0] <- 1
    }

    # Set probability of being in current state
    for (i in 1:m) {
      state_num <- (m - i) + 1
      if (count[i] > 0) {
        dist_state0[state_num+1,nunique_event_times0] <- 1 - Surv0[i,count[i]]
        # Subtract previously calculated probabilities
        if (state_num < m) {
          for (j in (state_num+1):m) {
            dist_state0[(state_num+1),nunique_event_times0] <- dist_state0[(state_num+1),nunique_event_times0] - dist_state0[(j+1),nunique_event_times0]
          }
        }
      }
    }
  }
  unique_event_times0=unique_event_times0[1:nunique_event_times0]
  conkm=conkm[1:m,1:nunique_event_times0]

  #---------------------------------------------------------------------
  # TREATMENT GROUP
  #
  #---------------------------------------------------------------------

  # Create KM survival curve matrices
  Time1 <- matrix(0,nrow=m,ncol=n)
  Surv1 <- matrix(0,nrow=m,ncol=n)
  trtkm <- matrix(0,nrow=m,ncol=m*n)

  for (k in 1:m) {
    # temp <- survfit(Surv(time_km[k,dataset1$trt==1],delta_km[k,dataset1$trt==1])~1)$time
    # temp2 <- survfit(Surv(time_km[k,dataset1$trt==1],delta_km[k,dataset1$trt==1])~1)$surv
    temp <- survival::survfit(survival::Surv(time_km[k,(n0+1):n],delta_km[k,(n0+1):n])~1)$time
    temp2 <- survival::survfit(survival::Surv(time_km[k,(n0+1):n],delta_km[k,(n0+1):n])~1)$surv
    Time1[k,1:length(temp)] <- temp
    Surv1[k,1:length(temp2)] <- temp2
  }

  nkm <- numeric(m)
  for (k in 1:m) {
    nkm[k] <- sum(Time1[k,] != 0)
  }

  dist_state1 <- matrix(0,nrow=(m+1),ncol=sum(nkm))
  unique_event_times1 <- rep(0,times=sum(nkm))

  # Get max follow time
  trt_follow_times <- rep(0,times=m)
  for (i in 1:m) {
    trt_follow_times[i] <- Time1[i,nkm[i]]
  }

  max_follow1 <- min(trt_follow_times)

  time <- 0
  count <- rep(0,m)
  max_times <- rep(0,m)
  nunique_event_times1 <- 0

  next_count <- rep(0,m)

  while (time < max_follow1) {
    for (i in 1:m) {
      next_count[i] <- count[i] + 1
    }
    for (i in 1:m) {
      if (next_count[i] > nkm[i]) {
        next_count[i] <- nkm[i]
        max_times[i] <- 1
      }
    }

    time <- max(Time1)
    for (k in 1:m) {
      if (max_times[k] == 0) {
        if (Time1[k,next_count[k]] < time) {
          time <- Time1[k,next_count[k]]
        }
      }
    }

    nunique_event_times1=nunique_event_times1+1
    unique_event_times1[nunique_event_times1]=time

    #  Update trtkm[]
    #
    for (k in 1:m) {
      if (Time1[k,next_count[k]] == time) {
        trtkm[k,nunique_event_times1] <- Surv1[k,next_count[k]]
      }
      else {
        if (nunique_event_times1==1) {
          trtkm[k,nunique_event_times1] <- 1.0
        }
        else {
          trtkm[k,nunique_event_times1] <- trtkm[k,nunique_event_times1-1]
        }
      }
    }

    for (k in 1:m) {
      if (Time1[k,next_count[k]] <= time) {
        count[k] <- next_count[k]
      }
    }

    # Set probability of being in state 0
    if (count[m] > 0) {
      dist_state1[1,nunique_event_times1] <- Surv1[m,count[m]]
    }
    else {
      dist_state1[1,nunique_event_times1] <- 1
    }

    # Set probability of being in current state
    for (i in 1:m) {
      state_num <- (m - i) + 1
      if (count[i] > 0) {
        dist_state1[state_num+1,nunique_event_times1] <- 1 - Surv1[i,count[i]]
        # Subtract previously calculated probabilities
        if (state_num < m) {
          for (j in (state_num+1):m) {
            dist_state1[(state_num+1),nunique_event_times1] <- dist_state1[(state_num+1),nunique_event_times1] - dist_state1[(j+1),nunique_event_times1]
          }
        }
      }
    }
  }
  unique_event_times1=unique_event_times1[1:nunique_event_times1]
  trtkm=trtkm[1:m,1:nunique_event_times1]

  #---------------------------------------------------------------------
  # COMBINED ARM
  #
  #---------------------------------------------------------------------

  # Create KM survival curve matrices
  Time2 <- matrix(0,nrow=m,ncol=n)
  Surv2 <- matrix(0,nrow=m,ncol=n)
  comkm <- matrix(0,nrow=m,ncol=m*n)

  for (k in 1:m) {
    # temp <- survfit(Surv(time_km[k,dataset1$trt==0],delta_km[k,dataset1$trt==0])~1)$time
    # temp2 <- survfit(Surv(time_km[k,dataset1$trt==0],delta_km[k,dataset1$trt==0])~1)$surv
    temp <- survival::survfit(survival::Surv(time_km[k,1:n],delta_km[k,1:n])~1)$time
    temp2 <- survival::survfit(survival::Surv(time_km[k,1:n],delta_km[k,1:n])~1)$surv
    Time2[k,1:length(temp)] <- temp
    Surv2[k,1:length(temp2)] <- temp2
  }
#  cat('length(Time2[1,])=',length(Time2[1,]),'\n')
#  cat('length(Time2[2,])=',length(Time2[2,]),'\n')
#  cat('length(Time2[3,])=',length(Time2[3,]),'\n')
#  cat('Time2[1,1:120]=',Time2[1,1:120],'\n')
#  cat('Surv2[1,1:15]=',Surv2[1,1:15],'\n')
#  cat('Time2[2,1:120]=',Time2[2,1:120],'\n')
#  cat('Surv2[2,1:15]=',Surv2[2,1:15],'\n')
#  cat('Time2[3,1:120]=',Time2[3,1:120],'\n')
#  cat('Surv2[3,1:15]=',Surv2[3,1:15],'\n')
#  cat('------------------------------------------','\n')

  nkm <- numeric(m)
  for (k in 1:m) {
    nkm[k] <- sum(Time2[k,] != 0)
  }

#  cat('nkm[1]=',nkm[1],'\n')
#  cat('nkm[2]=',nkm[2],'\n')
#  cat('nkm[3]=',nkm[3],'\n')
#  cat('------------------------------------------','\n')


  dist_state2 <- matrix(0,nrow=(m+1),ncol=sum(nkm))
  unique_event_times2 <- rep(0,times=sum(nkm))

  # Get max follow time
  com_follow_times <- rep(0,times=m)
  for (i in 1:m) {
    com_follow_times[i] <- Time2[i,nkm[i]]
  }
  max_follow2 <- min(com_follow_times)

#   cat('max_follow2=',max_follow2,'\n')
#   cat('------------------------------------------','\n')


  time <- 0
  count <- rep(0,m)
  max_times <- rep(0,m)
  nunique_event_times2 <- 0

  next_count <- rep(0,m)

  while (time < max_follow2) {
    for (i in 1:m) {
      next_count[i] <- count[i] + 1
    }
    for (i in 1:m) {
      if (next_count[i] > nkm[i]) {
        next_count[i] <- nkm[i]
        max_times[i] <- 1
      }
    }

    time <- max(Time2)
    for (k in 1:m) {
      if (max_times[k] == 0) {
        if (Time2[k,next_count[k]] < time) {
          time <- Time2[k,next_count[k]]
        }
      }
    }

    nunique_event_times2=nunique_event_times2+1
    unique_event_times2[nunique_event_times2]=time

    #  Update comkm[]
    #
    for (k in 1:m) {
      if (Time2[k,next_count[k]] == time) {
        comkm[k,nunique_event_times2] <- Surv2[k,next_count[k]]
      }
      else {
        if (nunique_event_times2==1) {
          comkm[k,nunique_event_times2] <- 1.0
        }
        else {
          comkm[k,nunique_event_times2] <- comkm[k,nunique_event_times2-1]
        }
      }
    }

    # if (time < Time2[3,15]) {
    # cat('time=',time,'\n')
    # cat('nunique_event_times2=',nunique_event_times2,'\n')
    # cat('next_count[1]=',next_count[1],'\n')
    # cat('next_count[2]=',next_count[2],'\n')
    # cat('next_count[3]=',next_count[3],'\n')
    # cat('comkm[1,1:15]=',comkm[1,1:15],'\n')
    # cat('comkm[2,1:15]=',comkm[2,1:15],'\n')
    # cat('comkm[3,1:15]=',comkm[3,1:15],'\n')
    # cat('------------------------------------------','\n')
    # }

    for (k in 1:m) {
      if (Time2[k,next_count[k]] <= time) {
        count[k] <- next_count[k]
      }
    }

    # Set probability of being in state 0
    if (count[m] > 0) {
      dist_state2[1,nunique_event_times2] <- Surv2[m,count[m]]
    }
    else {
      dist_state2[1,nunique_event_times2] <- 1
    }

    # Set probability of being in current state
    for (i in 1:m) {
      state_num <- (m - i) + 1
      if (count[i] > 0) {
        dist_state2[state_num+1,nunique_event_times2] <- 1 - Surv2[i,count[i]]
        # Subtract previously calculated probabilities
        if (state_num < m) {
          for (j in (state_num+1):m) {
            dist_state2[(state_num+1),nunique_event_times2] <- dist_state2[(state_num+1),nunique_event_times2] - dist_state2[(j+1),nunique_event_times2]
          }
        }
      }
    }
  }
  unique_event_times2=unique_event_times2[1:nunique_event_times2]
  dist_state2=dist_state2[1:(m+1),1:nunique_event_times2]
  comkm=comkm[1:m,1:nunique_event_times2]
#  cat('----------------------------------------------------','\n')
#  cat('nunique_event_times2=',nunique_event_times2,'\n')
#  cat('unique_event_times2=',unique_event_times2,'\n')
#  cat('dist_state2=','\n')
#  print(dist_state2)
#  cat('----------------------------------------------------','\n')

  return(list(dist_state0 = dist_state0,dist_state1 = dist_state1,unique_event_times0 = unique_event_times0,unique_event_times1 = unique_event_times1,
              nunique_event_times0 = nunique_event_times0,nunique_event_times1 = nunique_event_times1,max_follow0 = max_follow0,max_follow1 = max_follow1,
              dist_state2 = dist_state2, unique_event_times2 = unique_event_times2, nunique_event_times2 = nunique_event_times2, max_follow2 = max_follow2,
              comkm=comkm,trtkm=trtkm,conkm=conkm))
}
