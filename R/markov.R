#' Fit a Markov model
#'
#' This function fits an extended Markov model to calculate the state probabilities for each arm. In the wintime package, the returned state
#' probability distributions are used in all non-pairwise methods. The extended  Markov model is recommended for the Expected win time against
#' reference (EWTR) method and the EWTR-composite max test (MAX) method.
#'
#' @param n0 The number of participants in the control arm.
#' @param n1 The number of participants in the active treatment arm.
#' @param m The number of events in the hierarchy.
#' @param Time A `m x (n0 + n1)` matrix of event times (days). Rows should represent events and columns should represent participants. Event
#' rows should be in increasing order of clinical severity.
#' @param Delta A `m x (n0 + n1)` matrix of event indicators. Rows should represent events and columns should represent participants.
#' Event rows should be in increasing order of clinical severity.
#' @return A list containing: a matrix of control arm state probabilities, a matrix of treatment arm state probabilities,
#' a vector of unique control arm event times (days), a vector of unique treatment arm event times (days), the number of unique control
#' arm event times, the number of unique treatment arm event times, the control arm max follow time (days), the treatment arm
#' max follow time (days).
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
#' # -------------------------
#' # markov Examples
#' # -------------------------
#'
#' z <- markov(n0, n1, m, Time, Delta)
#' print(z)
#'
#' @export

# ---------------------
# Markov model
# ---------------------
markov <- function(n0,n1,m,Time,Delta) {
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

  time <- Time[m:1, ]
  delta <- Delta[m:1, ]

  #  Treat_Group==0
  event_times0=rep(NA,times=3*n0)
  trans <- array(data=0,dim=c(m,m+1,3*n0))
  nevent_times0=0

  # Loop 1
  i <- 1
  while (i <= n0) {
    current_state <- 0
    nevent <- 0
    for (num_event in m:2) {
      if (delta[num_event,i] == 1) {
        # If the current event time is greater than a higher level event, do not include it
        lt <- TRUE
        k <- 2
        while (k <= (num_event-1)) {
          if (time[k,i] <= time[num_event,i]) {
            lt <- FALSE
          }
          k <- k + 1
        }
        if (lt == TRUE) {
          nevent <- nevent + 1
          event_times0[nevent_times0+nevent] <- time[num_event,i]
          trans[current_state+1,(m-num_event)+1,nevent_times0+nevent] <- 1
          current_state <- (m-num_event)+1
        }
      }
    }
    if (delta[1,i] == 1) {
      nevent <- nevent + 1
      event_times0[nevent_times0+nevent] <- time[1,i]
      trans[current_state+1,m,nevent_times0+nevent] <- 1
      current_state <- m
    }
    else {
      nevent <- nevent + 1
      event_times0[nevent_times0+nevent] <- time[1,i]
      trans[current_state+1,m+1,nevent_times0+nevent] <- 1
    }
    nevent_times0 <- nevent_times0+nevent
    i <- i+1
  }


  event_times0=event_times0[1:nevent_times0]
  trans <- trans[ , ,1:nevent_times0]
  perm_event_times0=order(event_times0)
  ordered_event_times0=event_times0[perm_event_times0]
  trans <- trans[ , ,perm_event_times0]
  unique_event_times0=ordered_event_times0
  num_state <- matrix(data=0,nrow=m+1,ncol=nevent_times0)

  # Initialize all control patients into state 0
  for (k in 1:nevent_times0) {
    num_state[1,k] <- n0
  }

  # Initialize array holding unique state transitions
  unique_trans <- array(data=0,dim=c(m,m+1,nevent_times0))

  i = 1
  nunique_event_times0 = 0

  # Loop 2
  while (i <= nevent_times0) {
    nunique_event_times0 = nunique_event_times0 + 1
    if (i > 1) {
      for (k in 1:(m+1)) {
        num_state[k, nunique_event_times0] <- num_state[k, nunique_event_times0 - 1]
      }
    }
    unique_event_times0[nunique_event_times0] = ordered_event_times0[i]

    j = i
    end = 0
    while (j < nevent_times0 & end == 0) {
      if (ordered_event_times0[j + 1] == ordered_event_times0[j]) {
        j = j + 1
      } else {
        end = 1
      }
    }
    while (i <= j) {
      # Subtract transitions out of state 0
      num_state[1, nunique_event_times0] <- num_state[1, nunique_event_times0] - sum(trans[1, , i])
      # Add and subtract transitions between state 0 and m
      for (event_num in 2:m) {
        num_state[event_num,nunique_event_times0] <- num_state[event_num,nunique_event_times0] + sum(trans[ ,event_num-1,i]) - sum(trans[event_num, ,i])
      }
      # Add transitions into state m
      num_state[m+1,nunique_event_times0] <- num_state[m+1,nunique_event_times0] + sum(trans[ ,m,i])

      # Populate matrix with unique state transitions
      unique_trans[ , ,nunique_event_times0] <- unique_trans[ , ,nunique_event_times0] + trans[ , ,i]
      i = i + 1
    }
  }

  unique_event_times0=unique_event_times0[1:nunique_event_times0]
  tnum_state <- matrix(data=0,nrow=(m+1),ncol=nunique_event_times0)
  for (k in 1:(m+1)) {
    tnum_state[k,] <- num_state[k,1:nunique_event_times0]
  }

  num_state <- tnum_state
  dim_unique_trans <- dim(unique_trans)
  tunique_trans <- array(NA, dim=c(dim_unique_trans[1], dim_unique_trans[2], nunique_event_times0))

  for (i in 1:m) {
    for (j in 1:m+1) {
      tunique_trans[i,j,1:nunique_event_times0] <- unique_trans[i,j,1:nunique_event_times0]
    }
  }

  dist_state0 <- matrix(data=0,nrow=m+1,ncol=nunique_event_times0)
  dist_state0[1,1] <- 1
  trans_prob0 <- array(data=-1,dim=c(m,m,nunique_event_times0))

  from_state0 <- 0
  for (k in 1:m) {
    from_state0 <- from_state0 + unique_trans[1,k,1]
  }

  dist_state0[1,1] <- dist_state0[1,1] * (1 - from_state0/n0)
  for (k in 2:(m+1)) {
    dist_state0[k,1] <- dist_state0[k,1] + unique_trans[1,k-1,1] / n0
  }

  for (k in 1:m) {
    trans_prob0[1,k,1] <- unique_trans[1,k,1] / n0
  }

  max_follow0=0
  for (k in 1:m) {
    if (dist_state0[k,1] > 0 && num_state[k,1] == 0) {
      max_follow0 <- unique_event_times0[1]
      break
    }
  }


  # Loop 3
  i <- 2
  while(i <= nunique_event_times0) {
    for (event_num in 1:m) {
      if (num_state[event_num,i-1] > 0) {
        dist_state0[event_num,i] <- dist_state0[event_num,i-1] * (1-((sum(unique_trans[event_num, ,i]))-unique_trans[event_num,m+1,i])/num_state[event_num,i-1])
        for (trans_num in 1:m) {
          trans_prob0[event_num,trans_num,i] <- unique_trans[event_num,trans_num,i]/num_state[event_num,i-1]
        }
      }
      else {
        dist_state0[event_num,i] <- dist_state0[event_num,i-1]
      }
    }

    dist_state0[m+1,i] <- dist_state0[m+1,i-1]

    for (event_num in 1:m) {
      if (num_state[event_num,i-1] > 0) {
        for (trans_num in (event_num+1):(m+1)) {
          dist_state0[trans_num,i] <- dist_state0[trans_num,i] + dist_state0[event_num,i-1] * unique_trans[event_num,trans_num-1,i]/num_state[event_num,i-1]
        }
      }
      if (max_follow0 == 0) {
        if (dist_state0[event_num,i] > 0 && num_state[event_num,i] == 0) {
          max_follow0 <- unique_event_times0[i]
        }
      }
    }
    i <- i + 1
  }


  #  Treat_Group==1
  event_times1=rep(NA,times=3*n1)
  trans <- array(data=0,dim=c(m,m+1,3*n1))
  nevent_times1=0

  # Loop 1
  i <- n0+1
  while (i <= (n0+n1)) {
    current_state <- 0
    nevent <- 0
    for (num_event in m:2) {
      if (delta[num_event,i] == 1) {
        # If the current event time is greater than a higher level event, do not include it
        lt <- TRUE
        k <- 2
        while (k <= (num_event-1)) {
          if (time[k,i] <= time[num_event,i]) {
            lt <- FALSE
          }
          k <- k + 1
        }
        if (lt == TRUE) {
          nevent <- nevent + 1
          event_times1[nevent_times1+nevent] <- time[num_event,i]
          trans[current_state+1,(m-num_event)+1,nevent_times1+nevent] <- 1
          current_state <- (m-num_event)+1
        }
      }
    }
    if (delta[1,i] == 1) {
      nevent <- nevent + 1
      event_times1[nevent_times1+nevent] <- time[1,i]
      trans[current_state+1,m,nevent_times1+nevent] <- 1
      current_state <- m
    }
    else {
      nevent <- nevent + 1
      event_times1[nevent_times1+nevent] <- time[1,i]
      trans[current_state+1,m+1,nevent_times1+nevent] <- 1
    }
    nevent_times1 <- nevent_times1+nevent
    i <- i+1
  }


  event_times1=event_times1[1:nevent_times1]
  trans <- trans[ , ,1:nevent_times1]
  perm_event_times1=order(event_times1)
  ordered_event_times1=event_times1[perm_event_times1]
  trans <- trans[ , ,perm_event_times1]
  unique_event_times1=ordered_event_times1
  num_state <- matrix(data=0,nrow=m+1,ncol=nevent_times1)

  # Initialize all treatment patients into state 0
  for (k in 1:nevent_times1) {
    num_state[1,k] <- n1
  }

  # Initialize array holding unique state transitions
  unique_trans <- array(data=0,dim=c(m,m+1,nevent_times1))

  i = 1
  nunique_event_times1 = 0
  # Loop 2
  while (i <= nevent_times1) {
    nunique_event_times1 = nunique_event_times1 + 1
    if (i > 1) {
      for (k in 1:(m+1)) {
        num_state[k, nunique_event_times1] <- num_state[k, nunique_event_times1 - 1]
      }
    }
    unique_event_times1[nunique_event_times1] = ordered_event_times1[i]

    j = i
    end = 0
    while (j < nevent_times1 & end == 0) {
      if (ordered_event_times1[j + 1] == ordered_event_times1[j]) {
        j = j + 1
      } else {
        end = 1
      }
    }
    # Inner Loop
    while (i <= j) {
      # Subtract transitions out of state 0
      num_state[1, nunique_event_times1] <- num_state[1, nunique_event_times1] - sum(trans[1, , i])
      # Add and subtract transitions between state 0 and m
      for (event_num in 2:m) {
        num_state[event_num,nunique_event_times1] <- num_state[event_num,nunique_event_times1] + sum(trans[ ,event_num-1,i]) - sum(trans[event_num, ,i])
      }
      # Add transitions into state m
      num_state[m+1,nunique_event_times1] <- num_state[m+1,nunique_event_times1] + sum(trans[ ,m,i])

      # Populate array with unique state transitions
      unique_trans[ , ,nunique_event_times1] <- unique_trans[ , ,nunique_event_times1] + trans[ , ,i]
      i = i + 1
    }
  }

  unique_event_times1=unique_event_times1[1:nunique_event_times1]
  tnum_state <- matrix(data=0,nrow=(m+1),ncol=nunique_event_times1)
  for (k in 1:(m+1)) {
    tnum_state[k,] <- num_state[k,1:nunique_event_times1]
  }

  num_state <- tnum_state
  dim_unique_trans <- dim(unique_trans)
  tunique_trans <- array(NA, dim=c(dim_unique_trans[1], dim_unique_trans[2], nunique_event_times1))

  for (i in 1:m) {
    for (j in 1:m+1) {
      tunique_trans[i,j,1:nunique_event_times1] <- unique_trans[i,j,1:nunique_event_times1]
    }
  }

  dist_state1 <- matrix(data=0,nrow=m+1,ncol=nunique_event_times1)
  dist_state1[1,1] <- 1
  trans_prob1 <- array(data=-1,dim=c(m,m,nunique_event_times1))

  from_state0 <- 0
  for (k in 1:m) {
    from_state0 <- from_state0 + unique_trans[1,k,1]
  }

  dist_state1[1,1] <- dist_state1[1,1] * (1 - from_state0/n1)
  for (k in 2:(m+1)) {
    dist_state1[k,1] <- dist_state1[k,1] + unique_trans[1,k-1,1] / n1
  }

  for (k in 1:m) {
    trans_prob1[1,k,1] <- unique_trans[1,k,1] / n1
  }

  max_follow1=0
  for (k in 1:m) {
    if (dist_state1[k,1] > 0 && num_state[k,1] == 0) {
      max_follow1 <- unique_event_times1[1]
      break
    }
  }

  # Loop 3
  i <- 2
  while(i <= nunique_event_times1) {
    for (event_num in 1:m) {
      if (num_state[event_num,i-1] > 0) {
        dist_state1[event_num,i] <- dist_state1[event_num,i-1] * (1-((sum(unique_trans[event_num, ,i]))-unique_trans[event_num,m+1,i])/num_state[event_num,i-1])
        for (trans_num in 1:m) {
          trans_prob1[event_num,trans_num,i] <- unique_trans[event_num,trans_num,i]/num_state[event_num,i-1]
        }
      }
      else {
        dist_state1[event_num,i] <- dist_state1[event_num,i-1]
      }
    }

    dist_state1[m+1,i] <- dist_state1[m+1,i-1]

    for (event_num in 1:m) {
      if (num_state[event_num,i-1] > 0) {
        for (trans_num in (event_num+1):(m+1)) {
          dist_state1[trans_num,i] <- dist_state1[trans_num,i] + dist_state1[event_num,i-1] * unique_trans[event_num,trans_num-1,i]/num_state[event_num,i-1]
        }
      }
      if (max_follow1 == 0) {
        if (dist_state1[event_num,i] > 0 && num_state[event_num,i] == 0) {
          max_follow1 <- unique_event_times1[i]
        }
      }
    }
    i <- i + 1
  }
  return(list(dist_state0 = dist_state0,dist_state1 = dist_state1,unique_event_times0 = unique_event_times0,unique_event_times1 = unique_event_times1,
              nunique_event_times0 = nunique_event_times0,nunique_event_times1 = nunique_event_times1,max_follow0 = max_follow0,max_follow1 = max_follow1))
}
