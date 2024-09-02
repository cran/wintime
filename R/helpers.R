#' Helper functions for package functions
#'
#' Win time difference
#'
#' This function calculates the win time difference integral for a single pair. This function is used in all pairwise win time methods.
#'
#' @param m The number of events in the hierarchy.
#' @param etimes A sorted vector of event times (days) (returned from wintime::setEventTimes()).
#' @param time0 A vector containing the control person's event times (days).
#' @param time1 A vector containing the treatment person's event times (days).
#' @param delta0 A vector containing the control person's event indicators.
#' @param delta1 A vector containing the treatment person's event indicators.
#' @return The win time difference integral.

# --------------------------------
# Helper functions
# --------------------------------
getWintimeIntegral <- function(m,etimes,time0,time1,delta0,delta1) {
  integral <- 0
  # Iterate through event times and calculate integral
  for (k in 2:length(etimes)) {
    if (etimes[k] != etimes[k-1]) {
      for (event in m:1) {
        if (time0[event] <= etimes[k-1] && time1[event] > etimes[k-1]) {
          # Treatment person wins
          integral <- integral + etimes[k] - etimes[k-1]
          break
        }
        else if (time1[event] <= etimes[k-1] && time0[event] > etimes[k-1]) {
          # Treatment person loses
          integral <- integral - etimes[k] + etimes[k-1]
          break
        }
        # If there is a tie
        else if (time1[event] <= etimes[k-1] && time0[event] <= etimes[k-1]) {
          break
        }
      }
    }
  }

  # If an event is observed on the last day, update the integral value to reflect it
  for (event in m:1) {
    # If both the control and the treatment person have an event time equal to the last day
    # and the control person has an event while the treatment person does not
    if (time0[event] == etimes[length(etimes)] && time1[event] == etimes[length(etimes)]
        && delta0[event] == 1 && delta1[event] == 0) {
      integral <- integral + 1
      break
    }
    # If both the control and the treatment person have an event time equal to the last day
    # and the control person does not have an event while the treatment person does
    if (time0[event] == etimes[length(etimes)] && time1[event] == etimes[length(etimes)]
        && delta0[event] == 0 && delta1[event] == 1) {
      integral <- integral - 1
      break
    }
  }
  return(integral)
}

#' Set event times and indicators used in the Kaplan-Meier survival curve calculation
#'
#' This function creates the time_km and delta_km matrices used for wintime::km().
#'
#' @param n The total number of trial participants.
#' @param m The number of events in the hierarchy.
#' @param time The row reversal of the Time matrix (days) (created inside wintime::km()).
#' @param delta The row reversal of the Delta matrix (created inside wintime::km()).
#' @return A list containing the event time matrix and the event indicator matrix used in wintime::km().

setKM <- function(n,m,time,delta) {
  # Initialize matrices to hold KM times and indicators
  time_km <- matrix(data=0,nrow=m,ncol=n)
  delta_km <- matrix(data=0,nrow=m,ncol=n)


  time_km[1, ] <- time[1, ]
  delta_km[1, ] <- delta[1, ]

  for (event_num in 2:m) {
    time_km[event_num, ] <- pmin(time_km[event_num-1, ], time[event_num, ])
    delta_km[event_num, ] <- delta[event_num, ]
  }

  # Adjust delta_km matrix
  i <- 1
  while (i <= n) {
    for (current_event in 2:m) {
      for (observed_event in 1:(current_event-1)) {
        if (time[observed_event,i] == time_km[current_event,i]) {
          delta_km[current_event,i] <- delta[observed_event,i]
          break
        }
      }
    }
    i <- i + 1
  }
  return(list(time_km,delta_km))
}

#' Created a sorted vector of event times
#'
#' This function creates a sorted vector of event times for a pair. This function is used in all pairwise functions.
#'
#' @param m The number of events in the hierarchy.
#' @param delta0 A vector of event indicators for the control person.
#' @param delta1 A vector of event indicators for the treatment person.
#' @param time0 A vector of event times (days) for the control person.
#' @param time1 A vector of event times (days) for the treatment person.
#' @param follow The maximum follow up time (days).
#' @return A sorted vector of event times (days) for a given pair.

setEventTimes <- function(m,delta0,delta1,time0,time1,follow) {
  # Initialize vector that will hold event times for both groups
  etimes <- numeric(0)

  # Control arm
  for (event in m:1) {
    # If the current event is not observed in person j,
    # or their current event time is greater than the max follow-up time
    if (delta0[event] == 0 || time0[event] > follow) {
      # Set their current event time to the max follow-up time and their indicator to 0
      time0[event] <- follow
      delta0[event] <- 0
    }

    # Append current event time to the event times vector
    etimes <- c(etimes,time0[event])
  }

  #Treatment arm
  for (event in m:1) {
    # If the current event is not observed in person i,
    # or their current event time is greater than the max follow-up time
    if (delta1[event] == 0 || time1[event] > follow) {
      # Set their current event time to the follow-up time and their indicator to 0
      time1[event] <- follow
      delta1[event] <- 0
    }
    # Append current event time to the event times vector
    etimes <- c(etimes,time1[event])
  }

  etimes <- c(etimes,follow)

  # Sorted array of event times
  etimes <- sort(etimes)
  return(etimes)
}
