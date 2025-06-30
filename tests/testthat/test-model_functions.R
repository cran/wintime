library(testthat)
library(survival)
library(wintime)

# ----------------------------------
# Example inputs
# ----------------------------------

# Event time vectors
Time1 <- c(256,44,29,186,29,190,80,11,380,102,33,303)
Time2 <- c(128,44,95,186,69,190,66,153,380,117,33,131)

# Event time matrix
Time <- rbind(Time1,Time2)

# Event indicator vectors
Delta1 <- c(1,0,1,0,1,0,1,1,0,1,0,1)
Delta2 <- c(1,0,0,0,0,0,1,1,0,0,0,1)

# Event indicator matrix
Delta <- rbind(Delta1,Delta2)

# Treatment arm indicator vector
trt <- c(1,1,1,1,1,1,0,0,0,0,0,0)

# Number of control arm patients
n0 <- sum(trt == 0)

# Number of treatment arm patients
n1 <- sum(trt == 1)

# Number of events in the hierarchy
m <- nrow(Time)

# ---------------------------------
# Test wintime::markov()
# ---------------------------------
test_that("markov function handles valid inputs", {
  result <- markov(n0, n1, m, Time, Delta)
  expect_type(result, "list")
  expect_named(result, c("dist_state0", "dist_state1", "unique_event_times0", "unique_event_times1", "nunique_event_times0", "nunique_event_times1",
                         "max_follow0", "max_follow1","dist_state2","unique_event_times2","nunique_event_times2","max_follow2","trans_prob2",
                         "trans_prob1","trans_prob0"))
})

test_that("markov function handles invalid inputs", {
  # Non-numeric inputs
  expect_error(markov(n0, n1, "invalid", Time, Delta), "Input data must be numeric.")

  # Mismatched dimensions
  Delta_mismatch1 <- c(1,0,1,0,1,0,1,0,1,0,1)
  Delta_mismatch2 <- c(1,0,1,0,1,0,1,0,1,0,1)
  Delta_mismatch <- rbind(Delta_mismatch1, Delta_mismatch2)

  expect_error(markov(n0, n1, m, Time, Delta_mismatch), "Input matrices must have compatible dimensions.")

  # Incorrectly formatted matrices
  trans_Time <- t(Time)
  trans_Delta <- t(Delta)

  expect_error(markov(n0, n1, m, trans_Time, trans_Delta), "Rows of input matrices must correspond to events. Columns must correspond to participants. Ensure that the number of events
         is equal to the number of rows and that the number of participants is equal to the columns.")
})

test_that("markov function handles edge cases", {
  # Empty input matrices
  Time_empty <- matrix(numeric(0), nrow = 0, ncol = 0)
  Delta_empty <- matrix(numeric(0), nrow = 0, ncol = 0)

  expect_error(markov(n0, n1, m, Time_empty, Delta_empty), "Input matrices cannot be empty.")
})

# ------------------------------
# Test wintime::km()
# ------------------------------

test_that("km function handles valid inputs", {
  result <- km(n0, n1, m, Time, Delta)
  expect_type(result, "list")
  expect_named(result, c("dist_state0", "dist_state1", "unique_event_times0", "unique_event_times1", "nunique_event_times0", "nunique_event_times1",
                         "max_follow0", "max_follow1", "dist_state2","unique_event_times2","nunique_event_times2","max_follow2","comkm","trtkm","conkm"))
})

test_that("km function handles invalid inputs", {
  # Non-numeric inputs
  expect_error(km(n0, n1, "invalid", Time, Delta), "Input data must be numeric.")

  # Mismatched dimensions
  Delta_mismatch1 <- c(1,0,1,0,1,0,1,0,1,0,1)
  Delta_mismatch2 <- c(1,0,1,0,1,0,1,0,1,0,1)
  Delta_mismatch <- rbind(Delta_mismatch1, Delta_mismatch2)

  expect_error(km(n0, n1, m, Time, Delta_mismatch), "Input matrices must have compatible dimensions.")

  # Incorrectly formatted matrices
  trans_Time <- t(Time)
  trans_Delta <- t(Delta)

  expect_error(km(n0, n1, m, trans_Time, trans_Delta), "Rows of input matrices must correspond to events. Columns must correspond to participants. Ensure that the number of events
         is equal to the number of rows and that the number of participants is equal to the columns.")
})

test_that("km function handles edge cases", {
  # Empty input matrices
  Time_empty <- matrix(numeric(0), nrow = 0, ncol = 0)
  Delta_empty <- matrix(numeric(0), nrow = 0, ncol = 0)

  expect_error(km(n0, n1, m, Time_empty, Delta_empty), "Input matrices cannot be empty.")
})
