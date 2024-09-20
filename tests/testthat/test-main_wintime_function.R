library(testthat)
library(survival)
library(wintime)

# ------------------------------------------------
# Example matrices and vectors for testing
# ------------------------------------------------

# Event time vectors
TIME_1 <- c(256,44,29,186,29,190,292,26,11,858,47,33,11,230,394,58,179,303,186,110)
TIME_2 <- c(128,44,95,186,69,190,332,7,11,497,47,33,234,321,219,58,29,303,172,427)
TIME_3 <- c(435,44,95,186,69,190,332,325,11,1348,47,33,234,1085,442,58,219,303,186,544)

# Event time matrix
Time <- rbind(TIME_1, TIME_2, TIME_3)

# Event indicator vectors
DELTA_1 <- c(1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1)
DELTA_2 <- c(1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1)
DELTA_3 <- c(0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0)

# Event indicator matrix
Delta <- rbind(DELTA_1, DELTA_2, DELTA_3)

# Treatment arm indicator vector
trt <- c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0)

# Covariate vectors
cov1 <- c(66,67,54,68,77,65,55,66,77,54,23,13,87,105,32,5,22,19,66,92)
cov2 <- c(3,6,4,2,3,5,8,5,3,5,1,1,9,14,0,3,7,2,1,5)
cov3 <- c(34.6,543.6,45.8,54.7,44.3,55.6,65.9,54.7,77.9,31.2,11,55.8,79.0,323.5,54.2,28.6,452.1,77.9,41.6,89.3)

# Covariate matrix
cov <- cbind(cov1, cov2, cov3)

# ---------------------------------
# Test wintime::wintime()
# ---------------------------------

test_that("wintime function works with valid inputs", {
  # Test 'ewt' method with default settings
  result <- wintime("ewt", Time, Delta, trt)
  expect_type(result, "list")
  expect_named(result, c("data","resample_data","message","variance","p","wins","losses"))

  # Test 'ewtr' method with default settings
  result <- wintime("ewtr", Time, Delta, trt, cov = cov)
  expect_type(result, "list")
  expect_named(result, c("data","resample_data","message","variance","p","wins","losses"))

  # Test 'rmt' method with default settings and 10 resamples
  result <- wintime("rmt", Time, Delta, trt, rmst_restriction = 365, resample_num = 10)
  expect_type(result, "list")
  expect_named(result, c("data","resample_data","message","variance","p","wins","losses"))

  # Test 'max' method with default settings and 10 resamples
  result <- wintime("max", Time, Delta, trt, resample_num = 10)
  expect_type(result, "list")
  expect_named(result, c("data","resample_data","message","variance","p","wins","losses"))

  # Test 'wtr' method with default settings
  result <- wintime("wtr", Time, Delta, trt)
  expect_type(result, "list")
  expect_named(result, c("data","resample_data","message","variance","p","wins","losses"))

  # Test 'rwtr' method with default settings and 5 resamples
  result <- wintime("rwtr", Time, Delta, trt, resample_num = 5)
  expect_type(result, "list")
  expect_named(result, c("data","resample_data","message","variance","p","wins","losses"))

  # Test 'pwt' method with default settings
  result <- wintime("pwt", Time, Delta, trt)
  expect_type(result, "list")
  expect_named(result, c("data","resample_data","message","variance","p","wins","losses"))
})

test_that("wintime function handles invalid inputs", {
  # Invalid type
  expect_error(wintime("invalid_type", Time, Delta, trt), "Invalid type:")

  # Non-numeric input for Time
  expect_error(wintime("ewtr", "invalid", Delta, trt), "Input vectors and matrices must be numeric.")

  # Non-matrix input for Delta
  expect_error(wintime("ewtr", Time, Delta[1, ], trt), "Event times and event indicators must be organized into matrices. Treatment arm indicators must be contained in a single vector.")

  # Mismatched dimensions
  TIME_1_mismatch <- c(256,44,29,186,29,80,11,380,102)
  TIME_2_mismatch <- c(128,44,95,186,69,66,153,380,117)
  TIME_3_mismatch <- c(435,44,95,186,69,270,1063,380,117)

  # Event time matrix
  Time_mismatch <- rbind(TIME_1_mismatch, TIME_2_mismatch, TIME_3_mismatch)
  expect_error(wintime("ewtr", Time_mismatch, Delta, trt), "Input matrices and vectors must have compatible dimensions and lengths.")
})

test_that("wintime function handles edge cases", {
  # Test with empty inputs
  Time_empty <- matrix(numeric(0), nrow = 0, ncol = 0)
  Delta_empty <- matrix(numeric(0), nrow = 0, ncol = 0)
  trt_empty <- numeric(0)
  expect_error(wintime("ewtr", Time_empty, Delta_empty, trt_empty), "Input matrices and vectors cannot be empty.")
})

test_that("wintime function handles warnings", {
  # Test warnings for 'ewt' with Markov model
  expect_warning(result <- wintime("ewt", Time, Delta, trt, model = "markov"), "For this method, it is strongly recommended to use a KM model and to resample using permutations. These are set as defaults.")

  # Test warnings for 'ewtr' with KM model
  expect_warning(result <- wintime("ewtr", Time, Delta, trt, model = "km"), "For this method, it is strongly recommended to use a Markov model. This is set as a default.")
})
