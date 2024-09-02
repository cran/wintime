## ----- Load Packages----------------------------------------------------------
library(wintime)
library(survival)

## ----- Create Data------------------------------------------------------------

# Event time vectors
TIME_1 <- c(256,44,29,186,29,80,11,380,102,33)
TIME_2 <- c(128,44,95,186,69,66,153,380,117,33)
TIME_3 <- c(435,44,95,186,69,270,1063,380,117,33)

# Event time matrix
Time <- rbind(TIME_1, TIME_2, TIME_3)

# Event indicator vectors
DELTA_1 <- c(1,0,1,0,1,1,1,0,1,0)
DELTA_2 <- c(1,0,0,0,0,1,1,0,0,0)
DELTA_3 <- c(0,0,0,0,0,0,0,0,0,0)

# Event indicator matrix
Delta <- rbind(DELTA_1, DELTA_2, DELTA_3)

# Treatment arm indicator vector
trt <- c(1,1,1,1,1,0,0,0,0,0)

# Covariate vectors
cov1 <- c(66,67,54,68,77,65,55,66,77,54)
cov2 <- c(3,6,4,2,3,5,8,5,3,5)
cov3 <- c(34.6,543.6,45.8,54.7,44.3,55.6,65.9,54.7,77.9,31.2)

# Covariate matrix
cov <- cbind(cov1, cov2, cov3)

cat("Time =", "\n")
print(Time)
cat("Delta =", "\n")
print(Delta)
cat("trt", "\n")
print(trt)
cat("cov =", "\n")
print(cov)

## ----- wtr--------------------------------------------------------------------
# Run wtr
result <- wintime("wtr", Time, Delta, trt)
print(result)

## ----- rwtr-------------------------------------------------------------------
# Run rwtr
result <- wintime("rwtr", Time, Delta, trt)
print(result)

## ----- pwt--------------------------------------------------------------------
# Run pwt
result <- wintime("pwt", Time, Delta, trt)
print(result)

## ----- ewtr-------------------------------------------------------------------
# Run ewtr with covariates 
result <- wintime("ewtr", Time, Delta, trt, cov = cov)
print(result)

## ----- ewt--------------------------------------------------------------------
# Run ewt
result <- wintime("ewt", Time, Delta, trt)
print(result)

## ----- max--------------------------------------------------------------------
# Run max
result <- wintime("max", Time, Delta, trt, cov = cov)
print(result)

## ----- rmt--------------------------------------------------------------------
# Run rmt
result <- wintime("rmt", Time, Delta, trt, rmst_restriction = round(1.5*365.25))
print(result)

## ----- ewt resampling---------------------------------------------------------
# Run ewt with 10 resamples
result <- wintime("ewt", Time, Delta, trt, resample_num = 10, seed = 123)
print(result)

## ----- wtr resampling---------------------------------------------------------
# Run wtr with 5 resamples
result <- wintime("wtr", Time, Delta, trt, resample_num = 5)
print(result)

## ----- ewt bootstraps---------------------------------------------------------
# Run ewt on 10 bootstraps
result <- wintime("ewt", Time, Delta, trt, resample_num = 10, resample = "boot")
print(result)

## ----- ewtr w/ KM model-------------------------------------------------------
# Run ewtr with a KM model
result <- wintime("ewtr", Time, Delta, trt, cov = cov, model = "km")
print(result)

## ----- set parameters---------------------------------------------------------
# Number of control arm patients
n0 <- sum(trt == 0)

# Number of treatment arm patients
n1 <- sum(trt == 1)

# Number of endpoints
m <- nrow(Time)

## ----- markov-----------------------------------------------------------------
# Run markov
result <- markov(n0, n1, m, Time, Delta)

# Control arm probabilities
dist0 <- result[[1]]
cat("dist0 =", "\n")
print(dist0)

# Treatment arm probabilities
dist1 <- result[[2]]
cat("dist1 =", "\n")
print(dist1)

# Unique control arm event times
untimes0 <- result[[3]]
cat("untimes0 =", "\n")
print(untimes0)

# Unique treatment arm event times
untimes1 <- result[[4]]
cat("untimes1 =", "\n")
print(untimes1)

# Number of unique control arm event times
nuntimes0 <- result[[5]]
cat("nuntimes0 =", "\n")
print(nuntimes0)

# Number of unique treatment arm event times
nuntimes1 <- result[[6]]
cat("nuntimes1 =", "\n")
print(nuntimes1)

# Control arm max follow time
max_follow0 <- result[[7]]
cat("max_follow0 =", "\n")
print(max_follow0)

# Treatment arm max follow time
max_follow1 <- result[[8]]
cat("max_follow1 =", "\n")
print(max_follow1)

## ----- km---------------------------------------------------------------------
# Run km
result <- km(n0, n1, m, Time, Delta)

# Control arm probabilities
dist0 <- result[[1]]
cat("dist0 =", "\n")
print(dist0)

# Treatment arm probabilities
dist1 <- result[[2]]
cat("dist1 =", "\n")
print(dist1)

# Unique control arm event times
untimes0 <- result[[3]]
cat("untimes0 =", "\n")
print(untimes0)

# Unique treatment arm event times
untimes1 <- result[[4]]
cat("untimes1 =", "\n")
print(untimes1)

# Number of unique control arm event times
nuntimes0 <- result[[5]]
cat("nuntimes0 =", "\n")
print(nuntimes0)

# Number of unique treatment arm event times
nuntimes1 <- result[[6]]
cat("nuntimes1 =", "\n")
print(nuntimes1)

# Control arm max follow time
max_follow0 <- result[[7]]
cat("max_follow0 =", "\n")
print(max_follow0)

# Treatment arm max follow time
max_follow1 <- result[[8]]
cat("max_follow1 =", "\n")
print(max_follow1)

