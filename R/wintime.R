#' Run a win time calculation
#'
#' This function runs one of the win time methods on observed and resampled data.
#'
#' @param type A string value indicating the desired win time method. Methods include 'ewt', 'ewtr', 'rmt', 'max', 'wtr', 'rwtr', 'pwt', 'ewtp', 'rewtp', 'ewtpr','rewtpr', and 'rpwt'.
#' @param Time A `m x n` matrix of event times (days), where `m` is the number of events in the hierarchy, and `n` is the total number of trial participants.
#' Rows should represent events and columns should represent participants. Event rows should be in increasing order of clinical severity.
#' @param Delta A `m x n` matrix of event indicators, where `m` is the number of events in the hierarchy, and `n` is the total number of trial participants.
#' Rows should represent events and columns should represent participants. Event rows should be in increasing order of clinical severity.
#' @param trt A numeric vector containing treatment arm indicators (1 for treatment, 0 for control).
#' @param cov Optional. A `n x p` matrix of covariate values, where `n` is the total number of trial participants and `p` is the number of covariates.
#' Rows should represent participants and columns should represent covariate values.
#' @param model Optional. String value. The type of model used to calculate state distributions. Options include 'km' and 'markov'.  Default depends on `type`.
#' @param resample Optional. String value. The resampling method run after the observed data calculation. Options include 'boot' and 'perm'.  Default depends on `type`.
#' @param resample_num Optional. The number of desired resamples. Default is 0.
#' @param time_restriction Required only for `type` = 'rmt','rewtp','rewtpr', and 'rpwt'. The cutoff time (days).
#' @param seed Optional. Seed used for random number generation in resampling.
#' @param nimp Required only for `type` = 'ewtpr','rewtpr'. The number of random imputations for Redistribution-to-the-right.
#' @return A list containing: the observed treatment effect, a vector of length `resample_num` containing resampled treatment effects, a message
#' indicating the method ran and the type of resampling done, the variance, the p-value, the total wins on treatment (pairwise methods only),
#' the total losses on treatment (pairwise methods only), a vector of length 'm' with the components of the treatment effect,
#' a vector of length 'm' with the variance of the components. A warning message will be printed for combinations of `type` and `model`/`resample`
#' that are not recommended.
#' @details The type parameter specifies the procedure you would like to run.
#' 'ewt' is Expected Win Time.
#' 'ewtr' is Expected Win Time Against Reference (Control Arm).
#' 'rmt' is Restricted Mean Time in Favor of Treatment.
#' 'max' is the MAX procedure (max(COMP,EWTR)).
#' 'wtr' is Win Time Ratio.
#' 'rwtr' is Restricted Win Time Ratio.
#' 'pwt' is Pairwise Win Time.
#' 'ewtp' is Expected Win Time Against Trial Population.
#' 'ewtpr' is Expected Win Time Against Trial Population With Redistribution.
#' 'rewtp' is Time Restricted Expected Win Time Against Trial Population.
#' 'rewtpr'is Time Restricted Expected Win Time Against Trial Population With Redistribution.
#' 'rpwt' is Time Restricted Pairwise Win Time.
#' @examples
#' # ------------------------------
#' # Example Inputs
#' # ------------------------------
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
#' # Covariate vectors
#' cov1 <- c(54,53,55,61,73,65,63,63,82,58,66,66)
#' cov2 <- c(34.4,32.1,34.7,54.1,55.7,43.6,32.1,44.8,85.2,12.5,33.4,21.4)
#'
#' # Covariate vectors
#' cov1 <- c(66,67,54,68,77,65,55,66,77,54)
#' cov2 <- c(3,6,4,2,3,5,8,5,3,5)
#' cov3 <- c(34.6,543.6,45.8,54.7,44.3,55.6,65.9,54.7,77.9,31.2)
#'
#' # Covariate matrix
#' cov <- cbind(cov1, cov2, cov3)
#'
#' # ------------------------
#' # wintime Examples
#' # ------------------------
#'
#' # Run WTR
#' z <- wintime("wtr", Time, Delta, trt)
#' print(z)
#'
#' # Run EWT with default settings and 10 resamples
#' z <- wintime("ewt", Time, Delta, trt, resample_num = 10)
#' print(z)
#'
#' # Run EWTR with default settings
#' z <- wintime("ewtr", Time, Delta, trt, cov = cov)
#' print(z)
#'
#' # Run EWTR with KM model (This will produce a warning message)
#' z <- wintime("ewtr", Time, Delta, trt, cov = cov, model = "km")
#' print(z)
#'
#' @import survival
#' @export

# ---------------------------------------------
# Main calling function
# ---------------------------------------------
wintime <- function(type,Time,Delta,trt,cov = NULL,model = NULL,resample = NULL,resample_num = 0,time_restriction = NA,seed = NA,nimp = 0) {
#  cat('Wintime called with model=',model,'\n')
  # -------------------------------
  # Validate inputs
  # -------------------------------
  type <- tolower(type)
  if (!(type %in% c("ewt","ewtr","rmt","max","wtr","rwtr","pwt","ewtp","rewtp","ewtpr","rewtpr","rpwt"))) {
    stop(paste("Invalid type:", type, "- Please specify one of 'ewt', 'ewtr', 'rmt', 'max', 'wtr', 'rwtr', 'pwt', 'ewtp', 'rewtp', 'ewtpr', 'rewtpr', 'rpwt'"))
  }
  if (!is.numeric(Time) || !is.numeric(Delta) || !is.numeric(trt)) {
    stop("Input vectors and matrices must be numeric.")
  }
  if (!is.matrix(Time) || !is.matrix(Delta) || !is.vector(trt)) {
    stop("Event times and event indicators must be organized into matrices. Treatment arm indicators must be contained in a single vector.")
  }
  if (nrow(Time) != nrow(Delta) || ncol(Time) != ncol(Delta) || ncol(Delta) != length(trt)) {
    stop("Input matrices and vectors must have compatible dimensions and lengths.")
  }
  if (length(trt) == 0 || length(Time) == 0 || length(Delta) == 0) {
    stop("Input matrices and vectors cannot be empty.")
  }


  # cat('set parameters','\n')
  # -----------------------------
  # Set Parameters
  # -----------------------------
  # Number of control patients
  n0 <- sum(trt == 0)

  # Number of treatment patients
  n1 <- sum(trt == 1)

  # Total patients
  n <- n0 + n1

  # Number of hierarchy levels
  m <- nrow(Time)

  # Tau
  tau <- max(Time[m, ])


  # -------------------------------
  # Further validate inputs
  # -------------------------------
  for (k in 1:(m-1)) {
    if (any(Time[k,] < Time[m,] & Delta[k,]==0)) {
      stop("Input Time of nonfatal event with Delta=0 can't preceed cens/death Time.")
    }
  }
  for (k in 1:(m-1)) {
    if (any(Time[k,]*Delta[k,] > Time[m,])) {
      stop("Input Time of nonfatal event with Delta=1 can't exceed cens/death Time.")
    }
  }

  # Set seed for RNG
  if (!is.na(seed)) {
    set.seed(seed)
  }

  # ----------------------------
  # Sort Input Data
  # ----------------------------
  sorted_indices <- order(trt)
  trt <- trt[sorted_indices]

  Time <- Time[ ,sorted_indices]
  Delta <- Delta[ ,sorted_indices]
  cov <- cov[sorted_indices, ,drop = FALSE]

  # cat('set defaults','\n')
  # ----------------------------
  # Set defaults
  # ----------------------------
  type <- tolower(type)
  if (is.null(model)) {
    if (type %in% c("ewt","rmt","ewtp","rewtp")) {
      model <- "km"
    }
    else if (type %in% c("ewtr","max","ewtpr","rewtpr")) {
      model <- "markov"
    }
    else {
      model <- NULL
    }
  }

  if (is.null(resample)) {
    if (type == "ewt" || type == "rmt") {
      resample <- "perm"
    }
    else if (type == "wtr" || type == "rwtr" || type == "pwt" || type == "max" || type == "rpwt") {
      resample <- "boot"
    }
    else {
      resample <- NULL
    }
  }

#  cat('fit model','\n')
  # ---------------------------
  # Fit model
  # ---------------------------
  z <- NULL
  markov_ind <- FALSE

  if (!is.null(model)) {
    # If a Markov model is specified, call the function to fit a Markov model
    if (model == "markov") {
      z <- markov(n0,n1,m,Time,Delta)
      markov_ind <- TRUE
    }
    # If no model is specified, use the default (KM) model
    if (model == "km") {
#      cat('call km','\n')
      z <- km(n0,n1,m,Time,Delta)
    }
  }
#  cat('is.null(z)=',is.null(z),'\n')
#  cat('model=',model,'\n')

  # Parameters for type function calls
  if (!is.null(z)) {
    dist_state0 <- z[[1]]
    dist_state1 <- z[[2]]
    untimes0 <- z[[3]]
    untimes1 <- z[[4]]
    nuntimes0 <- z[[5]]
    nuntimes1 <- z[[6]]
    max_follow0 <- z[[7]]
    max_follow1 <- z[[8]]
    dist_state2 <- z[[9]]
    untimes2 <- z[[10]]
    nuntimes2 <- z[[11]]
    max_follow2 <- z[[12]]
    if (model == "km") {
      comkm <- z[[13]]
      trtkm <- z[[14]]
      conkm <- z[[15]]
#      cat('km was fit','\n')
#      cat('conkm=',conkm,'\n')
      trans_prob2 <- array(data=0,dim=c(m,m,nuntimes2))
      trans_prob1 <- array(data=0,dim=c(m,m,nuntimes1))
      trans_prob0 <- array(data=0,dim=c(m,m,nuntimes0))
    }
    if (model == "markov") {
      trans_prob2 <- z[[13]]
      trans_prob1 <- z[[14]]
      trans_prob0 <- z[[15]]
      comkm <- array(data=0,dim=c(m,nuntimes2))
      trtkm <- array(data=0,dim=c(m,nuntimes1))
      conkm <- array(data=0,dim=c(m,nuntimes0))
    }
    z_ewtr <- 0
    z_comp <- 0
    z_ewtp <- 0
    z_rewtp <- 0
    z_ewtpr <- 0
    z_rewtpr <- 0
  }
#   cat('----------------------------------------------','\n')
#   cat('----------------------------------------------','\n')
#   cat('nuntimes2=',nuntimes2,'\n')
#   cat('untimes2=',untimes2,'\n')
#   cat('----------------------------------------------','\n')

  #cat('type function calls','\n')
  # --------------------------
  # Type function calls
  # --------------------------
  obs_data <- NA
  wins <- NA
  losses <- NA
  components <- rep(NA,m)
  components_var <- rep(NA,m)

  if (type == "ewt") {
    obs_data <- EWT(m,dist_state0,dist_state1,untimes0,untimes1,nuntimes0,nuntimes1)
  }
  else if (type == "ewtr") {
    obs_data <- EWTR(n,m,nuntimes0,max_follow0,untimes0,Time,Delta,dist_state0,markov_ind,cov,trt)
  }
  else if (type == "ewtp") {
    obs_data <- EWTP(n,m,nuntimes2,max_follow2,untimes2,Time,Delta,dist_state2,markov_ind,cov,trt)
  }
  else if (type == "rewtp") {
    if (is.na(time_restriction)) {
      stop(paste("Please provide the time_restriction argument for this method."))
    }
    obs_data <- REWTP(n,m,nuntimes2,max_follow2,untimes2,Time,Delta,dist_state2,markov_ind,cov,trt,time_restriction)
  }
  else if (type == "ewtpr") {
    obs_data <- EWTPR(n,m,nuntimes2,max_follow2,untimes2,Time,Delta,dist_state2,markov_ind,cov,trt,comkm,trans_prob2,nuntimes1,max_follow1,untimes1,dist_state1,trtkm,trans_prob1,nuntimes0,max_follow0,untimes0,dist_state0,conkm,trans_prob0,nimp)
  }
  else if (type == "rewtpr") {
    if (is.na(time_restriction)) {
      stop(paste("Please provide the time_restriction argument for this method."))
    }
    obs_data <- REWTPR(n,m,nuntimes2,max_follow2,untimes2,Time,Delta,dist_state2,markov_ind,cov,trt,comkm,trans_prob2,time_restriction,nuntimes1,max_follow1,untimes1,dist_state1,trtkm,trans_prob1,nuntimes0,max_follow0,untimes0,dist_state0,conkm,trans_prob0,nimp)
  }
  else if (type == "rmt") {
    if (is.na(time_restriction)) {
      stop(paste("Please provide the time_restriction argument for this method."))
    }
    obs_data <- RMT(m,time_restriction,dist_state0,dist_state1,untimes0,untimes1,nuntimes0,nuntimes1)
  }
  else if (type == "wtr") {
#    cat('call WTR for obs data','\n')
#    temp <- WTR(n,m,tau,trt,Time,Delta)
#    cat('return WTR for obs data','\n')
    obs_data <- WTR(n,m,tau,trt,Time,Delta)
    wins <- obs_data[[2]]
    losses <- obs_data[[3]]
  }
  else if (type == "rwtr") {
    obs_data <- RWTR(n,m,tau,trt,Time,Delta)
    wins <- obs_data[[2]]
    losses <- obs_data[[3]]
  }
  else if (type == "pwt") {
#    cat('call PWT for obs data','\n')
    obs_data <- PWT(n,n0,n1,m,Time,Delta,trt,tau)
#    cat('return PWT for obs data','\n')
  }
  else if (type == "rpwt") {
    obs_data <- RPWT(n,n0,n1,m,Time,Delta,trt,tau,time_restriction)
  }
  # type == "max"
  else {
    z_ewtr <- EWTR(n,m,nuntimes0,max_follow0,untimes0,Time,Delta,dist_state0,markov_ind,cov,trt)[[3]]
#    z_ewtp <- EWTP(n,m,nuntimes2,max_follow2,untimes2,Time,Delta,dist_state2,markov_ind,cov,trt)[[3]]
#    z_rewtp <- REWTP(n,m,nuntimes2,max_follow2,untimes2,Time,Delta,dist_state2,markov_ind,cov,trt,time_restriction)[[3]]
#    z_ewtpr <- EWTPR(n,m,nuntimes2,max_follow2,untimes2,Time,Delta,dist_state2,markov_ind,cov,trt,comkm,trans_prob)[[3]]
#    z_rewtpr <- REWTPR(n,m,nuntimes2,max_follow2,untimes2,Time,Delta,dist_state2,markov_ind,cov,trt,comkm,trans_prob,time_restriction)[[3]]
    z_comp <- COMP(n,Time,Delta,cov,trt)[[1]]
    obs_data <- max(z_ewtr,z_comp)
  }

  #cat('set data used for p-value','\n')
  # Set data used for p-value calculation
  data <- NA
  if (type == "wtr" || type == "rwtr") {
    data <- obs_data[[1]]
  }
  else if (type %in% c("ewtr","ewtp","rewtp","ewtpr","rewtpr")) {
    data <- obs_data[[1]]
    components <- obs_data[[4]]
    components_var <- obs_data[[5]]
  }
  else if (type %in% c("pwt","rpwt","rmt","ewt")) {
    data <- obs_data[[1]]
    components <- obs_data[[2]]
  }
  else {
    data <- obs_data
  }

  #cat('hypothesis tests','\n')
  #cat('data=','\n')
  #print(data)
  #cat('components=','\n')
  #print(components)
  # ----------------------------
  # Hypothesis testing
  # ----------------------------
  message <- NULL
  p <- NA
  resample_data <- NA
  resample_components <- matrix(NA,nrow=m,ncol=resample_num)
  sd_type <- NA
  sd_components <- rep(NA,m)
  variance <- NA
  #components_var <- rep(NA,m)
  if (!is.null(resample)) {
    components_var <- rep(NA,m)
    # Run type function on bootstrap data
    if (resample == "boot" || resample == "bootstrap") {
      message <- paste("Resampling of ", type, " done on ", resample_num, " bootstraps.")
      temp <- bootstrap(type,time_restriction,model,n,m,Time,Delta,trt,cov,z_ewtr,z_comp,resample_num,seed)
      resample_data <- temp[[1]]
      if (type %in% c("ewtr","ewtp","rewtp","ewtpr","rewtpr","pwt","rpwt","rmt","ewt")) {
        resample_components <- temp[[2]]
      }
      # cat('In Main resample_components=','\n')
      # print(resample_components)
    }
    # Run type function on permuted data
    else if (resample == "perm" || resample == "permutation" || resample == "perms" || resample == "permutations") {
      message <- paste("Resampling of ", type, " done on ", resample_num, " permutations.")
      temp <- perm(type,time_restriction,model,n,m,Time,Delta,trt,cov,z_ewtr,z_comp,resample_num,seed)
      resample_data <- temp[[1]]
      if (type %in% c("ewtr","ewtp","rewtp","ewtpr","rewtpr","pwt","rpwt","rmt","ewt")) {
        resample_components <- temp[[2]]
      }
    }
    else {
      message <- paste("No resampling done. Analysis performed only on observed data.'",resample,"' is not a valid resampling technique. For resampling, please enter 'boot' for bootstraps, 'perm' for permutations, or allow the default resampling for this method.")
      resample <- NA
    }
  }
  #cat('resample_data=',resample_data,'\n')
  if (!is.null(resample)) {
    # Calculate SD and p-value
    sd_type <- sd(resample_data)
    variance <- sd_type^2
    for (k in 1:m) {
      temp <- resample_components[k,]
      sd_components[k] <- sd(temp)
      # cat('sd_components[k]=',sd_components[k],'\n')
      # cat('resample_components[k,]=',resample_components[k,],'\n')
    }
    components_var <- sd_components^2
    if (type == "wtr" || type == "rwtr") {
      p <- 2*(1-pnorm(abs(data-1)/sd_type,mean=0,sd=1))
    }
    else if (type == "ewt" || type == "rmt" || type == "max") {
      p <- (sum(data < resample_data) + 1)/(resample_num + 1)
    }
    else {
      p <- 2*(1-pnorm(abs(data)/sd_type,mean=0,sd=1))
    }
  }
  else {
    # message <- paste("No resampling method specified; analysis performed on original data.")
    if (type %in% c("ewtr","ewtp","rewtp","ewtpr","rewtpr")) {
      variance <- obs_data[[2]]
      #cat('variance=',variance,'\n')
      p <- 2*(1-pnorm(abs(data)/sqrt(variance),mean=0,sd=1))
    }
  }
#  cat('warnings','\n')
  # ---------------------------------------------------------------------------------------------------
  # Warning message for certain combinations of type and model/resampling method
  # ---------------------------------------------------------------------------------------------------
  if ((type %in% c("ewt","rmt")) && (model == "markov" || (resample %in% c("boot","bootstrap")))) {
    warning_message <- warning("For this method, it is strongly recommended to use a KM model and to resample using permutations. These are set as defaults.")
  }
  if ((type %in% c("ewtr","ewtpr","rewtpr")) && model == "km") {
    warning_message <- warning("For this method, it is strongly recommended to use a Markov model. This is set as a default.")
  }
  if (type %in% c("wtr","rwtr","pwt","rpwt") && resample %in% c("perm","perms","permutation","permutations")) {
    warning_message <- warning("For this method, it is strongly recommended to resample using bootstraps. This is set as a default.")
  }
  if (type == "max" && (model == "km" || resample %in% c("perm","perms","permutation","permutations"))) {
    warning_message <- warning("For this method, it is strongly recommended to use a Markov model and to resample using bootstraps. These are set as defaults.")
  }
    return(list(data = data, resample_data = resample_data, message = message, variance = variance, p = p, wins = wins, losses = losses, components=components, components_var=components_var))
}
