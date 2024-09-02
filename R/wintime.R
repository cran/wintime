#' Run a win time calculation
#'
#' This function runs one of the win time methods on observed and resampled data.
#'
#' @param type A string value indicating the desired win time method. Methods include 'ewt', 'ewtr', 'rmt', 'max', 'wtr', 'rwtr', and 'pwt'.
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
#' @param rmst_restriction Required only for `type` = 'rmt'. The RMT cutoff time (days).
#' @param seed Optional. Seed used for random number generation in resampling.
#' @return A list containing: the observed treatment effect, a vector of length `resample_num` containing resampled treatment effects, a message
#' indicating the method ran and the type of resampling done, the variance, the p-value, the total wins on treatment (pairwise methods only),
#' the total losses on treatment (pairwise methods only). A warning message will be printed for combinations of `type` and `model`/`resample`
#' that are not recommended.
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
wintime <- function(type,Time,Delta,trt,cov = NULL,model = NULL,resample = NULL,resample_num = 0,rmst_restriction = NA,seed = NA) {
  # -------------------------------
  # Validate inputs
  # -------------------------------
  type <- tolower(type)
  if (!(type %in% c("ewt","ewtr","rmt","max","wtr","rwtr","pwt"))) {
    stop(paste("Invalid type:", type, "- Please specify one of 'ewt', 'ewtr', 'rmt', 'wtr', 'rwtr', 'pwt', 'max'"))
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

  # ----------------------------
  # Set defaults
  # ----------------------------
  type <- tolower(type)
  if (is.null(model)) {
    if (type == "ewt" || type == "rmt") {
      model <- "km"
    }
    else if (type == "ewtr" || type == "max") {
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
    else if (type == "wtr" || type == "rwtr" || type == "pwt" || type == "max") {
      resample <- "boot"
    }
    else {
      resample <- NULL
    }
  }

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
      z <- km(n0,n1,m,Time,Delta)
    }
  }

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
    z_ewtr <- 0
    z_comp <- 0
  }

  # --------------------------
  # Type function calls
  # --------------------------
  obs_data <- NA
  wins <- NA
  losses <- NA
  if (type == "ewt") {
    obs_data <- EWT(m,dist_state0,dist_state1,untimes0,untimes1,nuntimes0,nuntimes1)
  }
  else if (type == "ewtr") {
    obs_data <- EWTR(n,m,nuntimes0,max_follow0,untimes0,Time,Delta,dist_state0,markov_ind,cov,trt)
  }
  else if (type == "rmt") {
    if (is.na(rmst_restriction)) {
      stop(paste("Please provide the rmst_restriction argument for this method."))
    }
    obs_data <- RMT(m,rmst_restriction,dist_state0,dist_state1,untimes0,untimes1,nuntimes0,nuntimes1)
  }
  else if (type == "wtr") {
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
    obs_data <- PWT(n,n0,n1,m,Time,Delta,trt,tau)
  }
  # type == "max"
  else {
    z_ewtr <- EWTR(n,m,nuntimes0,max_follow0,untimes0,Time,Delta,dist_state0,markov_ind,cov,trt)[[3]]
    z_comp <- COMP(n,Time,Delta,cov,trt)
    obs_data <- max(z_ewtr,z_comp)
  }


  # Set data used for p-value calculation
  data <- NA
  if (type == "wtr" || type == "rwtr") {
    data <- obs_data[[1]]
  }
  else if (type == "ewtr") {
    data <- obs_data[[1]]
  }
  else {
    data <- obs_data
  }

  # ----------------------------
  # Hypothesis testing
  # ----------------------------
  message <- NULL
  p <- NA
  resample_data <- NA
  sd_type <- NA
  variance <- NA
  if (!is.null(resample)) {
    # Run type function on bootstrap data
    if (resample == "boot" || resample == "bootstrap") {
      message <- paste("Resampling of ", type, " done on ", resample_num, " bootstraps.")
      resample_data <- bootstrap(type,rmst_restriction,model,n,m,Time,Delta,trt,cov,z_ewtr,z_comp,resample_num,seed)
    }
    # Run type function on permuted data
    else if (resample == "perm" || resample == "permutation" || resample == "perms" || resample == "permutations") {
      message <- paste("Resampling of ", type, " done on ", resample_num, " permutations.")
      resample_data <- perm(type,rmst_restriction,model,n,m,Time,Delta,trt,cov,z_ewtr,z_comp,resample_num,seed)
    }
    else {
      message <- paste("No resampling done. Analysis performed only on observed data.'",resample,"' is not a valid resampling technique. For resampling, please enter 'boot' for bootstraps, 'perm' for permutations, or allow the default resampling for this method.")
      resample <- NA
    }
  }
  if (!is.null(resample)) {
    # Calculate SD and p-value
    sd_type <- sd(resample_data)
    variance <- sd_type^2
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
    if (type == "ewtr") {
      variance <- obs_data[[2]]
      p <- 2*(1-pnorm(abs(data)/sqrt(variance),mean=0,sd=1))
    }
  }

  # ---------------------------------------------------------------------------------------------------
  # Warning message for certain combinations of type and model/resampling method
  # ---------------------------------------------------------------------------------------------------
  if ((type %in% c("ewt","rmt")) && (model == "markov" || (resample %in% c("boot","bootstrap")))) {
    warning_message <- warning("For this method, it is strongly recommended to use a KM model and to resample using permutations. These are set as defaults.")
  }
  if (type == "ewtr" && model == "km") {
    warning_message <- warning("For this method, it is strongly recommended to use a Markov model. This is set as a default.")
  }
  if (type %in% c("wtr","rwtr","pwt") && resample %in% c("perm","perms","permutation","permutations")) {
    warning_message <- warning("For this method, it is strongly recommended to resample using bootstraps. This is set as a default.")
  }
  if (type == "max" && (model == "km" || resample %in% c("perm","perms","permutation","permutations"))) {
    warning_message <- warning("For this method, it is strongly recommended to use a Markov model and to resample using bootstraps. These are set as defaults.")
  }
    return(list(data = data, resample_data = resample_data, message = message, variance = variance, p = p, wins = wins, losses = losses))
}
