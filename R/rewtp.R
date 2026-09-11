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
#' @param max_time_inc Optional. Maximal time increment for updating multi-state distribution when parametric exponential extension models are used.
#' If unspecified, updates done at each event time in combined trial.
#' @return A list containing: The estimated treatment effect from the linear regression model, the variance, the Z-statistic,
#' the components of the treatment effect, the variance of the components, and the maximum time used in comparisons.

# -------------------------------------------
# Expected win time against trial population
# -------------------------------------------
REWTP <- function(n,m,nunique,maxfollow,untimes,Time,Delta,dist,markov_ind,cov,trt,time_restriction,max_time_inc) {
  time <- Time[m:1, ]
  delta <- Delta[m:1, ]
  components <- rep(NA,m)
  components_var <- rep(NA,m)
  rewtp <- rep(0,n)
  rewtp_components <- matrix(0,nrow=m,ncol=n)
  #max_time <- 0

  # cat('----------------------------------------------','\n')
  # cat('----------------------------------------------','\n')
  # cat('nunique=',nunique,'\n')
  # cat('untimes=',untimes,'\n')
  # cat('time[,3]=',time[,3],'\n')
  # cat('delta[,3]=',delta[,3],'\n')
  # cat('----------------------------------------------','\n')

  #----------------------------------------------------------
  # Get extended grid of times and update dist
  #----------------------------------------------------------
  max_number_before_grid=length(untimes[untimes < maxfollow])
  if (!is.null(max_time_inc) & !is.na(max_time_inc)) {ceiling_after_grid=ceiling((untimes[nunique]-untimes[max_number_before_grid+1])/max_time_inc)}
  if (!is.null(max_time_inc) & !is.na(max_time_inc)) {max_number_times=nunique+ceiling((untimes[nunique]-untimes[max_number_before_grid+1])/max_time_inc)}
  if (is.null(max_time_inc) | is.na(max_time_inc)) {max_number_times=nunique}
  new_untimes=rep(0,times=max_number_times)
  new_untimes[1]=untimes[1]

  new_dist=rep(0,max_number_times*(m+1))
  dim(new_dist)=c(m+1,max_number_times)

  new_dist[,1]=dist[,1]
  count=1

  #  cat("-----------------------------------------------", "\n")
  #  cat("trt_count=",trt_count,"\n")
  j=1
  while (j <= nunique-1) {
    # cat("-----------------------------------------------", "\n")
    # cat("j=",j,"\n")
    # if (i<10) {
    #   cat("untimes[j+1]=",untimes[j+1],"\n")
    #   cat("count=",count,"\n")
    #   cat("-----------------------------------------------", "\n")
    # }

    if (!is.null(max_time_inc) & !is.na(max_time_inc)) {
      if (untimes[j+1]-new_untimes[count] > max_time_inc & new_untimes[count] >= maxfollow) {
        addnum=ceiling((untimes[j+1]-untimes[j])/max_time_inc)
        #   cat("addnum triggered with addnum=",1,"\n")
        for (addcount in 1:(addnum-1)) {
          count=count+1
          #    cat("count=",count,"\n")
          new_untimes[count]=new_untimes[count-1]+(untimes[j+1]-untimes[j])/addnum
          new_dist[,count]=new_dist[,count-1]
        }
      } else {
        #    cat("addnum not triggered=","\n")
        count=count+1
        #    cat("count=",count,"\n")
        new_untimes[count]=untimes[j+1]
        new_dist[,count]=dist[,j+1]
        j=j+1
      }
    } else {
      count=count+1
      new_untimes[count]=untimes[j+1]
      new_dist[,count]=dist[,j+1]
      j=j+1
    }
  }
  nunique=count
  untimes=new_untimes[1:nunique]
  new_dist=new_dist[,1:nunique]

  #-------------------------------------------------------
  # ESTIMATE TRANSITION RATES USING SIMPLE EXPONENTIAL
  #-------------------------------------------------------
  start_time=rep(0,times=n)
  end_time=rep(0,times=n)
  rate2=rep(0,times=m*m)
  dim(rate2)=c(m,m)
  for (prev_state in 1:m) {
    #     cat("-----------------------------------------------", "\n")
    #     cat("prev_state=",prev_state,"\n")
    # Combined Arms
    for (end_state in prev_state:m) {
      #       cat("end_state=",end_state,"\n")
      if (prev_state != 1) {
        start_time=Time[prev_state-1,]
        if (prev_state < m) {
          end_time=apply(Time[prev_state:m,],2,min)
        } else {
          end_time=Time[m,]
        }
        end_time[Delta[prev_state-1,]==0]=start_time[Delta[prev_state-1,]==0]
      } else {
        end_time=apply(Time[prev_state:m,],2,min)
      }
      #       cat("end_time=","\n")
      #       print(end_time)
      number_trans=length(end_time[end_time==Time[end_state,] & Delta[end_state,]==1 & end_time-start_time > 0])
      total_duration=sum(end_time-start_time)
      if (total_duration > 0) {rate2[prev_state,end_state]=number_trans/total_duration}
      #       cat("number_trans=",number_trans,"\n")
      #       cat("total_duration=",total_duration,"\n")
      #       cat("rate2=",rate2[prev_state,end_state],"\n")
    }
  }

  #-------------------------------------------------------
  # EXTEND NEW_DIST USING SIMPLE EXPONENTIAL
  #-------------------------------------------------------
  #old_dist=new_dist
  j=length(untimes[untimes <= maxfollow])+1
  while (j <= nunique) {
    for (state in 1:m) {
      new_dist[state,j]=new_dist[state,j-1]*exp(-1*sum(rate2[state,])*(untimes[j]-untimes[j-1]))
      if (state > 1) {
        for (prev_state in 1:(state-1)) {
          if (sum(rate2[prev_state,])>0) {
            new_dist[state,j]=new_dist[state,j]+new_dist[prev_state,j-1]*(1-exp(-1*sum(rate2[prev_state,])*(untimes[j]-untimes[j-1])))*
              rate2[prev_state,state-1]/sum(rate2[prev_state,])
          }
        }
      }
    }
    new_dist[m+1,j]=new_dist[m+1,j-1]
    for (prev_state in 1:m) {
      if (sum(rate2[prev_state,])>0) {
        new_dist[m+1,j]=new_dist[m+1,j]+new_dist[prev_state,j-1]*(1-exp(-1*sum(rate2[prev_state,])*(untimes[j]-untimes[j-1])))*
          rate2[prev_state,m]/sum(rate2[prev_state,])
      }
    }
    j=j+1
  }

  #-------------------------------------------------------
  #-------------------------------------------------------
  # Start main loop over subjects
  #-------------------------------------------------------
  #-------------------------------------------------------
  for (i in 1:n) {
#    cat("i=",i,"\n")
    for (j in 1:(nunique-1)) {
#      if (i==1) {cat("j=",j,"\n")}
      if (untimes[j]>=time_restriction) {break}
      if (untimes[j+1]>time_restriction) {
        time_inc=time_restriction-untimes[j]
      } else {
        time_inc=untimes[j+1]-untimes[j]
      }
      state=0
      for (k in 0:(m-1)) {
        if (Time[m-k,i] <= untimes[j] & Delta[m-k,i]==1) {
          state=m-k
          break
        }
      }
 #     if (i==1) {cat("state=",state,"\n")}
      # Calculate EWD for current time interval
      # Wins
      if (state < m) {
        for (state2 in (state+2):(m+1)) {
  #        if (i==1) {cat("wins:state2=",state2,"\n")}
  #        if (i==1) {cat("wins:new_dist[state2,j]=",new_dist[state2,j],"\n")}
          rewtp[i]=rewtp[i]+new_dist[state2,j]*time_inc
          rewtp_components[state2-1,i]=rewtp_components[state2-1,i]+new_dist[state2,j]*time_inc
        }
      }
      # Losses
      if (state > 0) {
        for (state2 in 1:state) {
  #        if (i==1) {cat("losses:state2=",state2,"\n")}
  #        if (i==1) {cat("wins:new_dist[state2,j]=",new_dist[state2,j],"\n")}
          rewtp[i]=rewtp[i]-new_dist[state2,j]*time_inc
          rewtp_components[state,i]=rewtp_components[state,i]-new_dist[state2,j]*time_inc
        }
      }
 #     if (i==1) {cat("ewtp[i]=",ewtp[i],"\n")}
    }
 #   cat("final ewtp[i]=",ewtp[i],"\n")
  }

  # Get treatment estimate and variance for Z statistic
  fit_comp <- vector("list",m)
  if (!is.null(cov)) {
    fite=lm(rewtp~trt+cov)
    for (k in 1:m) {
      outcome <- rewtp_components[k,]
      dim(outcome) <- c(n)
      fit_comp[[k]] <- lm(outcome~trt+cov)
    }
  } else {
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
