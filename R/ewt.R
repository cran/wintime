#' Expected win time
#'
#' Calculates the state space probabilities using a Kaplan-Meier model (recommended) or a Markov model. This function uses these probabilities
#' to compare both arms and calculate the expected win time of the treatment arm.
#'
#' @param m The number of events in the hierarchy.
#' @param Time A m x n matrix of event times (days). Rows should represent events and columns should represent participants. Rows should be
#' in increasing order of clinical severity.
#' @param Delta A m x n matrix of event indicators Rows should represent events and columns should represent participants. Rows should be
#' in increasing order of clinical severity.
#' @param trt A vector containing treatment arm indicators (1 for treatment, 0 for control).
#' @param dist_state0 A matrix of control arm state probabilities (returned from wintime::km() or wintime::markov()).
#' @param dist_state1 A matrix of treatment arm state probabilities (returned from wintime::km() or wintime::markov()).
#' @param unique_event_times0 A vector of unique control arm event times (days) (returned from wintime::km() or wintime::markov()).
#' @param unique_event_times1 A vector of unique treatment arm event times (days) (returned from wintime::km() or wintime::markov()).
#' @param nunique_event_times0 The number of unique control arm event times (returned from wintime::km() or wintime::markov()).
#' @param nunique_event_times1 The number of unique treatment arm event times (returned from wintime::km() or wintime::markov()).
#' @param maxfollow0 The max control arm time (days) for valid nonparametric state space estimation (returned from wintime::markov() or wintime::km()).
#' @param maxfollow1 The max trt arm time (days) for valid nonparametric state space estimation (returned from wintime::markov() or wintime::km()).
#' @param max_time_inc Optional. Maximal time increment for updating multi-state distribution when parametric exponential extension models are used.
#' If unspecified, updates done at each event time in combined trial.
#' @return A list of the expected win time of the treatment arm, the components of the treatment effect.

# ------------------------
# Expected win time
# ------------------------
EWT <- function(m,Time,Delta,trt,dist_state0,dist_state1,unique_event_times0,unique_event_times1,nunique_event_times0,nunique_event_times1,maxfollow0,maxfollow1,max_time_inc) {
  # cat("Start EWT", "\n")
  components <- rep(0,m)
  n0=length(trt[trt==0])
  n1=length(trt[trt==1])

  #cat('From ewt: n0=',n0,'\n')
  #cat('From ewt: n1=',n1,'\n')
  # cat('From ewt: nunique_event_times0=',nunique_event_times0,'\n')
  # cat('From ewt: unique_event_times0=',unique_event_times0,'\n')
  # cat('From ewt: control group max event time=',max_follow0,'\n')
  # cat('From ewt: nunique_event_times1=',nunique_event_times1,'\n')
  # cat('From ewt: unique_event_times1=',unique_event_times1,'\n')
  # cat('From ewt: trt group max event time=',max_follow1,'\n')
  # cat('------------------------------------------------------------','\n')

  #----------------------------------------------------------
  # Initialize untimes as combined event times
  #----------------------------------------------------------
  untimes=rep(0,times=nunique_event_times0+nunique_event_times1)
  temp=unique(c(unique_event_times0,unique_event_times1))
  temp=temp[temp != 0]
  untimes=untimes[1:length(temp)]=temp
  nunique=length(temp)
  untimes=untimes[1:nunique]
  untimes=sort(untimes)


  # cat('-------------------------------------','\n')
  # cat('Before restriction unique _event_times=','\n')
  # print(unique_event_times)
  # cat('-------------------------------------','\n')
  # cat('dist_state0[k,] for being in state k-1=','\n')
  # print(dist_state0)
   # cat('-------------------------','\n')
   # cat('Before extension of times nunique=',nunique,'\n')
   # cat('Before extension of times untimes=','\n')
   # print(untimes)
   # cat('unique_event_times1=','\n')
   # print(unique_event_times1)
   # cat('new_dist_state1[k,] for being in state k-1=','\n')
   # print(dist_state1)
   # cat('-------------------------','\n')

  #----------------------------------------------------------
  # Get extended grid of times and update dist0, dist1
  #----------------------------------------------------------
  max_number_before_grid=length(untimes[untimes < min(maxfollow0,maxfollow1)])
  if (!is.null(max_time_inc) & !is.na(max_time_inc)) {ceiling_after_grid=ceiling((untimes[nunique]-untimes[max_number_before_grid+1])/max_time_inc)}
  if (!is.null(max_time_inc) & !is.na(max_time_inc)) {max_number_times=nunique+ceiling((untimes[nunique]-untimes[max_number_before_grid+1])/max_time_inc)}
  if (is.null(max_time_inc) | is.na(max_time_inc)) {max_number_times=nunique}
  #cat('max_number_times=',max_number_times,'\n')
  #cat('max_time_inc=',max_time_inc,'\n')
  new_untimes=rep(0,times=max_number_times)
  new_untimes[1]=untimes[1]

  new_dist0=rep(0,max_number_times*(m+1))
  dim(new_dist0)=c(m+1,max_number_times)
  new_dist1=rep(0,max_number_times*(m+1))
  dim(new_dist1)=c(m+1,max_number_times)

  if (new_untimes[1] %in% unique_event_times0) {
    new_dist0[,1]=dist_state0[,1]
    con_count=1
  } else {
    new_dist0[,1]=c(1,rep(0,times=m))
    con_count=0
  }
  if (new_untimes[1] %in% unique_event_times1) {
    new_dist1[,1]=dist_state1[,1]
    trt_count=1
  } else {
    new_dist1[,1]=c(1,rep(0,times=m))
    trt_count=0
  }
  count=1

  #  cat("-----------------------------------------------", "\n")
  #  cat("trt_count=",trt_count,"\n")
  j=1
  while (j <= nunique-1) {
    # cat("-----------------------------------------------", "\n")
    # cat("j=",j,"\n")
    # if (i<10) {
    #   cat("con_count=",con_count,"\n")
    #   cat("untimes[j+1]=",untimes[j+1],"\n")
    #   cat("untimes0[con_count]=",untimes0[con_count],"\n")
    #   cat("count=",count,"\n")
    #   cat("trt_count=",trt_count,"\n")
    #   cat("untimes1[trt_count]=",untimes1[trt_count],"\n")
    #   cat("-----------------------------------------------", "\n")
    # }

    if (!is.null(max_time_inc) & !is.na(max_time_inc)) {
      if (untimes[j+1]-new_untimes[count] > max_time_inc & new_untimes[count] >= min(maxfollow0,maxfollow1)) {
        addnum=ceiling((untimes[j+1]-untimes[j])/max_time_inc)
  #      cat("addnum triggered with addnum=",1,"\n")
        for (addcount in 1:(addnum-1)) {
          count=count+1
  #        cat("count=",count,"\n")
          new_untimes[count]=new_untimes[count-1]+(untimes[j+1]-untimes[j])/addnum
          new_dist0[,count]=new_dist0[,count-1]
          new_dist1[,count]=new_dist1[,count-1]
        }
      } else {
 #       cat("addnum not triggered=","\n")
        count=count+1
 #       cat("count=",count,"\n")
        new_untimes[count]=untimes[j+1]
        if (untimes[j+1] %in% unique_event_times0) {
          con_count=con_count+1
          new_dist0[,count]=dist_state0[,con_count]
        } else {
          new_dist0[,count]=new_dist0[,count-1]
        }
        if (untimes[j+1] %in% unique_event_times1) {
          trt_count=trt_count+1
  #        cat("trt_count=",trt_count,"\n")
          new_dist1[,count]=dist_state1[,trt_count]
  #        cat("new_dist1[,count]=",new_dist1[,count],"\n")
        } else {
          new_dist1[,count]=new_dist1[,count-1]
  #        cat("new_dist1[,count]=",new_dist1[,count],"\n")
        }
        j=j+1
      }
    } else {
      count=count+1
      new_untimes[count]=untimes[j+1]
      if (untimes[j+1] %in% unique_event_times0) {
        con_count=con_count+1
        new_dist0[,count]=dist_state0[,con_count]
      } else {
        new_dist0[,count]=new_dist0[,count-1]
      }
      if (untimes[j+1] %in% unique_event_times1) {
        trt_count=trt_count+1
        new_dist1[,count]=dist_state1[,trt_count]
      } else {
        new_dist1[,count]=new_dist1[,count-1]
      }
      j=j+1
    }
  }
  nunique=count
  untimes=new_untimes[1:nunique]
  new_dist0=new_dist0[,1:nunique]
  new_dist1=new_dist1[,1:nunique]

  #-------------------------------------------------------
  # ESTIMATE TRANSITION RATES USING SIMPLE EXPONENTIAL
  #-------------------------------------------------------
  start_time0=rep(0,times=n0)
  end_time0=rep(0,times=n0)
  start_time1=rep(0,times=n1)
  end_time1=rep(0,times=n1)
  rate0=rep(0,times=m*m)
  dim(rate0)=c(m,m)
  rate1=rep(0,times=m*m)
  dim(rate1)=c(m,m)
  Time0=Time[,trt==0,drop=FALSE]
  Delta0=Delta[,trt==0,drop=FALSE]
  Time1=Time[,trt==1, drop=FALSE]
  Delta1=Delta[,trt==1,drop=FALSE]
  #cat('TRT: dim(Time0)=',dim(Time0),'\n')
  #cat('TRT: dim(Time1)=',dim(Time1),'\n')
  for (prev_state in 1:m) {
    #     cat("-----------------------------------------------", "\n")
    #     cat("prev_state=",prev_state,"\n")
    # Control Arm
    for (end_state in prev_state:m) {
            # cat("end_state=",end_state,"\n")
      if (prev_state != 1) {
        start_time0=Time0[prev_state-1,]
        if (prev_state < m) {
          end_time0=apply(Time0[prev_state:m,],2,min)
        } else {
          end_time0=Time0[m,]
        }
        end_time0[Delta0[prev_state-1,]==0]=start_time0[Delta0[prev_state-1,]==0]
      } else {
        end_time0=apply(Time0[prev_state:m,],2,min)
      }
      #       cat("end_time=","\n")
      #       print(end_time)
      number_trans=length(end_time0[end_time0==Time0[end_state,] & Delta0[end_state,]==1 & end_time0-start_time0 > 0])
      total_duration=sum(end_time0-start_time0)
      if (total_duration > 0) {rate0[prev_state,end_state]=number_trans/total_duration}
      #       cat("number_trans=",number_trans,"\n")
      #       cat("total_duration=",total_duration,"\n")
      #       cat("rate2=",rate2[prev_state,end_state],"\n")
    }
    # Trt Arm
    for (end_state in prev_state:m) {
         #    cat("end_state=",end_state,"\n")
      if (prev_state != 1) {
        start_time1=Time1[prev_state-1,]
        if (prev_state < m) {
          end_time1=apply(Time1[prev_state:m,],2,min)
        } else {
          end_time1=Time1[m,]
        }
        end_time1[Delta1[prev_state-1,]==0]=start_time1[Delta1[prev_state-1,]==0]
      } else {
        end_time1=apply(Time1[prev_state:m,],2,min)
      }
      #       cat("end_time=","\n")
      #       print(end_time)
      number_trans=length(end_time1[end_time1==Time1[end_state,] & Delta1[end_state,]==1 & end_time1-start_time1 > 0])
      total_duration=sum(end_time1-start_time1)
      if (total_duration > 0) {rate1[prev_state,end_state]=number_trans/total_duration}
      #       cat("number_trans=",number_trans,"\n")
      #       cat("total_duration=",total_duration,"\n")
      #       cat("rate2=",rate2[prev_state,end_state],"\n")
    }
  }
  # cat("-----------------------------------------------", "\n")
  # cat("rate0=", "\n")
  # print(rate0)
  # cat("rate1=", "\n")
  # print(rate1)
  # cat("-----------------------------------------------", "\n")

  #-------------------------------------------------------
  # EXTEND NEW_DIST0 USING SIMPLE EXPONENTIAL
  #-------------------------------------------------------
  j=length(untimes[untimes <= maxfollow0])+1
  while (j <= nunique) {
    for (state in 1:m) {
      new_dist0[state,j]=new_dist0[state,j-1]*exp(-1*sum(rate0[state,])*(untimes[j]-untimes[j-1]))
      if (state > 1) {
        for (prev_state in 1:(state-1)) {
          new_dist0[state,j]=new_dist0[state,j]+new_dist0[prev_state,j-1]*(1-exp(-1*sum(rate0[prev_state,])*(untimes[j]-untimes[j-1])))*
            rate0[prev_state,state-1]/sum(rate0[prev_state,])
        }
      }
    }
    new_dist0[m+1,j]=new_dist0[m+1,j-1]
    for (prev_state in 1:m) {
      if (sum(rate0[prev_state,])>0) {
        new_dist0[m+1,j]=new_dist0[m+1,j]+new_dist0[prev_state,j-1]*(1-exp(-1*sum(rate0[prev_state,])*(untimes[j]-untimes[j-1])))*
          rate0[prev_state,m]/sum(rate0[prev_state,])
      }
    }
    j=j+1
  }

  # cat("-----------------------------------------------", "\n")
  # cat("After extension new dist0 =", "\n")
  # print(new_dist0)
  # cat("-----------------------------------------------", "\n")

   # cat("-----------------------------------------------", "\n")
   # cat("Before exponential extension nunique=",nunique,"\n")
   # cat("Before exponential extension untimes=","\n")
   # print(untimes)
   # cat("Before exponential extension new dist1 =", "\n")
   # print(new_dist1)
   # cat("-----------------------------------------------", "\n")

  #-------------------------------------------------------
  # EXTEND NEW_DIST1 USING SIMPLE EXPONENTIAL
  #-------------------------------------------------------
  j=length(untimes[untimes <= maxfollow1])+1
  while (j <= nunique) {
 #   cat("j=",j,"\n")
    for (state in 1:m) {
 #     cat("state=",state,"\n")
      new_dist1[state,j]=new_dist1[state,j-1]*exp(-1*sum(rate1[state,])*(untimes[j]-untimes[j-1]))
 #     cat("Prob you start at state and stay: new_dist1[state,j]=",new_dist1[state,j],"\n")
      if (state > 1) {
        for (prev_state in 1:(state-1)) {
 #         cat("prev_state=",prev_state,"\n")
 #         cat("new_dist1[prev_state,j-1]=",new_dist1[prev_state,j-1],"\n")
 #         cat("sum(rate1[prev_state,])=",sum(rate1[prev_state,]),"\n")
          if (sum(rate1[prev_state,])>0) {
            new_dist1[state,j]=new_dist1[state,j]+new_dist1[prev_state,j-1]*(1-exp(-1*sum(rate1[prev_state,])*(untimes[j]-untimes[j-1])))*
              rate1[prev_state,state-1]/sum(rate1[prev_state,])
          }
        }
      }
    }
    new_dist1[m+1,j]=new_dist1[m+1,j-1]
    for (prev_state in 1:m) {
      if (sum(rate1[prev_state,])>0) {
        new_dist1[m+1,j]=new_dist1[m+1,j]+new_dist1[prev_state,j-1]*(1-exp(-1*sum(rate1[prev_state,])*(untimes[j]-untimes[j-1])))*
          rate1[prev_state,m]/sum(rate1[prev_state,])
      }
    }
    j=j+1
  }

   # cat("-----------------------------------------------", "\n")
   # cat("After extension new dist1 =", "\n")
   # print(new_dist1)
   # cat("-----------------------------------------------", "\n")

  # cat('-------------------------------------','\n')
  # cat('After restriction, # of times used in calculating EWT=',nunique_event_times,'\n')
  # cat('After restriction, Largest time used in calculating EWT=',unique_event_times[nunique_event_times],'\n')
  # cat('-------------------------------------','\n')


  #-------------------------------------------------------
  # CALCULATE NET WINTIME
  # COMPARES DISTRIBUTIONS ACROSS TIMES IN COMBINED LIST
  #-------------------------------------------------------
  ewt_time=0
  j=1

  # Loop 3
  while (j < nunique) {
 #   cat("j=",j,"\n")
 #   cat("untimes[j]=",untimes[j],"\n")
 #   cat("untimes[j+1]=",untimes[j+1],"\n")
    # Add wintime
    for (event_num in 1:m) {
 #     if (j==138) {cat("event_num=",event_num,"\n")}
 #     if (j==138) {cat("new_dist0[,j]=",new_dist0[,j],"\n")}
 #     if (j==138) {cat("new_dist1[event_num,j]=",new_dist1[event_num,j],"\n")}
      ewt_time <- ewt_time + new_dist1[event_num,j] * sum(new_dist0[((event_num+1):(m+1)),j]) * (untimes[j+1]-untimes[j])
      for (k in (event_num+1):(m+1)) {
        components[k-1] <- components[k-1] + new_dist1[event_num,j] * new_dist0[k,j] * (untimes[j+1]-untimes[j])
      }
    }
 #   cat('after wins: ewt_time=',ewt_time,'\n')

    # Subtract Wintime
    for (event_num in 2:(m+1)) {
 #     if (j==138) {cat("event_num=",event_num,"\n")}
 #     if (j==138) {cat("new_dist0[,j]=",new_dist0[,j],"\n")}
 #     if (j==138) {cat("new_dist1[event_num,j]=",new_dist1[event_num,j],"\n")}
      ewt_time <- ewt_time - new_dist1[event_num,j] * sum(new_dist0[(1:(event_num-1)),j]) * (untimes[j+1]-untimes[j])
      components[event_num-1] <- components[event_num-1] - new_dist1[event_num,j] * sum(new_dist0[(1:(event_num-1)),j]) * (untimes[j+1]-untimes[j])
    }
 #   cat('after losses:ewt_time=',ewt_time,'\n')
    j=j+1
  }

  # cat('-------------------------------------','\n')
  # cat('-------------------------------------','\n')
  # cat('Final ewt_time=',ewt_time,'\n')
  # cat('-------------------------------------','\n')

  return(list(ewt_time,components))
}
