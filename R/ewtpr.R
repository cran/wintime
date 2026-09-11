#' Expected win time against trial population With redistribution to the right
#'
#' Calculates the combined arm state space probabilities using a Markov model or a Kaplan-Meier model (recommended).
#' This function uses these probabilities to compare each participant's clinical state to a distribution of combined arm states.
#' Calculation is extended by redistribution-to-the-right principles
#'
#' @param n The total number of trial participants.
#' @param m The number of events in the hierarchy.
#' @param nunique2 The number of unique combined arm event times (returned from wintime::markov() or wintime::km()).
#' @param maxfollow2 The max combined arm time (days) for valid nonparametric state space estimation (returned from wintime::markov() or wintime::km()).
#' @param untimes2 A vector containing unique combined arm event times (days) (returned from wintime::markov() or wintime::km()).
#' @param Time A m x n matrix of event times (days). Rows should represent events and columns should represent participants. Rows should be
#' in increasing order of clinical severity.
#' @param Delta A m x n matrix of event indicators Rows should represent events and columns should represent participants. Rows should be
#' in increasing order of clinical severity.
#' @param dist2 A matrix of combined arm state probabilities (returned from wintime::markov() or wintime::km()).
#' @param markov_ind An indicator of the model type used (1 for Markov, 0 for Kaplan-Meier).
#' @param cov A n x p matrix of covariate values, where p is the number of covariates.
#' @param trt A vector of length n containing treatment arm indicators (1 for treatment, 0 for control).
#' @param trans_prob2 A (m x m x number of combined arm event times)
#' matrix where (i,j,k)'th value is transition probability from state i-1 to state j at k'th combined arm event time. (returned from wintime::markov() or wintime::km()).
#' @param nunique1 The number of unique trt arm event times (returned from wintime::markov() or wintime::km()).
#' @param maxfollow1 The max trt arm time (days) for valid nonparametric state space estimation (returned from wintime::markov() or wintime::km()).
#' @param untimes1 A vector containing unique trt arm event times (days) (returned from wintime::markov() or wintime::km()).
#' @param dist1 A matrix of trt arm state probabilities (returned from wintime::markov() or wintime::km()).
#' @param trans_prob1 A (m x m x number of trt arm event times)
#' matrix where (i,j,k)'th value is transition probability from state i-1 to state j at k'th trt arm event time. (returned from wintime::markov() or wintime::km()).
#' @param nunique0 The number of unique control arm event times (returned from wintime::markov() or wintime::km()).
#' @param maxfollow0 The max control arm time (days) for valid nonparametric state space estimation (returned from wintime::markov() or wintime::km()).
#' @param untimes0 A vector containing unique control arm event times (days) (returned from wintime::markov() or wintime::km()).
#' @param dist0 A matrix of control arm state probabilities (returned from wintime::markov() or wintime::km()).
#' @param trans_prob0 A (m x m x number of control arm event times)
#' matrix where (i,j,k)'th value is transition probability from state i-1 to state j at k'th control arm event time. (returned from wintime::markov() or wintime::km()).
#' @param max_time_inc Optional. Maximal time increment for updating multi-state distribution when parametric exponential extension models are used.
#' If unspecified, updates done at each event time in combined trial.
#' @return A list containing: The estimated treatment effect from the linear regression model, the variance, the Z-statistic, the components of the treatment effect, the variance of the components.

# -----------------------------------------------------------------------------
# Expected win time against trial population With Redistribution to the Right
# -----------------------------------------------------------------------------
EWTPR <- function(n,m,nunique2,maxfollow2,untimes2,Time,Delta,dist2,markov_ind,cov,trt,trans_prob2,nunique1,maxfollow1,untimes1,dist1,trans_prob1,nunique0,maxfollow0,untimes0,dist0,trans_prob0,max_time_inc) {
  time <- Time[m:1, ]
  delta <- Delta[m:1, ]
  components <- rep(NA,m)
  components_var <- rep(NA,m)
  #imp_components <- matrix(NA,nrow=m,ncol=nimp)
  #imp_components_var <- matrix(NA,nrow=m,ncol=nimp)
  #max_time <- NA
  n0=length(trt[trt==0])
  n1=length(trt[trt==1])

   # cat("-----------------------------------------------", "\n")
   # cat("Markov_ind=",markov_ind,"\n")
   # cat("-----------------------------------------------", "\n")
    # cat("nunique2 =", nunique2, "\n")
    # cat("nunique1 =", nunique1, "\n")
    # cat("nunique0 =", nunique0, "\n")
    # cat("maxfollow2 =", maxfollow2, "\n")
    # cat("maxfollow1 =", maxfollow1, "\n")
    # cat("maxfollow0 =", maxfollow0, "\n")
    #  # cat("-----------------------------------------------", "\n")
    #   cat("untimes0=", "\n")
    #   print(untimes0)
    #   cat("-----------------------------------------------", "\n")
    #   cat("untimes1=", "\n")
    #   print(untimes1)
    #   cat("-----------------------------------------------", "\n")
    #   cat("-----------------------------------------------", "\n")
    #   cat("trans_prob1[0,1,] =", "\n")
    #   print(trans_prob1[1,1,])
     # cat("untimes0=", "\n")
     # print(untimes0)
     # cat("-----------------------------------------------", "\n")
   # cat("-----------------------------------------------", "\n")
   # cat("dist2 =", "\n")
   # print(dist2)
   # cat("-----------------------------------------------", "\n")
   # cat("dist1 =", "\n")
   # print(dist1)
   # cat("-----------------------------------------------", "\n")
   # cat("dist0 =", "\n")
   # print(dist0)
   # cat("-----------------------------------------------", "\n")
#  cat("trans_prob2 =", "\n")
#  print(trans_prob2)
#  cat("-----------------------------------------------", "\n")
  # cat("trans_prob1 =", "\n")
  # print(trans_prob1)
  # cat("-----------------------------------------------", "\n")
  # cat("trans_prob0 =", "\n")
  # print(trans_prob0)
  # cat("-----------------------------------------------", "\n")

  #untimes=rep(0,nunique0+nunique1)
  #temp=unique(c(untimes0,untimes1))
  #temp=temp[temp!=0]
#  cat("temp=", "\n")
#  print(temp)
  #untimes[1:length(temp)]=temp
  #nunique=length(temp)
  #untimes=untimes[1:nunique]
  #untimes=sort(untimes)
  #----------------------------------------------------------
  # Set untimes and nunique
  #----------------------------------------------------------
  untimes=untimes2
  nunique=nunique2

   #  cat("-----------------------------------------------", "\n")
   #  cat("Odered Unique times, untimes=", "\n")
   #  print(untimes)
   #  cat("-----------------------------------------------", "\n")

  #----------------------------------------------------------
  # Get extended grid of times and update dist2
  # update trans_prob0, trans_prob1
  #----------------------------------------------------------
  max_number_before_grid=length(untimes[untimes < min(maxfollow0,maxfollow1,maxfollow2)])
  if (!is.null(max_time_inc) & !is.na(max_time_inc)) {ceiling_after_grid=ceiling((untimes[nunique]-untimes[max_number_before_grid+1])/max_time_inc)}
  if (!is.null(max_time_inc) & !is.na(max_time_inc)) {max_number_times=nunique+ceiling((untimes[nunique]-untimes[max_number_before_grid+1])/max_time_inc)}
  if (is.null(max_time_inc) | is.na(max_time_inc)) {max_number_times=nunique}
  new_untimes=rep(0,times=max_number_times)
  new_untimes[1]=untimes[1]

  new_dist2=rep(0,max_number_times*(m+1))
  dim(new_dist2)=c(m+1,max_number_times)

  new_trans_prob0=rep(0,max_number_times*m*m)
  dim(new_trans_prob0)=c(m,m,max_number_times)
  new_trans_prob1=rep(0,max_number_times*m*m)
  dim(new_trans_prob1)=c(m,m,max_number_times)

  if (new_untimes[1] %in% untimes0) {
    new_trans_prob0[,,1]=trans_prob0[,,1]
    con_count=1
  } else {
    new_trans_prob0[,,1]=0
    con_count=0
  }
  if (new_untimes[1] %in% untimes1) {
    new_trans_prob1[,,1]=trans_prob1[,,1]
    trt_count=1
  } else {
    new_trans_prob1[,,1]=0
    trt_count=0
  }
  new_dist2[,1]=dist2[,1]
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
      if (untimes[j+1]-new_untimes[count] > max_time_inc & new_untimes[count] >= min(maxfollow0,maxfollow1,maxfollow2)) {
        addnum=ceiling((untimes[j+1]-untimes[j])/max_time_inc)
     #   cat("addnum triggered with addnum=",1,"\n")
        for (addcount in 1:(addnum-1)) {
          count=count+1
     #    cat("count=",count,"\n")
          new_untimes[count]=new_untimes[count-1]+(untimes[j+1]-untimes[j])/addnum
          new_dist2[,count]=new_dist2[,count-1]
          new_trans_prob0[,,count]=0
          new_trans_prob1[,,count]=0
        }
      } else {
    #    cat("addnum not triggered=","\n")
        count=count+1
    #    cat("count=",count,"\n")
        new_untimes[count]=untimes[j+1]
        if (untimes[j+1] %in% untimes0) {
          con_count=con_count+1
          new_trans_prob0[,,count]=trans_prob0[,,con_count]
        } else {
          new_trans_prob0[,,count]=0
        }
        if (untimes[j+1] %in% untimes1) {
          trt_count=trt_count+1
          new_trans_prob1[,,count]=trans_prob1[,,trt_count]
        } else {
          new_trans_prob1[,,count]=0
        }
        if (untimes[j+1] %in% untimes0 | untimes[j+1] %in% untimes1) {
          new_dist2[,count]=dist2[,j+1]
        } else {
          new_dist2[,count]=new_dist2[,count-1]
        }
        j=j+1
      }
    } else {
      count=count+1
      new_untimes[count]=untimes[j+1]
      if (untimes[j+1] %in% untimes0) {
        con_count=con_count+1
        new_trans_prob0[,,count]=trans_prob0[,,con_count]
      } else {
        new_trans_prob0[,,count]=0
      }
      if (untimes[j+1] %in% untimes1) {
        trt_count=trt_count+1
        new_trans_prob1[,,count]=trans_prob1[,,trt_count]
      } else {
        new_trans_prob1[,,count]=0
      }
      if (untimes[j+1] %in% untimes0 | untimes[j+1] %in% untimes1) {
        new_dist2[,count]=dist2[,j+1]
      } else {
        new_dist2[,count]=new_dist2[,count-1]
      }
      j=j+1
    }
  }
  nunique=count
  untimes=new_untimes[1:nunique]
  new_dist2=new_dist2[,1:nunique]
  new_trans_prob0=new_trans_prob0[,,1:nunique]
  new_trans_prob1=new_trans_prob1[,,1:nunique]

   #
   # cat("finished putting on common set of times and adding in update times", "\n")
   # cat("-----------------------------------------------", "\n")
   # cat("nunique=",nunique, "\n")
   # cat("untimes=", "\n")
   # print(untimes)
   # cat("-----------------------------------------------", "\n")
   # cat("new_trans_prob1[0,1,] =", "\n")
   # print(trans_prob1[1,1,])
   # cat("new_trans_prob1[0,2,] =", "\n")
   # print(new_trans_prob1[1,2,])
   # cat("new_trans_prob1[0,3,] =", "\n")
   # print(new_trans_prob1[1,3,])
   # cat("new_trans_prob1[1,2,] =", "\n")
   # print(new_trans_prob1[2,2,])
   # cat("new_trans_prob1[1,3,] =", "\n")
   # print(new_trans_prob1[2,3,])
   # cat("new_trans_prob1[2,3,] =", "\n")
   # print(new_trans_prob1[3,3,])
   # cat("-----------------------------------------------", "\n")
   # cat("new_trans_prob0[0,1,] =", "\n")
   # print(new_trans_prob0[1,1,])
   # cat("new_trans_prob0[0,2,] =", "\n")
   # print(new_trans_prob0[1,2,])
   # cat("new_trans_prob0[0,3,] =", "\n")
   # print(new_trans_prob0[1,3,])
   # cat("new_trans_prob0[1,2,] =", "\n")
   # print(new_trans_prob0[2,2,])
   # cat("new_trans_prob0[1,3,] =", "\n")
   # print(new_trans_prob0[2,3,])
   # cat("new_trans_prob0[2,3,] =", "\n")
   # print(new_trans_prob0[3,3,])
   # cat("-----------------------------------------------", "\n")

  # nunique0=length(untimes[untimes <= maxfollow0])
  # nunique1=length(untimes[untimes <= maxfollow1])

   # cat("-----------------------------------------------", "\n")
   # cat("new dist2 =", "\n")
   # print(new_dist2)
   # cat("-----------------------------------------------", "\n")

  # cat("-----------------------------------------------", "\n")
  # cat("maxfollow2 =", maxfollow2, "\n")
  # cat("maxfollow0 =", maxfollow0, "\n")
  # cat("maxfollow1 =", maxfollow1, "\n")
   # cat("Combined unique times =", "\n")
   # print(untimes)
   # cat("-----------------------------------------------", "\n")
   # cat("nunique0=",nunique0,"\n")
   # cat("nunique1=",nunique1,"\n")
   # cat("-----------------------------------------------", "\n")
  # cat("unique control times =", "\n")
  # print(untimes0)
  # cat("-----------------------------------------------", "\n")
  #
  #
  # cat("-----------------------------------------------", "\n")
  # cat("Largest # of times combined arm for redist-to-the-right =", nunique2, "\n")
  # cat("Largest time combined arm for redist-to-the-right =", untimes2[nunique2], "\n")
  # cat("Largest # of times for same arm for redist-to-the-right con arm =", nunique0, "\n")
  # cat("Largest time for same arm for redist-to-the-right con arm =", untimes2[nunique0], "\n")
  # cat("Largest # of times for same arm for redist-to-the-right trt arm =", nunique1, "\n")
  # cat("Largest time for same arm for redist-to-the-right trt arm =", untimes2[nunique1], "\n")
  # cat("time_restriction =", time_restriction, "\n")
  # cat("-----------------------------------------------", "\n")

   # cat("-----------------------------------------------", "\n")
   # cat("Time[1,]=","\n")
   # print(Time[1,])
   # cat("Time[2,]=","\n")
   # print(Time[2,])
   # cat("Time[3,]=","\n")
   # print(Time[3,])
   # cat("-----------------------------------------------", "\n")
   # cat("Delta[1,]=","\n")
   # print(Delta[1,])
   # cat("-----------------------------------------------", "\n")

  #-------------------------------------------------------
  # ESTIMATE TRANSITION RATES USING SIMPLE EXPONENTIAL
  #-------------------------------------------------------
   start_time=rep(0,times=n)
   end_time=rep(0,times=n)
   start_time0=rep(0,times=n0)
   end_time0=rep(0,times=n0)
   start_time1=rep(0,times=n1)
   end_time1=rep(0,times=n1)
   rate2=rep(0,times=m*m)
   dim(rate2)=c(m,m)
   rate0=rep(0,times=m*m)
   dim(rate0)=c(m,m)
   rate1=rep(0,times=m*m)
   dim(rate1)=c(m,m)
   Time0=Time[,trt==0]
   Delta0=Delta[,trt==0]
   Time1=Time[,trt==1]
   Delta1=Delta[,trt==1]
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
# Control Arm
     for (end_state in prev_state:m) {
#       cat("end_state=",end_state,"\n")
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
       #       cat("end_state=",end_state,"\n")
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
   # cat("rate2=", "\n")
   # print(rate2)
   # cat("rate0=", "\n")
   # print(rate0)
   # cat("rate1=", "\n")
   # print(rate1)
   # cat("-----------------------------------------------", "\n")

  #-------------------------------------------------------
  # EXTEND NEW_DIST2 USING SIMPLE EXPONENTIAL
  #-------------------------------------------------------
   #old_dist2=new_dist2
   j=length(untimes[untimes <= maxfollow2])+1
   while (j <= nunique) {
     for (state in 1:m) {
       new_dist2[state,j]=new_dist2[state,j-1]*exp(-1*sum(rate2[state,])*(untimes[j]-untimes[j-1]))
       if (state > 1) {
         for (prev_state in 1:(state-1)) {
           if (sum(rate2[prev_state,])>0) {
             new_dist2[state,j]=new_dist2[state,j]+new_dist2[prev_state,j-1]*(1-exp(-1*sum(rate2[prev_state,])*(untimes[j]-untimes[j-1])))*
               rate2[prev_state,state-1]/sum(rate2[prev_state,])
           }
         }
       }
     }
     new_dist2[m+1,j]=new_dist2[m+1,j-1]
     for (prev_state in 1:m) {
       if (sum(rate2[prev_state,])>0) {
         new_dist2[m+1,j]=new_dist2[m+1,j]+new_dist2[prev_state,j-1]*(1-exp(-1*sum(rate2[prev_state,])*(untimes[j]-untimes[j-1])))*
           rate2[prev_state,m]/sum(rate2[prev_state,])
       }
     }
     j=j+1
   }

   # cat("-----------------------------------------------", "\n")
   # cat("After extension new dist2 =", "\n")
   # print(new_dist2)
   # cat("-----------------------------------------------", "\n")













  #--------------------------------------------------------
  # FOR COMPARISON WITH FORTRAN
  # set.seed(1092423368)
  # set.seed(18945611)
  # set.seed(99803332)
  # random=runif(n*nunique*nimp)
  # iran=1
  #------------------------------------------------------
  #cat("Next 10 random numbers =",random[1:10], "\n")
  #cat("-----------------------------------------------", "\n")


  #ewtpr_time=rep(0,nimp)
  #ewtpr_time_var=rep(0,nimp)

  # Set jfinalmark
  #jfinalmark=nunique-1


  # START MULTIPLE IMPUTATION LOOP
  #for (imp in 1:nimp) {

      # cat('---------------------------------','\n')
      # cat('imp=',imp,'\n')
      # cat('iran=',iran,'\n')
      # cat('next 10 random=',random[iran:(iran+9)],'\n')


    # Initialize temporary variables
  ewtpr <- rep(0,n)
  ewtpr_components <- matrix(0,nrow=m,ncol=n)

  #cat("nunique=" ,nunique, "\n")
  #cat("dim(new_dist2)=",dim(new_dist2),"\n")

  # START LOOP OVER SUBJECTS
  for (i in 1:n) {
    dist=rep(0,times=(m+1)*nunique)
    dim(dist)=c(m+1,nunique)
      # cat('---------------------------------','\n')
      # cat('---------------------------------','\n')
      # cat('i=',i,'\n')

#     if (i==100) {
       #cat("-----------------------------------------------", "\n")
       # cat("-----------------------------------------------", "\n")
       # cat("i=" ,i, "\n")
        # cat("Time[,i]" ,Time[,i], "\n")
        # cat("Delta[,i]" ,Delta[,i], "\n")
        # cat("trt[i]" ,trt[i], "\n")
        # cat("maxfollow1=" ,maxfollow1, "\n")
#       cat("nunique0=" ,nunique0, "\n")
#       cat("nunique1=" ,nunique1, "\n")
#       cat("nunique2=" ,nunique2, "\n")
# #      cat("trans_prob[1,1,]" ,trans_prob[1,1,], "\n")
# #      cat("trans_prob[1,2,]" ,trans_prob[1,2,], "\n")
# #      cat("trans_prob[1,3,]" ,trans_prob[1,3,], "\n")
# #      cat("trans_prob[2,2,]" ,trans_prob[2,2,], "\n")
# #      cat("trans_prob[2,3,]" ,trans_prob[2,3,], "\n")
# #      cat("trans_prob[3,3,]" ,trans_prob[3,3,], "\n")
#      cat("-----------------------------------------------", "\n")
#    }
    for (j in 1:(nunique-1)) {
      # if (i==100) {
      #   cat('---------------------------------','\n')
      #   cat('j=',j,'\n')
      # }
      # if (i==2 & j==12) {
      #   cat('untimes[j]=',untimes[j],'\n')
      #   cat('max(Time[,i])=',max(Time[,i]),'\n')
      # }
      if (untimes[j] >= max(Time[,i]) & Delta[m,i]==0) {
# Use RTTR: after censoring to update dist
        if (trt[i]==0) {
# Control Arm
          if (untimes[j] <= maxfollow0) {
# Use nonparametric RTTR
#            cat('Nonparametric RTTR in Control Arm','\n')
            if (j==1) {
              dist[1,j]=1-sum(new_trans_prob0[1,,j])
              for (state in 1:m) {
                dist[state+1,j]=new_trans_prob0[1,state,j]
              }
            } else {
              for (state in 1:m) {
 #               cat('state=',state,'\n')
                dist[state,j]=dist[state,j-1]*(1-sum(new_trans_prob0[state,,j]))
                if (state > 1) {
                  for (prev_state in 1:(state-1)) {
#                   cat('prev_state=',prev_state,'\n')
                    dist[state,j]=dist[state,j]+dist[prev_state,j-1]*new_trans_prob0[prev_state,state-1,j]
                  }
                }
              }
              dist[m+1,j]=dist[m+1,j-1]
              for (prev_state in 1:m) {
                dist[m+1,j]=dist[m+1,j]+dist[prev_state,j-1]*new_trans_prob0[prev_state,m,j]
              }
            }
          } else {
# Use exponential extension RTTR
            for (state in 1:m) {
              dist[state,j]=dist[state,j-1]*exp(-1*sum(rate0[state,])*(untimes[j]-untimes[j-1]))
              if (state > 1) {
                for (prev_state in 1:(state-1)) {
                  if (sum(rate0[prev_state,]) > 0) {
                    dist[state,j]=dist[state,j]+dist[prev_state,j-1]*(1-exp(-1*sum(rate0[prev_state,])*(untimes[j]-untimes[j-1])))*rate0[prev_state,state-1]/sum(rate0[prev_state,])
                  }
                }
              }
            }
            dist[m+1,j]=dist[m+1,j-1]
            for (prev_state in 1:m) {
              if (sum(rate0[prev_state,]) > 0) {
                dist[m+1,j]=dist[m+1,j]+dist[prev_state,j-1]*(1-exp(-1*sum(rate0[prev_state,])*(untimes[j]-untimes[j-1])))*rate0[prev_state,m]/sum(rate0[prev_state,])
              }
            }
          }
        } else {
# TRT Arm
          if (untimes[j] <= maxfollow1) {
# Use nonparametric RTTR
            # if (i==100) {
            #   cat('Nonparametric RTTR in Trt Arm','\n')
            # }
            if (j==1) {
              dist[1,j]=1-sum(new_trans_prob1[1,,j])
              for (state in 1:m) {
                dist[state+1,j]=new_trans_prob1[1,state,j]
              }
            } else {
              for (state in 1:m) {
                dist[state,j]=dist[state,j-1]*(1-sum(new_trans_prob1[state,,j]))
                if (state > 1) {
                  for (prev_state in 1:(state-1)) {
                    dist[state,j]=dist[state,j]+dist[prev_state,j-1]*new_trans_prob1[prev_state,state-1,j]
                  }
                }
              }
              dist[m+1,j]=dist[m+1,j-1]
              for (prev_state in 1:m) {
                dist[m+1,j]=dist[m+1,j]+dist[prev_state,j-1]*new_trans_prob1[prev_state,m,j]
              }
            }
          } else {
# Use exponential extension RTTR
            # if (i==100) {
            #   cat('Exponential RTTR in Trt Arm','\n')
            # }
            for (state in 1:m) {
              dist[state,j]=dist[state,j-1]*exp(-1*sum(rate1[state,])*(untimes[j]-untimes[j-1]))
              if (state > 1) {
                for (prev_state in 1:(state-1)) {
                  if (sum(rate1[prev_state,]) > 0) {
                    dist[state,j]=dist[state,j]+dist[prev_state,j-1]*(1-exp(-1*sum(rate1[prev_state,])*(untimes[j]-untimes[j-1])))*rate1[prev_state,state-1]/sum(rate1[prev_state,])
                  }
                }
              }
            }
            dist[m+1,j]=dist[m+1,j-1]
            for (prev_state in 1:m) {
              if (sum(rate1[prev_state,]) > 0) {
                dist[m+1,j]=dist[m+1,j]+dist[prev_state,j-1]*(1-exp(-1*sum(rate1[prev_state,])*(untimes[j]-untimes[j-1])))*rate1[prev_state,m]/sum(rate1[prev_state,])
              }
            }
          }
        }
      } else {
# Update dist based on time before censoring
        # if (i==2) {
        #   cat('untimes[j]=',untimes[j],'\n')
        #   cat('min(Time[,i])=',min(Time[,i]),'\n')
        # }
        if (untimes[j] < min(Time[,i])) {
          dist[1,j]=1
        } else {
          for (state in 1:m) {
            # if (i==2 & j==12) {
            #   cat('state=',state,'\n')
            #   cat('Time[state,i]=',Time[state,i],'\n')
            #   cat('Delta[state,i]=',Delta[state,i],'\n')
            # }
            if (untimes[j] >= Time[state,i] & Delta[state,i]==1) {
              dist[,j]=rep(0,times=m+1)
              dist[state+1,j]=1
            }
          }
          if (identical(dist[,j],rep(0,times=m+1))) {dist[1,j]=1}
        }
      }
      # if (i==100) {
      #   cat('---------------------------------','\n')
      #   cat('dist[,j]=','\n')
      #   print(dist[,j])
      #   cat('---------------------------------','\n')
      #   cat('dist2[,j]=','\n')
      #   print(dist2[,j])
      # }
# Calculate EWD for current time interval
# Wins
      for (state in 1:m) {
        for (state2 in (state+1):(m+1)) {
          ewtpr[i]=ewtpr[i]+dist[state,j]*new_dist2[state2,j]*(untimes[j+1]-untimes[j])
          ewtpr_components[state2-1,i]=ewtpr_components[state2-1,i]+dist[state,j]*new_dist2[state2,j]*(untimes[j+1]-untimes[j])
        }
      }
# Losses
      for (state in 2:(m+1)) {
        for (state2 in 1:(state-1)) {
          ewtpr[i]=ewtpr[i]-dist[state,j]*new_dist2[state2,j]*(untimes[j+1]-untimes[j])
          ewtpr_components[state-1,i]=ewtpr_components[state-1,i]-dist[state,j]*new_dist2[state2,j]*(untimes[j+1]-untimes[j])
        }
      }
      # if (i==100) {
      #   cat('---------------------------------','\n')
      #   cat('Updated ewtpr[i]=',ewtpr[i],'\n')
      # }
    }
# END Loop over event times
  }
# END Loop over Subjects

  # Get treatment estimate and variance for Z statistic
  fit_comp <- vector("list",m)
  if (!is.null(cov)) {
    fite=lm(ewtpr~trt+cov)
    for (k in 1:m) {
      outcome <- ewtpr_components[k,]
      dim(outcome) <- c(n)
      fit_comp[[k]] <- lm(outcome~trt+cov)
    }
    #    cat('----------------------------------------------------','\n')
    #    cat('fit of ewtpr with trt and baseage=','\n')
    #    print(fite)
    #    cat('----------------------------------------------------','\n')
  }
  else {
    fite <- lm(ewtpr~trt)
    for (k in 1:m) {
      outcome <- ewtpr_components[k,]
      dim(outcome) <- c(n)
      fit_comp[[k]] <- lm(outcome~trt)
    }
  }
  ewtpr_est=coef(fite)[2]
  ewtpr_var_est=vcov(fite)[2,2]

  for (k in 1:m) {
    components[k] <- coef(fit_comp[[k]])[2]
    components_var[k] <- vcov(fit_comp[[k]])[2,2]
  }

  z_ewtpr <- ewtpr_est/sqrt(ewtpr_var_est)

  # cat('----------------------------------------------------','\n')
  # cat('trt=','\n')
  # print(trt)
  # cat('ewd=','\n')
  # print(ewtpr)
  # cat('ewtpr_est=',ewtpr_est,'\n')
  # cat('ewtpr_var=',ewtpr_var_est,'\n')
  # cat('components of ewtpr_est=','\n')
  # print(components)
  # cat('----------------------------------------------------','\n')

  return(list(ewtpr_est,ewtpr_var_est,z_ewtpr,components,components_var))
}
