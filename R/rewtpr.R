#' Time Restricted Expected win time against trial population With redistribution to the right
#'
#' Calculates the combined arm state space probabilities using a Markov model or a Kaplan-Meier model (recommended).
#' This function uses these probabilities to compare each participant's clinical state to a distribution of combined arm states.
#' Calculation is extended by redistribution-to-the-right principles and truncated at the user-specified time_restriction (days).
#'
#' @param n The total number of trial participants.
#' @param m The number of events in the hierarchy.
#' @param nunique2 The number of unique combined arm event times (returned from wintime::markov() or wintime::km()).
#' @param maxfollow2 The max combined arm follow up time (days) (returned from wintime::markov() or wintime::km()).
#' @param untimes2 A vector containing unique combined arm event times (days) (returned from wintime::markov() or wintime::km()).
#' @param Time A m x n matrix of event times (days). Rows should represent events and columns should represent participants. Rows should be
#' in increasing order of clinical severity.
#' @param Delta A m x n matrix of event indicators Rows should represent events and columns should represent participants. Rows should be
#' in increasing order of clinical severity.
#' @param dist2 A matrix of combined arm state probabilities (returned from wintime::markov() or wintime::km()).
#' @param markov_ind An indicator of the model type used (1 for Markov, 0 for Kaplan-Meier).
#' @param cov A n x p matrix of covariate values, where p is the number of covariates.
#' @param trt A vector of length n containing treatment arm indicators (1 for treatment, 0 for control).
#' @param comkm A m x nunique matrix of combined arm survival probabilities (returned from wintime::markov() or wintime::km()).
#' @param trans_prob2 A (m x m x number of combined arm event times)
#' matrix where (i,j,k)'th value is transition probability from state i to state j at k'th combined arm event time. (returned from wintime::markov() or wintime::km()).
#' @param time_restriction The time restriction (days) for calculation.
#' @param nunique1 The number of unique trt arm event times (returned from wintime::markov() or wintime::km()).
#' @param maxfollow1 The max trt arm follow up time (days) (returned from wintime::markov() or wintime::km()).
#' @param untimes1 A vector containing unique trt arm event times (days) (returned from wintime::markov() or wintime::km()).
#' @param dist1 A matrix of trt arm state probabilities (returned from wintime::markov() or wintime::km()).
#' @param trtkm A m x nunique matrix of trt arm survival probabilities (returned from wintime::markov() or wintime::km()).
#' @param trans_prob1 A (m x m x number of trt arm event times)
#' matrix where (i,j,k)'th value is transition probability from state i to state j at k'th trt arm event time. (returned from wintime::markov() or wintime::km()).
#' @param nunique0 The number of unique control arm event times (returned from wintime::markov() or wintime::km()).
#' @param maxfollow0 The max control arm follow up time (days) (returned from wintime::markov() or wintime::km()).
#' @param untimes0 A vector containing unique control arm event times (days) (returned from wintime::markov() or wintime::km()).
#' @param dist0 A matrix of control arm state probabilities (returned from wintime::markov() or wintime::km()).
#' @param conkm A m x nunique matrix of control arm survival probabilities (returned from wintime::markov() or wintime::km()).
#' @param trans_prob0 A (m x m x number of control arm event times)
#' matrix where (i,j,k)'th value is transition probability from state i to state j at k'th control arm event time. (returned from wintime::markov() or wintime::km()).
#'@param nimp The number of random imputations.
#' @return A list containing: The estimated treatment effect from the linear regression model, the variance, the Z-statistic, the components of the treatment effect, and the variance of the components.

# -----------------------------------------------------------------------------
# Expected win time against trial population With Redistribution to the Right
# -----------------------------------------------------------------------------
REWTPR <- function(n,m,nunique2,maxfollow2,untimes2,Time,Delta,dist2,markov_ind,cov,trt,comkm,trans_prob2,time_restriction,nunique1,maxfollow1,untimes1,dist1,trtkm,trans_prob1,nunique0,maxfollow0,untimes0,dist0,conkm,trans_prob0,nimp) {
  time <- Time[m:1, ]
  delta <- Delta[m:1, ]
  trans_prob2[trans_prob2==-1] <- 0
  trans_prob1[trans_prob1==-1] <- 0
  trans_prob0[trans_prob0==-1] <- 0
  components <- rep(NA,m)
  components_var <- rep(NA,m)
  imp_components <- matrix(NA,nrow=m,ncol=nimp)
  imp_components_var <- matrix(NA,nrow=m,ncol=nimp)

  # Get conkn,trtkm on times from untimes2
  #
  new_conkm=rep(0,nunique2*m)
  dim(new_conkm)=c(m,nunique2)
  new_trtkm=rep(0,nunique2*m)
  dim(new_trtkm)=c(m,nunique2)
  con_count=1
  trt_count=1

  for (i in 1:nunique2) {
    #    cat("-----------------------------------------------", "\n")
    #    cat("i=",i,"\n")
    #    cat("trt_count=",trt_count,"\n")
    #    cat("untimes2[i]=",untimes2[i],"\n")
    #    cat("untimes1[trt_count]=",untimes1[trt_count],"\n")
    #    cat("-----------------------------------------------", "\n")
    if (con_count < nunique0) {
      if (untimes2[i]==untimes0[con_count]) {
        new_conkm[1:m,i]=conkm[1:m,con_count]
        con_count=con_count+1
      } else {
        if (i==1) {
          new_conkm[1:m,i]=1
        } else {
          new_conkm[1:m,i]=new_conkm[1:m,i-1]
        }
      }
    }
    if (trt_count < nunique1) {
      if (untimes2[i]==untimes1[trt_count]) {
        new_trtkm[1:m,i]=trtkm[1:m,trt_count]
        trt_count=trt_count+1
      } else {
        if (i==1) {
          new_trtkm[1:m,i]=1
        } else {
          new_trtkm[1:m,i]=new_trtkm[1:m,i-1]
        }
      }
    }
  }
  # cat("-----------------------------------------------", "\n")
  # cat("new conkm and trtkm finished", "\n")
  # cat("-----------------------------------------------", "\n")
  # Get dist0,dist1 on times from untimes2
  #
  new_dist0=rep(0,nunique2*(m+1))
  dim(new_dist0)=c(m+1,nunique2)
  new_dist1=rep(0,nunique2*(m+1))
  dim(new_dist1)=c(m+1,nunique2)
  con_count=1
  trt_count=1

  for (i in 1:nunique2) {
    if (con_count < nunique0) {
      if (untimes2[i]==untimes0[con_count]) {
        new_dist0[1:(m+1),i]=dist0[1:(m+1),con_count]
        con_count=con_count+1
      } else {
        if (i==1) {
          new_dist0[1,i]=1
          new_dist0[2:(m+1),i]=0
        } else {
          new_dist0[1:(m+1),i]=new_dist0[1:(m+1),i-1]
        }
      }
    }
    if (trt_count < nunique1) {
      if (untimes2[i]==untimes1[trt_count]) {
        new_dist1[1:(m+1),i]=dist1[1:(m+1),trt_count]
        trt_count=trt_count+1
      } else {
        if (i==1) {
          new_dist1[1,i]=1
          new_dist1[2:(m+1),i]=0
        } else {
          new_dist1[1:(m+1),i]=new_dist1[1:(m+1),i-1]
        }
      }
    }
  }

  # cat("-----------------------------------------------", "\n")
  # cat("new dist0 and dist1 finished", "\n")
  # cat("-----------------------------------------------", "\n")
  # Get trans_prob0,trans_prob1 on times from untimes2
  #
  new_trans_prob0=rep(0,nunique2*m*m)
  dim(new_trans_prob0)=c(m,m,nunique2)
  new_trans_prob1=rep(0,nunique2*m*m)
  dim(new_trans_prob1)=c(m,m,nunique2)
  con_count=1
  trt_count=1

  for (i in 1:nunique2) {
    if (con_count < nunique0) {
      if (untimes2[i]==untimes0[con_count]) {
        new_trans_prob0[1:m,1:m,i]=trans_prob0[1:m,1:m,con_count]
        con_count=con_count+1
      } else {
        new_trans_prob0[1:m,1:m,i]=0
      }
    }
    if (trt_count < nunique1) {
      if (untimes2[i]==untimes1[trt_count]) {
        new_trans_prob1[1:m,1:m,i]=trans_prob1[1:m,1:m,trt_count]
        trt_count=trt_count+1
      } else {
        new_trans_prob1[1:m,1:m,i]=0
      }
    }
  }

  nunique0=length(untimes2[untimes2 <= maxfollow0])
  nunique1=length(untimes2[untimes2 <= maxfollow1])

  # cat("-----------------------------------------------", "\n")
  # cat("maxfollow2 =", maxfollow2, "\n")
  # cat("maxfollow0 =", maxfollow0, "\n")
  # cat("maxfollow1 =", maxfollow1, "\n")
  # cat("Combined unique times =", "\n")
  # print(untimes2)
  # cat("-----------------------------------------------", "\n")
  # cat("Control unique times =", "\n")
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

#--------------------------------------------------------
# FOR COMPARISON WITH FORTRAN
  random=runif(n*nunique2*nimp)
  iran=1
#--------------------------------------------------------

  rewtpr_time=rep(0,nimp)
  rewtpr_time_var=rep(0,nimp)

  # Set jfinalmark
  jfinalmark=nunique2-1


  # START MULTIPLE IMPUTATION LOOP
  for (imp in 1:nimp) {
    # cat("-----------------------------------------------", "\n")
    # cat("imp=" ,imp, "\n")

    # Initialize temporary variables
    rewtpr <- rep(0,n)
    rewtpr_components <- matrix(0,nrow=m,ncol=n)


  # START LOOP OVER SUBJECTS
  for (i in 1:n) {

#if (i<=13) {
  # cat("-----------------------------------------------", "\n")
  # cat("-----------------------------------------------", "\n")
  # cat("i=" ,i, "\n")
  # cat('iran=',iran,'\n')
  # cat("time[,i]" ,time[,i], "\n")
  # cat("delta[,i]" ,delta[,i], "\n")
  # cat("trans_prob[1,1,]" ,trans_prob[1,1,], "\n")
  # cat("trans_prob[1,2,]" ,trans_prob[1,2,], "\n")
  # cat("trans_prob[1,3,]" ,trans_prob[1,3,], "\n")
  # cat("trans_prob[2,2,]" ,trans_prob[2,2,], "\n")
  # cat("trans_prob[2,3,]" ,trans_prob[2,3,], "\n")
  # cat("trans_prob[3,3,]" ,trans_prob[3,3,], "\n")
  #cat("-----------------------------------------------", "\n")
#}

    # Set jmax
    jmax <- 0
    for (j in 1:nunique2) {
      if (untimes2[j] < time[1,i]) {
        jmax <- jmax + 1
      }
    }
    if (markov_ind == FALSE) {
      jmax <- min(maxfollow2,jmax)
    }
    if (delta[1,i] == 1 | jmax >= nunique2) {
      jmax <- nunique2 - 1
    }

    # Set jsamemark
    #if (trt[i]==1) {jsamemark=nunique1-1}
    #if (trt[i]==0) {jsamemark=nunique0-1}
    jsamemark=min(nunique0,nunique1)
    if (jsamemark == nunique2) {jsamemark=nunique2-1}
    if (jsamemark < jmax) {jsamemark=jmax}

  #  if (i<=13) {
  #    cat('i=',i,'\n')
  #    cat('jmax=',jmax,'\n')
  #     cat('jsamemark=',jsamemark,'\n')
  #    cat('jfinalmark=',jfinalmark,'\n')
  #    cat('time[,i]=',time[,i],'\n')
  #    cat('delta[,i]=',delta[,i],'\n')
  # }

    state=0
    if (jmax != 0) {
# Comparison of known state to combined arm dist
      for (j in 1:jmax) {
        if (untimes2[j]>=time_restriction) {break}
        if (untimes2[j+1]>time_restriction) {
          time_inc=time_restriction-untimes2[j]
        } else {
          time_inc=untimes2[j+1]-untimes2[j]
        }
        state <- 0
        for (current_state in m:1) {
          temp_time_index <- m - current_state + 1
          if (time[temp_time_index,i] <= untimes2[j] && delta[temp_time_index,i] == 1) {
            state <- current_state
            break
          }
        }

      # if (i==13) {
      #    cat("-----------------------------------------------", "\n")
      #    cat('j=',j,'\n')
      #    cat('state=',state,'\n')
      #    cat("untimes2[j]" ,untimes2[j], "\n")
      #    cat("untimes2[j+1]" ,untimes2[j+1], "\n")
      #    cat('dist2[,j]=','\n')
      #    print(dist2[,j])
      #    cat("-----------------------------------------------", "\n")
      #  }

        # Calculate rewtpr

        # Calculate wins
        for (state_num in 0:(m-1)) {
          if (state_num == state) {
            # Add probabilities from higher states
            for (k in (state_num+1):m) {
              rewtpr[i] <- rewtpr[i] + dist2[k+1,j] * (time_inc)
              rewtpr_components[k,i] <- rewtpr_components[k,i] + dist2[k+1,j] * (time_inc)
            }
            break
          }
        }

        # Calculate Losses
        for (current in 1:m) {
          if (current == state) {
            # Subtract probabilities from lower states
            for (k in 1:current) {
              rewtpr[i] <- rewtpr[i] - dist2[k,j] * (time_inc)
              rewtpr_components[current,i] <- rewtpr_components[current,i] - dist2[k,j] * (time_inc)
            }
            break
          }
        }
        # if (i==13) {
        #   cat('rewtpr[13]=',rewtpr[13],'\n')
        # }
      }
    }

# Initialize dist_state
    state_dist=rep(0,m+1)
    for (k in 0:m) {
      if (state==k) {state_dist[k+1]=1}
    }
    new_state_dist=state_dist

   # if (i<=13) {
   #   cat('After jmax rewtpr=',rewtpr[i],'\n')
   #   cat('state=',state,'\n')
   #   cat('state_dist=',state_dist,'\n')
   # }
    #-------------------------------------------------------------
    # Start Redistribution-to-the-right using same arm
    #-------------------------------------------------------------

    if (jmax < jfinalmark) {
    for (j in (jmax+1):jsamemark) {


      if (untimes2[j]>=time_restriction) {break}
      if (untimes2[j+1]>time_restriction) {
        time_inc=time_restriction-untimes2[j]
      } else {
        time_inc=untimes2[j+1]-untimes2[j]
      }

     # if (i==1) {
     #   cat("-----------------------------------------------", "\n")
     #   cat('j> max=',j,'\n')
     #   cat('state_dist=',state_dist,'\n')
     #   cat("untimes2[j]" ,untimes2[j], "\n")
     #   cat("untimes2[j+1]" ,untimes2[j+1], "\n")
     #   cat("trans_prob0[1,1,j]" ,new_trans_prob0[1,1,j], "\n")
     #   cat("trans_prob0[1,2,j]" ,new_trans_prob0[1,2,j], "\n")
     #   cat("-----------------------------------------------", "\n")
     # }

      # Update state_dist
      if (markov_ind == 0) {
      # KM Model
        if (trt[i]==0) {
          # RTTR by con arm

          sum=state_dist[1]
          # if (i==20 & j==37) {
          #   cat('sum=',sum,'\n')
          # }

          if (j !=1) {
            for (k in 1:m) {
              if (new_conkm[k,j-1] != 0) {
                new_state_dist[k]=sum*new_conkm[k,j]/new_conkm[k,j-1]
                sum=sum+state_dist[k]
              } else {
                new_state_dist[k]=sum*new_conkm[k,j]
                sum=sum+state_dist[k]
              }
            }
          } else {
            for (k in 1:m) {
              new_state_dist[k]=sum*new_conkm[k,j]
              sum=sum+state_dist[k]
            }
          }
          # Enforce monitonicity
          for (k in 2:m) {
            if (new_state_dist[k] < new_state_dist[k-1]) {new_state_dist[k]=new_state_dist[k-1]}
          }
          for (k in 1:(m+1)) {
            if (k==m+1) {
              state_dist[k]=1-new_state_dist[m]
            } else {
              if (k==1) {
                state_dist[k]=new_state_dist[1]
              } else {
                state_dist[k]=new_state_dist[k]-new_state_dist[k-1]
              }
            }
          }
        } else {
          # RTTR by trt arm

          sum=state_dist[1]
          # if (i==20 & j==37) {
          #   cat('sum=',sum,'\n')
          # }

          if (j !=1) {
            for (k in 1:m) {
              if (new_trtkm[k,j-1] != 0) {
                new_state_dist[k]=sum*new_trtkm[k,j]/new_trtkm[k,j-1]
                sum=sum+state_dist[k]
              } else {
                new_state_dist[k]=sum*new_trtkm[k,j]
                sum=sum+state_dist[k]
              }
            }
          } else {
            for (k in 1:m) {
              new_state_dist[k]=sum*new_trtkm[k,j]
              sum=sum+state_dist[k]
            }
          }
          # Enforce monitonicity
          for (k in 2:m) {
            if (new_state_dist[k] < new_state_dist[k-1]) {new_state_dist[k]=new_state_dist[k-1]}
          }
          for (k in 1:(m+1)) {
            if (k==m+1) {
              state_dist[k]=1-new_state_dist[m]
            } else {
              if (k==1) {
                state_dist[k]=new_state_dist[1]
              } else {
                state_dist[k]=new_state_dist[k]-new_state_dist[k-1]
              }
            }
          }
        }
      } else {
        # Markov Model

        if (trt[i]==0) {
          # RTTR by con arm

          trans_out <- array(data=0,dim=c(m))
          for (l in 1:m) {
            for (k in l:m) {
              #           if (i==1) {
              #             cat('trans_out[l] calculation for l=',l,' with k=',k,'\n')
              #           }
              trans_out[l]=trans_out[l]+new_trans_prob0[l,k,j]
            }
          }
          trans_in <- array(data=0,dim=c(m,m))
          for (l in 1:m) {
            for (k in 1:l) {
              #            if (i==1) {
              #              cat('trans_in[l] calculation for l=',l,' with k=',k,'\n')
              #            }
              trans_in[k,l]=trans_in[k,l]+new_trans_prob0[k,l,j]
            }
          }

          #        if (i==1) {
          #          cat(' trans_in=',trans_in,'\n')
          #          cat(' trans_out=',trans_out,'\n')
          #        }


          for (k in 1:(m+1)) {
            if (k <= m) {
              new_state_dist[k]=state_dist[k]*(1-trans_out[k])
            } else {
              new_state_dist[k]=state_dist[k]
            }
            if (k > 1) {
              for (l in 1:(k-1)) {
                new_state_dist[k]=new_state_dist[k]+state_dist[l]*trans_in[l,k-1]
              }
            }
          }
          state_dist=new_state_dist
        } else {
          # RTTR by trt arm

          trans_out <- array(data=0,dim=c(m))
          for (l in 1:m) {
            for (k in l:m) {
              #           if (i==1) {
              #             cat('trans_out[l] calculation for l=',l,' with k=',k,'\n')
              #           }
              trans_out[l]=trans_out[l]+new_trans_prob1[l,k,j]
            }
          }
          trans_in <- array(data=0,dim=c(m,m))
          for (l in 1:m) {
            for (k in 1:l) {
              #            if (i==1) {
              #              cat('trans_in[l] calculation for l=',l,' with k=',k,'\n')
              #            }
              trans_in[k,l]=trans_in[k,l]+new_trans_prob1[k,l,j]
            }
          }

          #        if (i==1) {
          #          cat(' trans_in=',trans_in,'\n')
          #          cat(' trans_out=',trans_out,'\n')
          #        }


          for (k in 1:(m+1)) {
            if (k <= m) {
              new_state_dist[k]=state_dist[k]*(1-trans_out[k])
            } else {
              new_state_dist[k]=state_dist[k]
            }
            if (k > 1) {
              for (l in 1:(k-1)) {
                new_state_dist[k]=new_state_dist[k]+state_dist[l]*trans_in[l,k-1]
              }
            }
          }
          state_dist=new_state_dist
        }
      }
      #END State Dist Update via RTTR using Same Arm

      # if (i==1) {
      #   cat(' Before Randomness','\n')
      #   cat(' new state_dist=',state_dist,'\n')
      #   cat('iran=',iran,'\n')
      # }

      #
      # Use Randomness to determine a state
      #
      sum=state_dist[1]
      for (k in 1:(m+1)) {
        # if (i==13 & j==31) {
        #   cat('k=',k,'\n')
        #   cat('sum=',sum,'\n')
        #   cat('iran=',iran,'\n')
        #   cat('random=',random[iran],'\n')
        # }
        #if (runif(1) < sum | k==m+1) {
        #--------------------------------------------------------
        # FOR COMPARISON WITH FORTRAN
        if (random[iran] < sum | k==m+1) {

          state_dist[1:(m+1)]=0
          state_dist[k]=1
          iran=iran+1
          break
        #} else {
        #  iran=iran+1
        }
        sum=sum+state_dist[k+1]
      }

      # if (i==1) {
      #   cat(' After Randomness','\n')
      #   cat(' new state_dist=',state_dist,'\n')
      #   cat('iran=',iran,'\n')
      # }

# Update REWTPR for Wins
      sum1=sum(dist2[2:(m+1),j])
      for (k in 1:m) {
        rewtpr[i] <- rewtpr[i] + state_dist[k] * sum1 * (time_inc)
        sum1=sum1-dist2[k+1,j]
        for (l in (k+1):(m+1)) {
          rewtpr_components[l-1,i] <- rewtpr_components[l-1,i] + state_dist[k] * dist2[l,j] * (time_inc)
        }
      }

      # Update REWTPR for Losses
      sum1=sum(dist2[1:m,j])
      for (k in (m+1):2) {
        rewtpr[i] <- rewtpr[i] - state_dist[k] * sum1 * (time_inc)
        rewtpr_components[k-1,i] <- rewtpr_components[k-1,i] - state_dist[k] * sum1 * (time_inc)
        sum1=sum1-dist2[k-1,j]
      }

      # if (i==13 & j<=40) {
      #   cat('rewtpr[13]=',rewtpr[13],'\n')
      #   cat("-----------------------------------------------", "\n")
      # }

    }
    }
    #-------------------------------------------------------------
    # End Redistribution-to-the-right using same arm
    #-------------------------------------------------------------
    # if (i==13) {
    #   cat('After jsamemax','\n')
    #   cat('rewtpr[i]=',rewtpr[i],'\n')
    #   cat('state_dist=',state_dist,'\n')
    #   cat('untimes2[j]' ,untimes2[j], "\n")
    #   cat('untimes2[j+1]' ,untimes2[j+1], "\n")
    #   cat('dist2[,j]=','\n')
    #   print(dist2[,j])
    #   cat("-----------------------------------------------", "\n")
    # }
    #------------------------------------------------------
    # Start Redistribution-to-the-right using combined arms
    #-------------------------------------------------------------

    if (jsamemark < jfinalmark) {
      for (j in (jsamemark+1):jfinalmark) {

        if (untimes2[j]>=time_restriction) {break}
        if (untimes2[j+1]>time_restriction) {
          time_inc=time_restriction-untimes2[j]
        } else {
          time_inc=untimes2[j+1]-untimes2[j]
        }

        # if (i==13) {
        #   cat('j> jsamemark=',j,'\n')
        #   cat('untimes2[j]' ,untimes2[j], "\n")
        #   cat('untimes2[j+1]' ,untimes2[j+1], "\n")
        #   cat('state_dist=',state_dist,'\n')
        #  #cat('markov_ind=',markov_ind,'\n')
        # }

        # Update state_dist
        if (markov_ind == 0) {

          # KM Model

          sum=state_dist[1]
          # if (i==20 & j==37) {
          #   cat('sum=',sum,'\n')
          # }

          if (j !=1) {
            for (k in 1:m) {
              if (comkm[k,j-1] != 0) {
                new_state_dist[k]=sum*comkm[k,j]/comkm[k,j-1]
                sum=sum+state_dist[k]
              } else {
                new_state_dist[k]=sum*comkm[k,j]
                sum=sum+state_dist[k]
              }
            }
          } else {
            for (k in 1:m) {
              new_state_dist[k]=sum*comkm[k,j]
              sum=sum+state_dist[k]
            }
          }
          # Enforce monitonicity
          for (k in 2:m) {
            if (new_state_dist[k] < new_state_dist[k-1]) {new_state_dist[k]=new_state_dist[k-1]}
          }
          for (k in 1:(m+1)) {
            if (k==m+1) {
              state_dist[k]=1-new_state_dist[m]
            } else {
              if (k==1) {
                state_dist[k]=new_state_dist[1]
              } else {
                state_dist[k]=new_state_dist[k]-new_state_dist[k-1]
              }
            }
          }
        } else {
          # Markov Model

          trans_out <- array(data=0,dim=c(m))
          for (l in 1:m) {
            for (k in l:m) {
              #           if (i==1) {
              #             cat('trans_out[l] calculation for l=',l,' with k=',k,'\n')
              #           }
              trans_out[l]=trans_out[l]+trans_prob2[l,k,j]
            }
          }
          trans_in <- array(data=0,dim=c(m,m))
          for (l in 1:m) {
            for (k in 1:l) {
              #            if (i==1) {
              #              cat('trans_in[l] calculation for l=',l,' with k=',k,'\n')
              #            }
              trans_in[k,l]=trans_in[k,l]+trans_prob2[k,l,j]
            }
          }

          #        if (i==1) {
          #          cat(' trans_in=',trans_in,'\n')
          #          cat(' trans_out=',trans_out,'\n')
          #        }


          for (k in 1:(m+1)) {
            if (k <= m) {
              new_state_dist[k]=state_dist[k]*(1-trans_out[k])
            } else {
              new_state_dist[k]=state_dist[k]
            }
            if (k > 1) {
              for (l in 1:(k-1)) {
                new_state_dist[k]=new_state_dist[k]+state_dist[l]*trans_in[l,k-1]
              }
            }
          }
          state_dist=new_state_dist

        }
        #END State Dist Update

        #cat('i=',i,'\n')
        # if (i==13) {
        #   cat('After update state_dist=',state_dist,'\n')
        #   cat('-----------------------------------------------', '\n')
        # }

        #
        # Use Randomness to determine a state
        #
        sum=state_dist[1]
        for (k in 1:(m+1)) {

          #if (runif(1) < sum | k==m) {
          # FOR COMPARISON WITH FORTRAN
          if (random[iran] < sum | k==m+1) {

            state_dist[1:(m+1)]=0
            state_dist[k]=1
            iran=iran+1
            break
          #} else {
          #  iran=iran+1
          }
          sum=sum+state_dist[k+1]
        }

        # Update EWTPR for Wins
        sum1=sum(dist2[2:(m+1),j])
        for (k in 1:m) {
          rewtpr[i] <- rewtpr[i] + state_dist[k] * sum1 * (time_inc)
          sum1=sum1-dist2[k+1,j]
          for (l in (k+1):(m+1)) {
            rewtpr_components[l-1,i] <- rewtpr_components[l-1,i] + state_dist[k] * dist2[l,j] * (time_inc)
          }
        }

        # Update EWTPR for Losses
        sum1=sum(dist2[1:m,j])
        for (k in (m+1):2) {
          rewtpr[i] <- rewtpr[i] - state_dist[k] * sum1 * (time_inc)
          rewtpr_components[k-1,i] <- rewtpr_components[k-1,i] - state_dist[k] * sum1 * (time_inc)
          sum1=sum1-dist2[k-1,j]
        }

         # if (i==13) {
         #   cat('ewtpr[i]=',ewtpr[i],'\n')
         #   cat('state_dist=',state_dist,'\n')
        #   cat('untimes[j]' ,untimes[j], "\n")
        #   cat('untimes[j+1]' ,untimes[j+1], "\n")
        #   cat('dist[,j]=','\n')
        #   print(dist[,j])
         #   cat("-----------------------------------------------", "\n")
         # }
      }
    }
    # End Redistribution-to-the-right using combined arms

  }
# END LOOP OVER SUBJECTS

  # cat('----------------------------------------------------','\n')
  # cat('rewtpr=','\n')
  # print(rewtpr[1:30])
  # cat('rewtpr_components=','\n')
  # print(rewtpr_components)
  # cat('trt=','\n')
  # print(trt)
  # cat('covariate=','\n')
  # print(cov)
  # cat('----------------------------------------------------','\n')

  # Get treatment estimate and variance for Z statistic
  fit_comp <- vector("list",m)
  if (!is.null(cov)) {
    fite=lm(rewtpr~trt+cov)
    for (k in 1:m) {
       outcome <- rewtpr_components[k,]
       dim(outcome) <- c(n)
      fit_comp[[k]] <- lm(outcome~trt+cov)
    }
#    cat('----------------------------------------------------','\n')
#    cat('fit of rewtpr with trt and baseage=','\n')
#    print(fite)
#    cat('----------------------------------------------------','\n')
  }
  else {
    fite <- lm(rewtpr~trt)
    for (k in 1:m) {
       outcome <- rewtpr_components[k,]
       dim(outcome) <- c(n)
      fit_comp[[k]] <- lm(outcome~trt)
    }
  }
  rewtpr_time[imp] <- coef(fite)[2]
  rewtpr_time_var[imp] <- vcov(fite)[2,2]
  for (k in 1:m) {
    imp_components[k,imp] <- coef(fit_comp[[k]])[2]
    imp_components_var[k,imp] <- vcov(fit_comp[[k]])[2,2]
  }

  }
  # END MULTIPLE IMPUTATION LOOP

  # cat('----------------------------------------------------','\n')
  # cat('rewtpr_time=',rewtpr_time,'\n')

  rewtpr_est=mean(rewtpr_time)
  rewtpr_var_est=mean(rewtpr_time_var)+((nimp+1)/nimp)*var(rewtpr_time)
  z_rewtpr <- rewtpr_est/sqrt(rewtpr_var_est)
  for (k in 1:m) {
    components[k] <- mean(imp_components[k,])
    components_var[k] <- mean(imp_components_var[k,])+((nimp+1)/nimp)*var(imp_components[k,])
  }

  return(list(rewtpr_est,rewtpr_var_est,z_rewtpr,components,components_var))
}
