rm(list=ls())

library(dplyr)
library(survival)
library(MASS)
library(splines)
library(dfoptim)
library(orthogonalsplinebasis) # substitute R package (splines2) and function (bSpline)

select <- dplyr::select


###### Generate the simulated data ###### 

inner_ARC_SIR_sim_data_logistic <- function(seed,
                                            N = 200,
                                            group = "control"){  
  
  ####  size of simulated data for one arm
  
  ####  period over which we are conducting inference
  spline_seq <- seq(60, 460, by=1)
  
  library(dplyr)
  library(survival)
  library(MASS)
  library(splines)
  library(dfoptim)
  library(orthogonalsplinebasis) # substitute R package (splines2) and function (bSpline)
  
  select = dplyr::select
  
  #############################################################
  #################   Functions  ##############################
  #############################################################
  {
    
    ################################################################
    #####   This performs additional formatting on the cleaned ARC data
    #####   Input:  the cleaned ARC data & the day at which to place all censored observations
    #####   Output:  a list of 3 data frames:  baseline visits;  post-baseline visits;
    #####      and a data frame for use with the `survival' package, which has rows indicating 
    #####      censoring for subjects with fewer than 4 post-baseline assessments
    formatting_fn <- function(df, last_day){
      
      K=dim(df)[1]
      
      Prev_outcome=rep(0,K)
      Prev_time=rep(0,K)
      
      for(i in 1:K){
        
        if(df$Visit_number[i]==0){Prev_outcome[i]=NA}
        if(df$Visit_number[i] !=0){Prev_outcome[i]=df$Asthma_control[i-1]}
        
      }
      
      for(i in 1:K){
        
        if(df$Visit_number[i]==0){Prev_time[i]=NA}
        if(df$Visit_number[i] !=0){Prev_time[i]=df$time[i-1]}
        
      }
      
      df=mutate(df,Prev_outcome,Prev_time)
      df=mutate(df,Lag_time=time-Prev_time,Event=1)
      
      
      ####  Each subject's baseline visit
      Baseline_df=filter(df,Visit_number==0)%>%
        select(-Event,-Prev_outcome,-Prev_time,-Lag_time)
      
      
      ####  Each subject's post-baseline visits
      Visits_df=filter(df,Visit_number !=0)%>%select(-Event)
      
      
      Survival_df=df
      
      u=sort(unique(Survival_df$elig_pid))
      H=length(u)
      
      for(h in 1:H){
        
        df_h=filter(Survival_df,elig_pid==u[h])
        v=dim(df_h)[1]
        
        if(v<5){
          
          temp_df=data.frame(elig_pid=u[h],Asthma_control=NA,
                             time=last_day,Trt=df_h$Trt[v],Visit_number=v,
                             Prev_outcome=df_h$Asthma_control[v],
                             Prev_time=df_h$time[v],Lag_time=last_day-df_h$time[v],
                             Event=0)
          
          Survival_df=rbind(Survival_df,temp_df)
          
        }
      }
      
      Survival_df=arrange(Survival_df,elig_pid)%>%filter(Visit_number != 0)
      
      return(list(Baseline_df,Visits_df,Survival_df))
      
    }
    ################################################
    
    
    ####################################################
    # Estimate the baseline intensity function in lambda(t,O(t)), the function if on page 12 
    
    #######    Function to kernel-smooth a cumulative intensity using an 
    #######        Epanechnikov kernel 
    #######    Inputs are time t, bandwidth b, and a 'survfit' object generated
    #######        by the 'survival' package
    lambda0_fn <- function(t,b,surv){
      
      times=surv$time
      cumhaz=surv$cumhaz
      
      increment=rep(0,length(times))
      increment[1]=cumhaz[1]
      for(j in 2:length(times)){ increment[j]=cumhaz[j]-cumhaz[j-1] }
      
      f_t=(1/b)*sum( 0.75*(1-( (t-times)/b )^2 )*( abs(t-times)<b )*increment)
      return(f_t)
      
    }
    ####################################################
    
    
    ####################################################
    ######  Function to generate times of a given post-baseline assessment 
    ######  Uses Ogata's Thinning Algorithm to generate times that will have
    ######    a prescribed intensity -- in this case, the estimated intensity function
    ######    from the ARC data
    
    ######  Inputs are the Visit Number, vectors of subjects' previous assessment times 
    ######    & outcomes, a value lambda_star for the thinning, and 
    ######    baseline intensity vector lambda0_t and fitted cox model for desired 
    ######    intensity function
    ######  Outputs a data frame with the assessment time, event indicator, 
    ######     previous time & outcome, lag time, and visit number for each subject
    
    ######  Lambda_star should be a value > the maximum for the prescribed intensity
    Times_gen_fn <- function(Start_vec=rep(0,N),
                              Visit=1,
                              Outcomes_vec, 
                              lambda0_t, 
                              lambda_star,
                              cox_model){ 
          
          Prev_outcome_df=data.frame(Prev_outcome=c(Outcomes_vec), 
                                     Visit_number = rep(Visit, length(Outcomes_vec)))
          
          predict_vec=predict(cox_model, newdata=Prev_outcome_df, type="lp",reference="zero")
          # "lp" means linear predictor
          # A reference of "zero" causes no centering to be done
          
          Times_vec=rep(NA,N)
          Event_vec=rep(1,N)
          
          for(i in 1:N){
            
            start=Start_vec[i]
            
            done=0
            while(done==0){
              
              ####  no kth visit if they had no (k-1)st visit
              if(Start_vec[i]==End+Visit-5){ 
                
                Times_vec[i]=End+Visit-4
                done=1
                
              }else{
                
                #####  generate a potential assessment time Tstar 
                S=rexp(1,rate=lambda_star)
                Tstar=ceiling(start+S)
                
                #####  if the candidate assessment time is too late, end with no additional
                #####     assessment time for that subject 
                if(Tstar >=  End+Visit-5){
                  
                  Times_vec[i]=End+Visit-4
                  done=1
                  
                }else{
                  
                  ######  if the candidate assessment time is not too late, accept it with 
                  ######    a probability based on the the value of the intent
                  lambda=lambda0_t[Tstar]*exp(predict_vec[i])
                  
                  U=runif(1)
                  accept=1*(U < lambda/lambda_star)
                  
                  if(accept==1){
                    
                    Times_vec[i]=Tstar
                    done=1
                    
                  }else{ start=Tstar }
                  
                }
              } 
            }
            
            #####  Mark as having no new assessment time if applicable
            if( Times_vec[i]==(End+Visit-4) ){ Event_vec[i]=0 } 
            
          }
          
          df=data.frame(id=1:N,
                        time=Times_vec,
                        event=Event_vec,
                        Prev_outcome=Outcomes_vec,
                        prev_time=Start_vec,
                        Lag_time=Times_vec-Start_vec,
                        visit=Visit)
          
          return(df)
          
        }
    ###############################
    
    Y_draw_fn_binary <- function(prob){
      
      y = c(0, 1)
      pmf = c(1 - prob, prob)
      Y = sample(y, size = 1, prob=pmf)
      
      return(Y)
      
    }
    
    # f_transfer_old <- function(data){
    #   data$record <- 1
    #   for(i in 2:nrow(data)){
    #     temp <- data[i, ]
    #     if((temp$time == temp$prev_time + 1) & is.na(temp$Prev_outcome) == T){data[i, ]$record <- 0}
    #   }
    #   data <- filter(data, record == 1)
    #   return(data)
    # } 
    
    f_transfer <- function(data){
      
      # ids <- unique(data$id)
      data2 <- data[0, ]
      
      for(i in 1:N){
        
        temp <- filter(data, id == i)
        temp1 <- filter(temp, !is.na(outcome))
        
        if(nrow(temp1) < 4 && nrow(temp1) >= 1){
          temp1 <- rbind(temp1,
                         c(i, End, 0, 
                           temp1$outcome[nrow(temp1)], # previous outcome
                           temp1$time[nrow(temp1)], # previous time
                           End - temp1$time[nrow(temp1)], # lag time
                           temp1$visit[nrow(temp1)] + 1, 
                           NA))
        }else if(nrow(temp1) == 0){
          print(i)
          temp1 <- temp[1, ]
          temp1$time <- End
          temp1$Lag_time <- End - temp$prev_time[1]
        }
        
        data2 <- rbind(data2, temp1)
        
      }
      
      return(data2)
      
    }
    
  }
  
  
  ############################################################
  ##############  Generate simulated data ##############
  ############################################################
  {
    
    ############################################################
    #######   Fit models to the ARC data
    ############################################################
    
    ARC_data <- readRDS("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_data_binary_15.rds")
    
    ## Generate binary data
    {

      # ARC_data <- readRDS("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_data.rds")
      # ARC_data_binary <- ARC_data
      # ARC_data_binary$Asthma_control <- ifelse(ARC_data_binary$Asthma_control < 1.5, 1, 0)
      # # asthma control score less than 1.5, means it can be well or partly controlled;
      # #  larger than 1.5 means it is poorly controlled.
      # saveRDS(ARC_data_binary, "/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_data_binary_15.rds")

      # ARC_data_binary$Asthma_control <- ifelse(ARC_data_binary$Asthma_control < 0.75, 1, 0)
      # # asthma control score less than 0.75, means it can be well controlled;
      # #  larger than 0.75 means it is partly or poorly controlled.
      # saveRDS(ARC_data_binary, "/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_data_binary_075.rds")

    }
    
    End <- 830
    ARC_data2 = filter(ARC_data, time <= End)
    
    ARC_formatted = formatting_fn(df=ARC_data2, last_day = End)
    
    if(group == "control"){
      
      ARC_baseline=ARC_formatted[[1]]
      ARC_baseline_hv=filter(ARC_baseline, Trt=="control") # for control group
      n_hv=dim(ARC_baseline_hv)[1]
      
      ARC_visits=ARC_formatted[[2]]
      ARC_visits_hv=filter(ARC_visits, Trt=="control")
      
      ARC_survival=ARC_formatted[[3]]
      ARC_survival_hv=filter(ARC_survival, Trt=="control")
      
    }else{
      
      ARC_baseline=ARC_formatted[[1]]
      ARC_baseline_hv=filter(ARC_baseline,Trt=="home_visits") # for treatment group
      
      n_hv=dim(ARC_baseline_hv)[1]
      ARC_visits=ARC_formatted[[2]]
      ARC_visits_hv=filter(ARC_visits,Trt=="home_visits")
      
      ARC_survival=ARC_formatted[[3]]
      ARC_survival_hv=filter(ARC_survival,Trt=="home_visits")
      
    }
    
    ###########################    Intensity model  ######################
    AG_model_hv = coxph(Surv(Prev_time,time,Event) ~ Prev_outcome + strata(Visit_number),
                      id=elig_pid,
                      data=ARC_survival_hv)
    
    ############  Estimated baseline intensities ############
    {
      AG_surv_hv = survfit(AG_model_hv,newdata=data.frame(Prev_outcome=0))
      strata = AG_surv_hv$strata
      
      v1 = strata[1]
      cumhaz_v1 = data.frame(time=AG_surv_hv$time[1:v1],cumhaz=AG_surv_hv$cumhaz[1:v1])
      base_intens_v1 = sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v1)
      
      
      v2 = strata[2]
      cumhaz_v2 = data.frame(time=AG_surv_hv$time[(v1+1):(v1+v2)],
                           cumhaz=AG_surv_hv$cumhaz[(v1+1):(v1+v2)])
      base_intens_v2 = sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v2)
      
      
      v3 = strata[3]
      cumhaz_v3 = data.frame(time=AG_surv_hv$time[(v1+v2+1):(v1+v2+v3)],
                           cumhaz=AG_surv_hv$cumhaz[(v1+v2+1):(v1+v2+v3)])
      base_intens_v3 = sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v3)
      
      
      v4 = strata[4]
      cumhaz_v4 = data.frame(time=AG_surv_hv$time[(v1+v2+v3+1):(v1+v2+v3+v4)],
                           cumhaz=AG_surv_hv$cumhaz[(v1+v2+v3+1):(v1+v2+v3+v4)])
      base_intens_v4 = sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v4)
    }
    
    
    ###### Outcome model - logistic regression model (binary)  ###### 
    {
      Outcome_model_hv <- glm(Asthma_control ~ Prev_outcome + time + Lag_time,
                              family = binomial(link = "logit"),
                              data = ARC_visits_hv)
      # check the prediction
      # pred_out <- ifelse(predict(Outcome_model_hv, type = "response") >= 0.5, 1, 0)
      # table(pred_out, ARC_visits_hv$Asthma_control)
    }
    
    
    #############################################################
    ##############     Generate simulated data    ###############
    #############################################################
    {
      
      samp = sample(n_hv, size = N, replace=TRUE)
      
      ################################################
      Y0_sim = ARC_baseline_hv[samp, ]
      Y0_sim = data.frame(Y0_sim)%>%select(Asthma_control)%>%
        rename(Prev_outcome=Asthma_control)
      
      # There is no Visit_number in the dataset
      df1_sim = Times_gen_fn(lambda0_t = base_intens_v1,
                           Outcomes_vec = Y0_sim,
                           lambda_star = 0.02,
                           cox_model = AG_model_hv)

      newdat_1 = data.frame(Prev_outcome=Y0_sim,
                            time=df1_sim$time,
                            Lag_time=df1_sim$Lag_time)
      
      ############## The logistic regression model ##############
      {
        
        prob1 <- predict(Outcome_model_hv, newdata=newdat_1, type="response")
        
        Y1_sim = rep(NA,N)
        for(n in 1:N){
          Y1_sim[n] = Y_draw_fn_binary(prob = prob1[n])
        }
        
        # Y1_sim <- ifelse(prob1 >= 0.5, 1, 0)
        
        df1_sim = mutate(df1_sim,outcome=Y1_sim)
        
        for(i in 1:N){
          if(df1_sim$event[i] == 0){ df1_sim$outcome[i] = NA }
        }
        
      }
      
      
      ################################################
      df2_sim=Times_gen_fn(Start_vec=df1_sim$time, 
                           Visit=2,
                           lambda0_t=base_intens_v2,
                           Outcomes_vec=df1_sim$outcome,
                           lambda_star=0.02,
                           cox_model=AG_model_hv)
      
      newdat_2=data.frame(Prev_outcome=Y1_sim,
                          time=df2_sim$time,
                          Lag_time=df2_sim$Lag_time)
      
      ############## The negative binomial model ##############
      {
        
        prob2 <- predict(Outcome_model_hv,newdata=newdat_2,type="response")
        
        Y2_sim = rep(NA,N)
        for(n in 1:N){
          Y2_sim[n] = Y_draw_fn_binary(prob = prob2[n])
        }
        
        # Y2_sim <- ifelse(prob2 >= 0.5, 1, 0)
        
        df2_sim=mutate(df2_sim,outcome=Y2_sim)
        
        for(i in 1:N){
          if(df2_sim$event[i]==0){ df2_sim$outcome[i]=NA }
        }
        
      }
      
      
      ################################################
      df3_sim=Times_gen_fn(Start_vec=df2_sim$time, 
                           Visit=3,
                           lambda0_t=base_intens_v3,
                           Outcomes_vec=df2_sim$outcome,
                           lambda_star=0.02,
                           cox_model=AG_model_hv)
      
      newdat_3=data.frame(Prev_outcome=Y2_sim,
                          time=df3_sim$time,
                          Lag_time=df3_sim$Lag_time)
      
      
      ############## The negative binomial model ##############
      {
        
        prob3 <- predict(Outcome_model_hv,newdata=newdat_3,type="response")
        
        Y3_sim = rep(NA,N)
        for(n in 1:N){
          Y3_sim[n] = Y_draw_fn_binary(prob = prob3[n])
        }
        
        # Y3_sim <- ifelse(prob3 >= 0.5, 1, 0)
        
        df3_sim=mutate(df3_sim,outcome=Y3_sim)
        
        for(i in 1:N){
          if(df3_sim$event[i]==0){ df3_sim$outcome[i]=NA }
        }
        
      }
      
      
      ################################################
      df4_sim=Times_gen_fn(Start_vec=df3_sim$time, 
                           Visit=4,
                           lambda0_t=base_intens_v4,
                           Outcomes_vec=df3_sim$outcome,
                           lambda_star=0.02,
                           cox_model=AG_model_hv)
      
      newdat_4=data.frame(Prev_outcome=Y3_sim,
                          time=df4_sim$time,
                          Lag_time=df4_sim$Lag_time)
      
      ############## The negative binomial model ##############
      {
        
        prob4 <- predict(Outcome_model_hv,newdata=newdat_4,type="response")
        
        Y4_sim = rep(NA,N)
        for(n in 1:N){
          Y4_sim[n] = Y_draw_fn_binary(prob = prob4[n])
        }
        
        # Y4_sim <- ifelse(prob4 >= 0.5, 1, 0)
        
        df4_sim=mutate(df4_sim,outcome=Y4_sim)
        
        for(i in 1:N){
          if(df4_sim$event[i]==0){ df4_sim$outcome[i]=NA }
        }
        
      }
      
      
      ##########################################################################
      df_sim_full = rbind(df1_sim, df2_sim, df3_sim, df4_sim) %>% arrange(id, time)
      
    }
    
  }
  
  data_trans <- f_transfer(df_sim_full)
  
  for(i in 1:N){
    if(i == 1){
      temp <- filter(data_trans, id == i)
      # print(paste0(i,",", nrow(temp)))
      df_sim_temp <- data.frame(id = i, 
                                time = temp$prev_time[1], 
                                outcome = temp$Prev_outcome[1], 
                                Visit_number = 0)
      for(j in 1:nrow(temp)){
        df_sim_temp <- rbind(df_sim_temp,
                             data.frame(id = i, 
                                        time = temp$time[j], 
                                        outcome = temp$outcome[j], 
                                        Visit_number = temp$visit[j]))
      }
    }
    else{
      temp <- filter(data_trans, id == i)
      print(paste0(i,",", nrow(temp)))
      df_sim_temp1 <- data.frame(id = i, 
                                 time = temp$prev_time[1], 
                                 outcome = temp$Prev_outcome[1], 
                                 Visit_number = 0)
      for(j in 1:nrow(temp)){
        df_sim_temp1 <- rbind(df_sim_temp1,
                              data.frame(id = i, 
                                         time = temp$time[j], 
                                         outcome = temp$outcome[j], 
                                         Visit_number = temp$visit[j]))
      }
      
      df_sim_temp <- rbind(df_sim_temp, df_sim_temp1)
      
    }
  }
  
  df_sim_temp <- filter(df_sim_temp, is.na(outcome) == F)
  
  Results <- list(data_trans, df_sim_full, df_sim_temp)
  
  return(Results) 
  
}



# sim_data <- inner_ARC_SIR_sim_data_binary(N = 30)
# saveRDS(sim_data, "/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/test_R_pack/sim_data_30_trt_binary.RData")


# Generate 500 simulated data for the treatment group and the control group separately
# Set different random seed for the treatment group and the control group

# # size 200
# set.seed(783400) # for treatment group
# set.seed(sample(1:500000000, size = 1))

set.seed(101112) # for control group
set.seed(sample(1:500000000, size = 1))

B <- 500
Index_all <- sample(1:5000000, B, replace = F)
Index <- Index_all[1:500]


library(parallel)

start <- Sys.time()
result_trt <- parallel::mclapply(Index,
                                 inner_ARC_SIR_sim_data_logistic,
                                 N = 200,
                                 group = "control", 
                                 mc.cores = 64)
end <- Sys.time()
end - start

setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/parametric_sim_data")
# saveRDS(result_trt, file = 'SIM_sim_data_treat_logit_200_logistic')
saveRDS(result_trt, file = 'SIM_sim_data_control_logit_200_logistic')

