rm(list=ls())

library(dplyr)
library(survival)
# library(MAVE)
library(MASS)
library(dfoptim)
library(orthogonalsplinebasis) 
library(BB)

select <- dplyr::select

###### 1. Logistic regression - only for binary outcome
{
  
  
  ture_value_ARC_data_binary <- function(group = "control",
                                         alpha = -0.3, 
                                         seed  = 234,
                                         N = 2*10^5){  
    
    set.seed(seed)
    
    ####  period over which we are conducting inference
    spline_seq <- seq(60, 460, by=1)
    
    
    library(dplyr)
    library(survival)
    # library(MAVE)
    library(MASS)
    library(dfoptim)
    library(orthogonalsplinebasis) 
    library(BB)
    
    select <- dplyr::select
    
    
    #############################################################
    #################   Functions  ##############################
    #############################################################
    {
      
      formatting_fn=function(df, last_day){
        
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
      
      lambda0_fn=function(t,b,surv){
        
        times=surv$time
        cumhaz=surv$cumhaz
        
        increment=rep(0,length(times))
        increment[1]=cumhaz[1]
        for(j in 2:length(times)){ increment[j]=cumhaz[j]-cumhaz[j-1] }
        
        f_t=(1/b)*sum( 0.75*(1-( (t-times)/b )^2 )*( abs(t-times)<b )*increment)
        return(f_t)
        
      }
      
      Times_gen_fn=function(Start_vec=rep(0,N),
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
      
      Y_draw_fn_binary <- function(prob){
        
        y = c(0, 1)
        pmf = c(1 - prob, prob)
        Y = sample(y, size = 1, prob=pmf)
        
        return(Y)
        
      }
      
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
      
      Cond_mean_fn_binary <- function(prob, alpha){
        
        y = c(0, 1)
        pmf = c(1 - prob, prob)
        
        E_exp_alphaY = sum( exp(alpha*y)*pmf )
        
        E_Yexp_alphaY = sum( y*exp(alpha*y)*pmf )
        
        E_Y_past = E_Yexp_alphaY / E_exp_alphaY
        
        return(list(E_Y_past, E_exp_alphaY))
        
      }
      
      s_inverse <- function(mu){
        return(log(mu/(1-mu)))
      }
      
      s_fun <- function(z){
        return(1 / (1 + exp(-z)))
      }
    }
    
    
    ############################################################
    ##############  Generate simulated data ####################
    ############################################################
    # simulate data with a sample size of n=200 for treatment group
    {
      ############################################################
      #######   Fit models to the HAP2 data
      ############################################################
      
      # ## generate count or continuous outcome
      # ARC_data=readRDS("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_data.rds")
      
      ## generate binary outcome
      ARC_data <- readRDS("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_data_binary_15.rds")
      
      End=830
      ARC_data2=filter(ARC_data,time <= End)
      
      ARC_formatted=formatting_fn(df=ARC_data2,last_day = End)
      
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
      AG_model_hv=coxph(Surv(Prev_time,time,Event) ~ Prev_outcome + strata(Visit_number),
                        id=elig_pid, 
                        data=ARC_survival_hv)
      
      ############  Estimated baseline intensities ############
      {
        AG_surv_hv=survfit(AG_model_hv,newdata=data.frame(Prev_outcome=0))
        strata=AG_surv_hv$strata
        
        v1=strata[1]
        cumhaz_v1=data.frame(time=AG_surv_hv$time[1:v1],cumhaz=AG_surv_hv$cumhaz[1:v1])
        base_intens_v1=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v1)
        
        
        v2=strata[2]
        cumhaz_v2=data.frame(time=AG_surv_hv$time[(v1+1):(v1+v2)],
                             cumhaz=AG_surv_hv$cumhaz[(v1+1):(v1+v2)])
        base_intens_v2=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v2)
        
        
        v3=strata[3]
        cumhaz_v3=data.frame(time=AG_surv_hv$time[(v1+v2+1):(v1+v2+v3)],
                             cumhaz=AG_surv_hv$cumhaz[(v1+v2+1):(v1+v2+v3)])
        base_intens_v3=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v3)
        
        
        v4=strata[4]
        cumhaz_v4=data.frame(time=AG_surv_hv$time[(v1+v2+v3+1):(v1+v2+v3+v4)],
                             cumhaz=AG_surv_hv$cumhaz[(v1+v2+v3+1):(v1+v2+v3+v4)])
        base_intens_v4=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v4)
      }
      
      ###### Outcome model - logistic regression model (binary)  ###### 
      {
        
        Outcome_model_hv <- glm(Asthma_control~Prev_outcome+time+Lag_time,
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
        
        samp=sample(n_hv, size = N, replace=TRUE)
        
        ################################################
        
        Y0_sim=ARC_baseline_hv[samp, ]
        Y0_sim=data.frame(Y0_sim)%>%select(Asthma_control)%>%
          rename(Prev_outcome=Asthma_control)
        
        # There is no Visit_number in the dataset
        df1_sim=Times_gen_fn(lambda0_t=base_intens_v1,
                             Outcomes_vec=Y0_sim,
                             lambda_star=0.02,
                             cox_model = AG_model_hv)
        
        newdat_1=data.frame(Prev_outcome=Y0_sim,
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
            if(df1_sim$event[i]==0){ df1_sim$outcome[i]=NA }
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
        
        
        ################################################
        
        df_sim_full=rbind(df1_sim,df2_sim,df3_sim,df4_sim) %>% arrange(id,time)
        
      }
      
    }
    
    
    ############################################################
    ############## Calculate the True value of beta ############
    ############################################################
    {
      
      df_sim <- f_transfer(df_sim_full) 
      
      All_visits_df_sim <- filter(df_sim, event==1) %>% rename(Asthma_control=outcome)
      Visits_df_sim <- filter(All_visits_df_sim,time >= min(spline_seq),time <= max(spline_seq))
      coefs = Outcome_model_hv$coefficients
      a = coefs[3] + coefs[4]
      
      #####  the p x 401 matrix with columns B(t), 
      basis <- SplineBasis(knots = c(60,60,60,60,260,460,460,460,460))
      B_t_matrix <- evaluate(basis, spline_seq)
      
      # the length of beta
      p <- dim(B_t_matrix)[2]
      
      #####  we approximated the integral V=\int_t B(t)B(t)'dt using sums of rectangles of width 1
      V <- GramMatrix(basis)
      V_inverse <- solve(V)
      Weights_term2 <- V_inverse%*%t(B_t_matrix)
      
      # Term2_mat <- matrix(nrow = p, ncol = N)
      EYt_true <- matrix(0, nrow = length(spline_seq), ncol = 1)
      
      start <- Sys.time()
      for(i in 1:N){
        
        print(i)
        df_i1 <- filter(df_sim, id==i)
        
        ####  predictors for time = 1
        newdat_i <- data.frame(time=1, Lag_time=1, Prev_outcome=df_i1$Prev_outcome[1])
        
        df_i2 <- filter(Visits_df_sim, id==i)
        
        times_i <- df_i2$time + 1
        
        if(length(times_i)>0){
          
          #####  predictors for times just after post-baseline visits in [a,b]
          newdat_i2 <- data.frame(time=times_i, Lag_time=1, Prev_outcome=df_i2$Asthma_control)
          
          newdat_i <- rbind(newdat_i, newdat_i2)
          
        }
        
        #####  predicted means mu(t*) just after each of participant i's assessments
        prob_i <- predict.glm(Outcome_model_hv, newdata=newdat_i, type="response")
        
        times_i = c(1, times_i)
        
        time_means_fn  <- function(t){
          
          #####  for each t, mu(t) = mu(t*)exp(a Lag), where t* is the t* is the latest
          #####    day-after-an-assessment day on or before day t
          m <- max(which(times_i<=t))
          s <- t - times_i[m]
          delta <- a * s
          
          temp = Cond_mean_fn_binary(prob = (exp(delta) * prob_i[m]) / (1 -  prob_i[m] + exp(delta) * prob_i[m]),
                                     alpha = alpha)
          
          return(temp[[1]])
          
        }
        
        Time_means <- matrix(sapply(spline_seq, time_means_fn), ncol = 1)
        
        EYt_true <- EYt_true + Time_means / N 
        
      }
      end <- Sys.time()
      end - start
      
    }
    
    Beta_true_logit <- V_inverse %*% (1/2 * (t(B_t_matrix[-length(spline_seq), ]) %*% s_inverse(EYt_true[-length(spline_seq), ]) +
                                               t(B_t_matrix[-1, ]) %*% s_inverse(EYt_true[-1, ])))
    
    
    B3 <- evaluate(basis, 90)
    B6 <- evaluate(basis, 180)
    B9 <- evaluate(basis, 270)
    B12 <- evaluate(basis, 360)
    
    month3_orig <- EYt_true[which(spline_seq == 90)]
    month6_orig <- EYt_true[which(spline_seq == 180)]
    month9_orig <- EYt_true[which(spline_seq == 270)]
    month12_orig <- EYt_true[which(spline_seq == 360)]
    
    month3_logit  <- s_fun(B3 %*% Beta_true_logit)
    month6_logit  <- s_fun(B6 %*% Beta_true_logit)
    month9_logit  <- s_fun(B9 %*% Beta_true_logit)
    month12_logit  <- s_fun(B12 %*% Beta_true_logit)
    
    
    return(list(alpha, 
                Beta_true_logit,
                c(month3_orig, month6_orig, month9_orig, month12_orig),
                c(month3_logit, month6_logit, month9_logit, month12_logit)
    ))
    
  }
  
  
  ### treatment group
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_logit_logistic")
  
  result_n08 <-  ture_value_ARC_data_binary(group = "treat", alpha = -0.8, seed  = 234, N = 2*10^5)
  saveRDS(result_n08, file = 'treat_true_neg08_200000')
  
  result_n04 <-  ture_value_ARC_data_binary(group = "treat", alpha = -0.4, seed  = 345, N = 2*10^5)
  saveRDS(result_n04, file = 'treat_true_neg04_200000')
  
  result_0 <-  ture_value_ARC_data_binary(group = "treat", alpha = 0, seed  = 456, N = 2*10^5)
  saveRDS(result_0, file = 'treat_true_0_200000')
  
  result_04 <-  ture_value_ARC_data_binary(group = "treat", alpha = 0.4, seed  = 567, N = 2*10^5)
  saveRDS(result_04, file = 'treat_true_04_200000')

  result_08 <-  ture_value_ARC_data_binary(group = "treat", alpha = 0.8, seed  = 678, N = 2*10^5)
  saveRDS(result_08, file = 'treat_true_08_200000')
  
  result_12 <-  ture_value_ARC_data_binary(group = "treat", alpha = 1.2, seed  = 123, N = 2*10^5)
  saveRDS(result_12, file = 'treat_true_12_200000')
  
  
  
  ### control group
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_logit_logistic")
  
  result_n08 <-  ture_value_ARC_data_binary(group = "control", alpha = -0.8, seed  = 234, N = 2*10^5)
  saveRDS(result_n08, file = 'control_true_neg08_200000')
  
  result_n04 <-  ture_value_ARC_data_binary(group = "control", alpha = -0.4, seed  = 345, N = 2*10^5)
  saveRDS(result_n04, file = 'control_true_neg04_200000')
  
  result_0 <-  ture_value_ARC_data_binary(group = "control", alpha = 0, seed  = 456, N = 2*10^5)
  saveRDS(result_0, file = 'control_true_0_200000')
  
  result_04 <-  ture_value_ARC_data_binary(group = "control", alpha = 0.4, seed  = 567, N = 2*10^5)
  saveRDS(result_04, file = 'control_true_04_200000')
  
  result_08 <-  ture_value_ARC_data_binary(group = "control", alpha = 0.8, seed  = 678, N = 2*10^5)
  saveRDS(result_08, file = 'control_true_08_200000')
  
  result_12 <-  ture_value_ARC_data_binary(group = "control", alpha = 1.2, seed  = 123, N = 2*10^5)
  saveRDS(result_12, file = 'control_true_12_200000')
  
  
  
  #####
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_logit_logistic")
  
  trt_neg08 <- readRDS(file = 'treat_true_neg08_200000') 
  trt_neg04 <- readRDS(file = 'treat_true_neg04_200000') 
  trt_0 <- readRDS(file = 'treat_true_0_200000') 
  trt_04 <- readRDS(file = 'treat_true_04_200000')
  trt_08 <- readRDS(file = 'treat_true_08_200000')
  trt_12 <- readRDS(file = 'treat_true_12_200000')
  
  crl_neg08 <- readRDS(file = 'control_true_neg08_200000') 
  crl_neg04 <- readRDS(file = 'control_true_neg04_200000') 
  crl_0 <- readRDS(file = 'control_true_0_200000') 
  crl_04 <- readRDS(file = 'control_true_04_200000')
  crl_08 <- readRDS(file = 'control_true_08_200000')
  crl_12 <- readRDS(file = 'control_true_12_200000')
  
  
  # link
  temp1 <- trt_neg08[[4]] |> round(5) |> t()
  temp2 <- trt_neg04[[4]] |> round(5) |> t()
  temp3 <- trt_0[[4]] |> round(5) |> t()
  temp4 <- trt_04[[4]] |> round(5) |> t()
  temp5 <- trt_08[[4]] |> round(5) |> t()
  temp6 <- trt_12[[4]] |> round(5) |> t()
  
  # month <- 3
  c(temp1[1], temp2[1], temp3[1], temp4[1], temp5[1], temp6[1])
  month3_treat <- c(0.17528, 0.23253, 0.29979, 0.37766, 0.46339, 0.55269)

  # month <- 6
  c(temp1[2], temp2[2], temp3[2], temp4[2], temp5[2], temp6[2])
  month6_treat <- c(0.19083, 0.25123, 0.31990, 0.39857, 0.48426, 0.57133)

  # month <- 9
  c(temp1[3], temp2[3], temp3[3], temp4[3], temp5[3], temp6[3])
  month9_treat <- c(0.19837, 0.26035, 0.33039, 0.40918, 0.49422, 0.58054)

  # month <- 12
  c(temp1[4], temp2[4], temp3[4], temp4[4], temp5[4], temp6[4])
  month12_treat <- c(0.20356, 0.26597, 0.33736, 0.41639, 0.50118, 0.58660)

  
  temp1 <- crl_neg08[[4]] |> round(5) |> t()
  temp2 <- crl_neg04[[4]] |> round(5) |> t()
  temp3 <- crl_0[[4]] |> round(5) |> t()
  temp4 <- crl_04[[4]] |> round(5) |> t()
  temp5 <- crl_08[[4]] |> round(5) |> t()
  temp6 <- crl_12[[4]] |> round(3) |> t()
  
  c(temp1[1], temp2[1], temp3[1], temp4[1], temp5[1], temp6[1])
  month3_control <- c(0.13761, 0.18518, 0.24377, 0.31265, 0.39025, 0.47700)
  
  c(temp1[2], temp2[2], temp3[2], temp4[2], temp5[2], temp6[2])
  month6_control <- c(0.14696, 0.19688, 0.25781, 0.32844, 0.40780, 0.49500)
  
  c(temp1[3], temp2[3], temp3[3], temp4[3], temp5[3], temp6[3])
  month9_control <- c(0.15750, 0.20958, 0.27300, 0.34542, 0.42643, 0.51300)
  
  c(temp1[4], temp2[4], temp3[4], temp4[4], temp5[4], temp6[4])
  month12_control <- c(0.17020, 0.22499, 0.29071, 0.36504, 0.44775, 0.53400)
  
  
  month3_treat <- c(0.17528, 0.23253, 0.29979, 0.37766, 0.46339, 0.55269)
  month6_treat <- c(0.19083, 0.25123, 0.31990, 0.39857, 0.48426, 0.57133)
  month9_treat <- c(0.19837, 0.26035, 0.33039, 0.40918, 0.49422, 0.58054)
  month12_treat <- c(0.20356, 0.26597, 0.33736, 0.41639, 0.50118, 0.58660)
  
  month3_control <- c(0.13761, 0.18518, 0.24377, 0.31265, 0.39025, 0.47700)
  month6_control <- c(0.14696, 0.19688, 0.25781, 0.32844, 0.40780, 0.49500)
  month9_control <- c(0.15750, 0.20958, 0.27300, 0.34542, 0.42643, 0.51300)
  month12_control <- c(0.17020, 0.22499, 0.29071, 0.36504, 0.44775, 0.53400)
  
}


###### 2. Negative binomial model  - only for log outcome
{
  
  
  ture_value_ARC_data_neg <- function(group = "control",
                                      alpha = -0.3, 
                                      seed  = 234,
                                      N = 2*10^5){  
    
    set.seed(seed)
    
    ####  period over which we are conducting inference
    spline_seq <- seq(60, 460, by=1)
    
    
    library(dplyr)
    library(survival)
    # library(MAVE)
    library(MASS)
    library(dfoptim)
    library(orthogonalsplinebasis) 
    library(BB)
    
    select <- dplyr::select
    
    
    #############################################################
    #################   Functions  ##############################
    #############################################################
    {
      
      formatting_fn=function(df, last_day){
        
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
      
      lambda0_fn=function(t,b,surv){
        
        times=surv$time
        cumhaz=surv$cumhaz
        
        increment=rep(0,length(times))
        increment[1]=cumhaz[1]
        for(j in 2:length(times)){ increment[j]=cumhaz[j]-cumhaz[j-1] }
        
        f_t=(1/b)*sum( 0.75*(1-( (t-times)/b )^2 )*( abs(t-times)<b )*increment)
        return(f_t)
        
      }
      
      Times_gen_fn=function(Start_vec=rep(0,N),
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
      
      Y_draw_fn=function(mu,theta){
        
        pmf_fn=function(y){
          
          f_y=dnbinom(y,size=theta,prob=theta/(theta+mu))
          return(f_y)
          
        }
        pmf=matrix(sapply(0:36,pmf_fn),ncol=1)
        
        y=seq(0,6,by=1/6)
        Y=sample(y, size=1, prob=pmf)
        
        return(Y)
        
      }
      
      Cond_mean_fn=function(mu,theta,alpha){
        
        pmf_fn=function(z){
          
          f_z=dnbinom(z,size=theta,prob=theta/(theta+mu))
          return(f_z)
          
        }
        
        #####  pmf truncated at z=36
        pmf=matrix(sapply(0:36,pmf_fn),ncol=1)
        
        const=sum(pmf)
        pmf=pmf/const
        
        y=seq(0,6,by=1/6)
        
        E_exp_alphaY=sum( exp(alpha*y)*pmf )
        
        E_Yexp_alphaY=sum( y*exp(alpha*y)*pmf )
        
        E_Y_past=E_Yexp_alphaY/E_exp_alphaY
        
        return(list(E_Y_past,E_exp_alphaY))
        
      }
      
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
    ##############  Generate simulated data ####################
    ############################################################
    # simulate data with a sample size of n=200 for treatment group
    {
      ############################################################
      #######   Fit models to the HAP2 data
      ############################################################
      
      ## generate count or continuous outcome
      ARC_data=readRDS("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_data.rds")
      
      End=830
      ARC_data2=filter(ARC_data,time <= End)
      
      ARC_formatted=formatting_fn(df=ARC_data2,last_day = End)
      
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
      AG_model_hv=coxph(Surv(Prev_time,time,Event) ~ Prev_outcome + strata(Visit_number),
                        id=elig_pid, 
                        data=ARC_survival_hv)
      
      
      ############  Estimated baseline intensities ############
      {
        AG_surv_hv=survfit(AG_model_hv,newdata=data.frame(Prev_outcome=0))
        strata=AG_surv_hv$strata
        
        v1=strata[1]
        cumhaz_v1=data.frame(time=AG_surv_hv$time[1:v1],cumhaz=AG_surv_hv$cumhaz[1:v1])
        base_intens_v1=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v1)
        
        
        v2=strata[2]
        cumhaz_v2=data.frame(time=AG_surv_hv$time[(v1+1):(v1+v2)],
                             cumhaz=AG_surv_hv$cumhaz[(v1+1):(v1+v2)])
        base_intens_v2=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v2)
        
        
        v3=strata[3]
        cumhaz_v3=data.frame(time=AG_surv_hv$time[(v1+v2+1):(v1+v2+v3)],
                             cumhaz=AG_surv_hv$cumhaz[(v1+v2+1):(v1+v2+v3)])
        base_intens_v3=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v3)
        
        
        v4=strata[4]
        cumhaz_v4=data.frame(time=AG_surv_hv$time[(v1+v2+v3+1):(v1+v2+v3+v4)],
                             cumhaz=AG_surv_hv$cumhaz[(v1+v2+v3+1):(v1+v2+v3+v4)])
        base_intens_v4=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v4)
      }
      
      
      #################    Outcome model - negative binomial model  #################
      {
        
        Outcome_model_hv <- glm.nb((6*Asthma_control)~Prev_outcome+time+Lag_time,
                                   data=ARC_visits_hv)
        
        theta_ARC <- Outcome_model_hv$theta
        
      }
      
      
      #############################################################
      ##############     Generate simulated data    ###############
      #############################################################
      {
        
        samp=sample(n_hv, size = N, replace=TRUE)
        
        ################################################
        Y0_sim=ARC_baseline_hv[samp, ]
        Y0_sim=data.frame(Y0_sim)%>%select(Asthma_control)%>%
          rename(Prev_outcome=Asthma_control)
        
        # There is no Visit_number in the dataset
        df1_sim=Times_gen_fn(lambda0_t=base_intens_v1,
                             Outcomes_vec=Y0_sim,
                             lambda_star=0.02,
                             cox_model = AG_model_hv)
        
        newdat_1=data.frame(Prev_outcome=Y0_sim,
                            time=df1_sim$time,
                            Lag_time=df1_sim$Lag_time)
        
        ############## The negative binomial model ##############
        {
          mu1=predict(Outcome_model_hv,newdata=newdat_1,type="response")
          
          Y1_sim=rep(NA,N)
          
          for(n in 1:N){
            Y1_sim[n]=Y_draw_fn(mu=mu1[n],theta=theta_ARC)
          }
          
          df1_sim=mutate(df1_sim,outcome=Y1_sim)
          
          for(i in 1:N){
            if(df1_sim$event[i]==0){ df1_sim$outcome[i]=NA }
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
          mu2=predict(Outcome_model_hv,newdata=newdat_2,type="response")
          Y2_sim=rep(NA,N)
          
          for(n in 1:N){
            Y2_sim[n]=Y_draw_fn(mu=mu2[n],theta=theta_ARC)
          }
          
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
          mu3=predict(Outcome_model_hv,newdata=newdat_3,type="response")
          
          Y3_sim=rep(NA,N)
          for(n in 1:N){
            Y3_sim[n]=Y_draw_fn(mu=mu3[n],theta=theta_ARC)
          }
          
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
          mu4=predict(Outcome_model_hv,newdata=newdat_4,type="response")
          
          Y4_sim=rep(NA,N)
          
          for(n in 1:N){
            Y4_sim[n]=Y_draw_fn(mu=mu4[n],theta=theta_ARC)
          }
          
          df4_sim=mutate(df4_sim,outcome=Y4_sim)
          
          for(i in 1:N){
            if(df4_sim$event[i]==0){ df4_sim$outcome[i]=NA }
          }
        }
        
        
        ################################################
        
        df_sim_full=rbind(df1_sim,df2_sim,df3_sim,df4_sim)%>%arrange(id,time)
        
        # df_sim=filter(df_sim_full, time > (prev_time+1))
      }
    }
    
    
    ############################################################
    ############## Calculate the True value of beta ############
    ############################################################
    {
      
      data_true <- f_transfer(df_sim_full)
      
      df_sim <- data_true
      All_visits_df_sim <- filter(df_sim, event == 1) %>% rename(Asthma_control = outcome)
      Visits_df_sim <- filter(All_visits_df_sim,time >= min(spline_seq),time <= max(spline_seq))
      theta_sim = Outcome_model_hv$theta
      coefs = Outcome_model_hv$coefficients
      a = coefs[3] + coefs[4]
      
      #####  the p x 401 matrix with columns B(t), 
      basis <- SplineBasis(knots = c(60,60,60,60,260,460,460,460,460))
      B_t_matrix <- evaluate(basis, spline_seq)
      
      # the length of beta
      p <- dim(B_t_matrix)[2]
      
      #####  we approximated the integral V=\int_t B(t)B(t)'dt using sums of rectangles of width 1
      V <- GramMatrix(basis)
      V_inverse <- solve(V)
      Weights_term2 <- V_inverse%*%t(B_t_matrix)
      
      # Term2_mat <- matrix(nrow = p, ncol = N)
      EYt_true <- matrix(0, nrow = length(spline_seq), ncol = 1)
      
      start <- Sys.time()
      for(i in 1:N){
        print(i)
        df_i1 <- filter(df_sim, id==i)
        
        ####  predictors for time = 1
        newdat_i <- data.frame(time=1, Lag_time=1, Prev_outcome=df_i1$Prev_outcome[1])
        
        df_i2 <- filter(Visits_df_sim, id==i)
        
        times_i <- df_i2$time+1
        
        if(length(times_i)>0){
          
          #####  predictors for times just after post-baseline visits in [a,b]
          newdat_i2 <- data.frame(time=times_i, Lag_time=1, Prev_outcome=df_i2$Asthma_control)
          
          newdat_i=rbind(newdat_i,newdat_i2)
          
        }
        
        #####  predicted means mu(t*) just after each of participant i's assessments
        mu_i <- predict.glm(Outcome_model_hv, 
                            newdata=newdat_i, 
                            type="response")
        
        times_i=c(1,times_i)
        
        time_means_fn=function(t){
          
          #####  for each t, mu(t) = mu(t*)exp(a Lag), where t* is the t* is the latest
          #####    day-after-an-assessment day on or before day t
          m=max(which(times_i<=t))
          s=t-times_i[m]
          
          temp=Cond_mean_fn(mu=mu_i[m]*exp(a*s),theta=theta_sim,alpha=alpha)
          
          return(temp[[1]])
          
        }
        
        Time_means <- matrix(sapply(spline_seq, time_means_fn), ncol=1)
        
        EYt_true <- EYt_true + Time_means / N 
        
      }
      end <- Sys.time()
      end - start
      
    }
    
    
    Beta_true_identity <- V_inverse %*% (1/2 * (t(B_t_matrix[-length(spline_seq), ]) %*% EYt_true[-length(spline_seq), ] +
                                                  t(B_t_matrix[-1, ]) %*% EYt_true[-1, ]))
    
    Beta_true_log <- V_inverse %*% (1/2 * (t(B_t_matrix[-length(spline_seq), ]) %*% log(EYt_true[-length(spline_seq), ]) +
                                             t(B_t_matrix[-1, ]) %*% log(EYt_true[-1, ])))
    
    B3 <- evaluate(basis, 90)
    B6 <- evaluate(basis, 180)
    B9 <- evaluate(basis, 270)
    B12 <- evaluate(basis, 360)
    
    month3_orig <- EYt_true[which(spline_seq == 90)]
    month6_orig <- EYt_true[which(spline_seq == 180)]
    month9_orig <- EYt_true[which(spline_seq == 270)]
    month12_orig <- EYt_true[which(spline_seq == 360)]
    
    month3_identity  <- B3 %*% Beta_true_identity
    month6_identity  <- B6 %*% Beta_true_identity
    month9_identity  <- B9 %*% Beta_true_identity
    month12_identity  <- B12 %*% Beta_true_identity
    
    month3_log  <- exp(B3 %*% Beta_true_log)
    month6_log  <- exp(B6 %*% Beta_true_log)
    month9_log  <- exp(B9 %*% Beta_true_log)
    month12_log  <- exp(B12 %*% Beta_true_log)
    
    
    return(list(alpha, 
                Beta_true_identity,
                Beta_true_log,
                c(month3_orig, month6_orig, month9_orig, month12_orig),
                c(month3_identity, month6_identity, month9_identity, month12_identity),
                c(month3_log, month6_log, month9_log, month12_log)
    ))
    
  }
  
  
  ### new
  {
    
    ### treatment group
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_log_neg")
    
    result_n12 <-  ture_value_ARC_data_neg(group = "treat", alpha = -1.2, seed  = 123, N = 2*10^5)
    saveRDS(result_n12, file = 'treat_true_neg12_200000')
    
    result_n08 <-  ture_value_ARC_data_neg(group = "treat", alpha = -0.8, seed  = 234, N = 2*10^5)
    saveRDS(result_n08, file = 'treat_true_neg08_200000')
    
    result_n04 <-  ture_value_ARC_data_neg(group = "treat", alpha = -0.4, seed  = 345, N = 2*10^5)
    saveRDS(result_n04, file = 'treat_true_neg04_200000')
    
    result_0 <-  ture_value_ARC_data_neg(group = "treat", alpha = 0, seed  = 456, N = 2*10^5)
    saveRDS(result_0, file = 'treat_true_0_200000')
    
    result_04 <-  ture_value_ARC_data_neg(group = "treat", alpha = 0.4, seed  = 567, N = 2*10^5)
    saveRDS(result_04, file = 'treat_true_04_200000')
    
    result_08 <-  ture_value_ARC_data_neg(group = "treat", alpha = 0.8, seed  = 678, N = 2*10^5)
    saveRDS(result_08, file = 'treat_true_08_200000')
    
    
    
    ### control group
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_log_neg")
    
    result_n12 <-  ture_value_ARC_data_neg(group = "control", alpha = -1.2, seed  = 123, N = 2*10^5)
    saveRDS(result_n12, file = 'control_true_neg12_200000')
    
    result_n08 <-  ture_value_ARC_data_neg(group = "control", alpha = -0.8, seed  = 234, N = 2*10^5)
    saveRDS(result_n08, file = 'control_true_neg08_200000')
    
    result_n04 <-  ture_value_ARC_data_neg(group = "control", alpha = -0.4, seed  = 345, N = 2*10^5)
    saveRDS(result_n04, file = 'control_true_neg04_200000')
    
    result_0 <-  ture_value_ARC_data_neg(group = "control", alpha = 0, seed  = 456, N = 2*10^5)
    saveRDS(result_0, file = 'control_true_0_200000')
    
    result_04 <-  ture_value_ARC_data_neg(group = "control", alpha = 0.4, seed  = 567, N = 2*10^5)
    saveRDS(result_04, file = 'control_true_04_200000')
    
    result_08 <-  ture_value_ARC_data_neg(group = "control", alpha = 0.8, seed  = 678, N = 2*10^5)
    saveRDS(result_08, file = 'control_true_08_200000')
    
    
    ##### 
    
    
    #####
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_log_neg")
    
    trt_neg12 <- readRDS(file = 'treat_true_neg12_200000') 
    trt_neg08 <- readRDS(file = 'treat_true_neg08_200000') 
    trt_neg04 <- readRDS(file = 'treat_true_neg04_200000') 
    trt_0 <- readRDS(file = 'treat_true_0_200000') 
    trt_04 <- readRDS(file = 'treat_true_04_200000')
    trt_08 <- readRDS(file = 'treat_true_08_200000')
    
    crl_neg12 <- readRDS(file = 'control_true_neg12_200000') 
    crl_neg08 <- readRDS(file = 'control_true_neg08_200000') 
    crl_neg04 <- readRDS(file = 'control_true_neg04_200000') 
    crl_0 <- readRDS(file = 'control_true_0_200000') 
    crl_04 <- readRDS(file = 'control_true_04_200000')
    crl_08 <- readRDS(file = 'control_true_08_200000')
    
    
    ## log link
    temp1 <- trt_neg12[[6]] |> round(5) |> t()
    temp2 <- trt_neg08[[6]] |> round(5) |> t()
    temp3 <- trt_neg04[[6]] |> round(5) |> t()
    temp4 <- trt_0[[6]] |> round(5) |> t()
    temp5 <- trt_04[[6]] |> round(5) |> t()
    temp6 <- trt_08[[6]] |> round(5) |> t()
    
    
    ## Below is the single index model
    # month <- 3
    c(temp1[1], temp2[1], temp3[1], temp4[1], temp5[1], temp6[1])
    month3_treat <- c(1.16734, 1.39203, 1.69211, 2.09621, 2.61513, 3.22004)
    
    # month <- 6
    c(temp1[2], temp2[2], temp3[2], temp4[2], temp5[2], temp6[2])
    month6_treat <- c(1.13348, 1.34842, 1.63375, 2.01638, 2.50784, 3.09132)

    # month <- 9
    c(temp1[3], temp2[3], temp3[3], temp4[3], temp5[3], temp6[3])
    month9_treat <-  c(1.11376, 1.32203, 1.59900, 1.96872, 2.44638, 3.01779)

    # month <- 12
    c(temp1[4], temp2[4], temp3[4], temp4[4], temp5[4], temp6[4])
    month12_treat <- c(1.09788, 1.30106, 1.57155, 1.93054, 2.39640, 2.95888)
    
    
    temp1 <- crl_neg12[[6]] |> round(5) |> t()
    temp2 <- crl_neg08[[6]] |> round(5) |> t()
    temp3 <- crl_neg04[[6]] |> round(5) |> t()
    temp4 <- crl_0[[6]] |> round(5) |> t()
    temp5 <- crl_04[[6]] |> round(5) |> t()
    temp6 <- crl_08[[6]] |> round(5) |> t()
    
    c(temp1[1], temp2[1], temp3[1], temp4[1], temp5[1], temp6[1])
    month3_control <- c(1.02777, 1.27702, 1.63975, 2.16540, 2.86166, 3.63791)
    
    c(temp1[2], temp2[2], temp3[2], temp4[2], temp5[2], temp6[2])
    month6_control <- c(1.00005, 1.23851, 1.58236, 2.08165, 2.74729, 3.51121)
    
    c(temp1[3], temp2[3], temp3[3], temp4[3], temp5[3], temp6[3])
    month9_control <- c(0.98010, 1.21039, 1.54180, 2.02261, 2.67015, 3.42586)
    
    c(temp1[4], temp2[4], temp3[4], temp4[4], temp5[4], temp6[4])
    month12_control <- c(0.96400, 1.18720, 1.50757, 1.97351, 2.60517, 3.35527)
    
    
    month3_treat <- c(1.16734, 1.39203, 1.69211, 2.09621, 2.61513, 3.22004)
    month6_treat <- c(1.13348, 1.34842, 1.63375, 2.01638, 2.50784, 3.09132)
    month9_treat <-  c(1.11376, 1.32203, 1.59900, 1.96872, 2.44638, 3.01779)
    month12_treat <- c(1.09788, 1.30106, 1.57155, 1.93054, 2.39640, 2.95888)
    
    month3_control <- c(1.02777, 1.27702, 1.63975, 2.16540, 2.86166, 3.63791)
    month6_control <- c(1.00005, 1.23851, 1.58236, 2.08165, 2.74729, 3.51121)
    month9_control <- c(0.98010, 1.21039, 1.54180, 2.02261, 2.67015, 3.42586)
    month12_control <- c(0.96400, 1.18720, 1.50757, 1.97351, 2.60517, 3.35527)
    
  }
  
  
}


###### 3. Single index model
{
  
  
  ture_value_ARC_data_single <- function(group = "control",
                                         type = "logit", # logit, log or identity
                                         alpha = -0.6,
                                         seed = 123,
                                         # N = 2*10^5
                                         N = 500,
                                         SIM_model = "fixed_coef", # "norm1coef", "fixed_coef", "fixed_band"
                                         SIM_tol = 1e-7
  ){  
    
    set.seed(seed)
    
    ####  period over which we are conducting inference
    spline_seq <- seq(60, 460, by=1)
    
    
    library(dplyr)
    library(survival)
    # library(MAVE)
    library(MASS)
    library(dfoptim)
    library(orthogonalsplinebasis) 
    library(BB)
    
    select <- dplyr::select
    
    
    #############################################################
    #################   Functions  ##############################
    #############################################################
    {
      
      formatting_fn=function(df, last_day){
        
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
      
      lambda0_fn=function(t,b,surv){
        
        times=surv$time
        cumhaz=surv$cumhaz
        
        increment=rep(0,length(times))
        increment[1]=cumhaz[1]
        for(j in 2:length(times)){ increment[j]=cumhaz[j]-cumhaz[j-1] }
        
        f_t=(1/b)*sum( 0.75*(1-( (t-times)/b )^2 )*( abs(t-times)<b )*increment)
        return(f_t)
        
      }
      
      Times_gen_fn=function(Start_vec=rep(0,N),
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
      
      
      NW <- function(Xb, Y, xb, y, h, kernel = "dnorm"){
        
        if(kernel == "dnorm"){
          K <- function(x, h){dnorm(x/h, 0, 1)} # Gaussian 
        } else if(kernel == "K2_Biweight"){
          K <- function(x, h){15/16*(1-(x/h)^2)^2 * (abs(x) <= h)} # K2_biweight
        } else if(kernel=="K4_Biweight"){
          K <- function(x, h){105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) }# K4_biweight
        }
        
        Kxb <- sapply(xb, function(x, Xb) K(Xb-x, h), Xb=Xb)
        
        Ylty <- sapply(y, function(x, Y) 1*(Y <= x), Y=Y)
        
        denom <- colSums(Kxb)
        
        fyxb <- (denom!=0)*crossprod(Kxb, Ylty)/(denom + (denom==0))
        
        return(fyxb)
        
      }
      
      Y_draw_fn_single <- function(Fhat, y){
        
        # density function
        Fhat1 <- c(0, Fhat[1:(length(Fhat) - 1)]) 
        pmf <- Fhat - Fhat1
        
        pmf[which(pmf < 1e-10)] <- 0 
        pmf <- pmf / sum(pmf)
        
        Y_new <- sample(y, size = 1, replace = T, prob = pmf)
        return(Y_new)
        
      }
      
      # The general function in the package
      Cond_mean_fn_single <- function(alpha, X, Y, x, beta, bandwidth){
        
        y <- sort(unique(Y))
        # conditional distribution
        Fhat <- NW(Xb = X %*% beta, Y = Y,
                   xb = x %*% beta, y = y,
                   h = bandwidth)
        
        # density function
        Fhat1 <- c(0, Fhat[1:(length(y) - 1)]) 
        pmf <- Fhat - Fhat1
        
        E_exp_alphaY <- sum( exp(alpha*y)*pmf )
        
        E_Yexp_alphaY <- sum( y*exp(alpha*y)*pmf )
        
        E_Y_past <- E_Yexp_alphaY/E_exp_alphaY
        
        return(list(E_Y_past, E_exp_alphaY))
        
      }
      
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
    ##############  Generate simulated data ####################
    ############################################################
    # simulate data with a sample size of n=200 for treatment group
    {
      ############################################################
      #######   Fit models to the HAP2 data
      ############################################################
      
      if(type == "logit"){
        ## generate binary outcome
        ARC_data <- readRDS("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_data_binary_15.rds")
        
        s_inverse <- function(mu){
          return(log(mu/(1-mu)))
        }
        
        s_fun <- function(z){
          return(1 / (1 + exp(-z)))
        }
        
      }else if(type == "log"){
        ## generate count or continuous outcome
        ARC_data = readRDS("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_data.rds")
        
        s_inverse <- function(mu){
          return(log(mu))
        }
        
        s_fun <- function(z){
          return(exp(z))
        }
        
      }else{
        ## generate count or continuous outcome
        ARC_data = readRDS("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_data.rds")
        
        s_inverse <- function(mu){
          return(mu)
        }
        
        s_fun <- function(z){
          return(z)
        }
      }
      
      y_unique <- unique(ARC_data$Asthma_control)
      y_unique <- y_unique[is.na(y_unique) == F] |> sort()
      
      End=830
      ARC_data2=filter(ARC_data,time <= End)
      
      ARC_formatted=formatting_fn(df=ARC_data2,last_day = End)
      
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
      AG_model_hv=coxph(Surv(Prev_time,time,Event) ~ Prev_outcome + strata(Visit_number),
                        id=elig_pid, 
                        data=ARC_survival_hv)
      
      
      ############  Estimated baseline intensities ############
      {
        AG_surv_hv=survfit(AG_model_hv,newdata=data.frame(Prev_outcome=0))
        strata=AG_surv_hv$strata
        
        v1=strata[1]
        cumhaz_v1=data.frame(time=AG_surv_hv$time[1:v1],cumhaz=AG_surv_hv$cumhaz[1:v1])
        base_intens_v1=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v1)
        
        
        v2=strata[2]
        cumhaz_v2=data.frame(time=AG_surv_hv$time[(v1+1):(v1+v2)],
                             cumhaz=AG_surv_hv$cumhaz[(v1+1):(v1+v2)])
        base_intens_v2=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v2)
        
        
        v3=strata[3]
        cumhaz_v3=data.frame(time=AG_surv_hv$time[(v1+v2+1):(v1+v2+v3)],
                             cumhaz=AG_surv_hv$cumhaz[(v1+v2+1):(v1+v2+v3)])
        base_intens_v3=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v3)
        
        
        v4=strata[4]
        cumhaz_v4=data.frame(time=AG_surv_hv$time[(v1+v2+v3+1):(v1+v2+v3+v4)],
                             cumhaz=AG_surv_hv$cumhaz[(v1+v2+v3+1):(v1+v2+v3+v4)])
        base_intens_v4=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v4)
      }
      
      
      #################    Outcome model - negative binomial model  #################
      {
        
        Xi <- data.frame(ARC_visits_hv$Prev_outcome,
                         ARC_visits_hv$time,
                         ARC_visits_hv$Lag_time)
        Xi <- as.matrix(Xi)
        Yi <- as.matrix(ARC_visits_hv$Asthma_control, ncol = 1)
        
        
        
        if(SIM_model == "norm1coef"){
          res1 <- SensIAT::fit_SensIAT_single_index_norm1coef_model(formula = Asthma_control ~ -1 + Prev_outcome + time + Lag_time,
                                                                    data = ARC_visits_hv,
                                                                    id = elig_pid, 
                                                                    kernel = "dnorm",
                                                                    mave.method = "meanMAVE",
                                                                    bw.selection = "mse",
                                                                    bw.method = "optim",
                                                                    bw.range = c(0.01, 1.5),
                                                                    reestimate.coef = 0)
          SIM_coef <- res1$coefficients
          SIM_band <- res1$bandwidth
          
        }else if(SIM_model ==  "fixed_coef"){
          
          res2 <- SensIAT::fit_SensIAT_single_index_fixed_coef_model(formula = Asthma_control ~ -1 + Prev_outcome + time + Lag_time,
                                                                     data = ARC_visits_hv,
                                                                     id = elig_pid, 
                                                                     kernel = "dnorm",
                                                                     method = "nmk",
                                                                     abs.tol = 1e-7,
                                                                     initial = NULL)
          SIM_coef <- res2$coef
          SIM_band <- res2$bandwidth
          
        }else{
          
          res3 <- SensIAT::fit_SensIAT_single_index_fixed_bandwidth_model(formula = Asthma_control ~ -1 + Prev_outcome + time + Lag_time,
                                                                          data = ARC_visits_hv,
                                                                          id = elig_pid, 
                                                                          kernel = "dnorm",
                                                                          method = "nmk",
                                                                          abs.tol = SIM_tol,
                                                                          initial = NULL)
          SIM_coef <- res3$coefficients
          SIM_band <- res3$bandwidth
          
        }
        
        Xb = Xi %*% SIM_coef
        
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
        
        df1_sim=Times_gen_fn(lambda0_t=base_intens_v1,
                             Outcomes_vec=Y0_sim,
                             lambda_star=0.02,
                             cox_model = AG_model_hv)
        
        newdat_1=data.frame(Prev_outcome=Y0_sim,
                            time=df1_sim$time,
                            Lag_time=df1_sim$Lag_time)
        
        ############## The single index model ##############
        {
          
          # generate the probability of Yi = 1 at the first assessment time
          Y1_sim <- rep(NA, N)
          for(n_index in 1:N){
            temp <- NW(Xb = Xb, Y = Yi,
                       xb = as.matrix(newdat_1[n_index, ]) %*% SIM_coef, y = y_unique,
                       h = SIM_band)
            Y1_sim[n_index] <- Y_draw_fn_single(Fhat = temp, y = y_unique)
          }
          
          df1_sim = mutate(df1_sim,outcome=Y1_sim)
          for(i in 1:N){
            if(df1_sim$event[i]==0){ df1_sim$outcome[i]=NA }
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
        
        ############## The single index model ##############
        {
          # generate the probability of Yi = 1 at the first assessment time
          Y2_sim <- rep(NA, N)
          for(n_index in 1:N){
            temp <- NW(Xb = Xb, Y = Yi,
                       xb = as.matrix(newdat_2[n_index, ]) %*% SIM_coef, y = y_unique,
                       h = SIM_band)
            Y2_sim[n_index] <- Y_draw_fn_single(Fhat = temp, y = y_unique)
          }
          
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
        
        ############## The single index model ##############
        {
          # generate the probability of Yi = 1 at the first assessment time
          Y3_sim <- rep(NA, N)
          for(n_index in 1:N){
            temp <- NW(Xb = Xb, Y = Yi,
                       xb = as.matrix(newdat_3[n_index, ]) %*% SIM_coef, y = y_unique,
                       h = SIM_band)
            Y3_sim[n_index] <- Y_draw_fn_single(Fhat = temp, y = y_unique)
          }
          
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
        
        ############## The single index model ##############
        {
          # generate the probability of Yi = 1 at the first assessment time
          Y4_sim <- rep(NA, N)
          for(n_index in 1:N){
            temp <- NW(Xb = Xb, Y = Yi,
                       xb = as.matrix(newdat_4[n_index, ]) %*% SIM_coef, y = y_unique,
                       h = SIM_band)
            Y4_sim[n_index] <- Y_draw_fn_single(Fhat = temp, y = y_unique)
          }
          
          df4_sim=mutate(df4_sim,outcome=Y4_sim)
          for(i in 1:N){
            if(df4_sim$event[i]==0){ df4_sim$outcome[i]=NA }
          }
        }
        
        ################################################
        
        df_sim_full <- rbind(df1_sim,df2_sim,df3_sim,df4_sim) %>% arrange(id,time)
        
      }
      
      
    }
    
    
    ############################################################
    ############## Calculate the True value of beta ############
    ############################################################
    {
      
      df_sim <- f_transfer(df_sim_full) 
      
      All_visits_df_sim <- filter(df_sim, event==1) %>% rename(Asthma_control=outcome)
      Visits_df_sim <- filter(All_visits_df_sim,
                              time >= min(spline_seq),
                              time <= max(spline_seq))
      
      f_data_impute <- function(t, 
                                df_i1,
                                min_time,
                                max_time,
                                df_i1_outcome_max,
                                df_i1_Prev_outcome_1){
        
        df_est <- c(NA, t, NA)
        
        if(t %in% df_i1$time){
          index <- which(df_i1$time == t)
          df_est <- c(df_i1$Prev_outcome[index], t, df_i1$Lag_time[index])
        }
        else if(t < min_time){
          df_est <- c(df_i1_Prev_outcome_1, t, t)
        }
        else if(t <= max_time){
          t_index1 <- max(which(df_i1$time < t))
          temp <- df_i1$outcome[t_index1]
          df_est <- c(temp, t, t - df_i1$time[t_index1])
        }
        else{
          temp <- df_i1_outcome_max
          df_est <- c(temp, t, t - max_time) 
        }
        
        return(df_est)
        
      }
      
      # EYt_true <- rep(0, 4)
      # Time_all <- c(3, 6, 9, 12) * 30
      # length_time <- 4
      
      # start <- Sys.time()
      # for(i in 1:N){
      #   
      #   print(i)
      #   df_i1 <- filter(df_sim, id==i)
      #   
      #   min_time <- min(df_i1$time)
      #   max_time <- max(df_i1$time)
      #   df_i1_outcome_max <- df_i1$outcome[nrow(df_i1)]
      #   df_i1_Prev_outcome_1 <- df_i1$Prev_outcome[1]
      # 
      #   temp <- lapply(Time_all, 
      #                  f_data_impute, 
      #                  df_i1 = df_i1,
      #                  min_time = min_time,
      #                  max_time = max_time,
      #                  df_i1_outcome_max = df_i1_outcome_max,
      #                  df_i1_Prev_outcome_1 = df_i1_Prev_outcome_1)
      #   df_est <- as.data.frame(do.call(rbind, temp))
      #   
      #   Time_means_single  <- matrix(0, nrow = length_time, ncol = 1)
      #   for(l in 1:length_time){
      #     temp <- Cond_mean_fn_single(alpha,
      #                                 X = Xi,
      #                                 Y = Yi,
      #                                 x = as.matrix(df_est[l, ]),
      #                                 beta = SIM_coef,
      #                                 bandwidth = SIM_band)
      #     Time_means_single[l] <- temp[[1]]
      #   }
      #   EYt_true <- EYt_true + Time_means_single / N 
      #   
      # }
      # end <- Sys.time()
      # end - start
      
      
      #####  the p x 401 matrix with columns B(t), 
      basis <- SplineBasis(knots = c(60,60,60,60,260,460,460,460,460))
      B_t_matrix <- evaluate(basis, spline_seq)
      p <- dim(B_t_matrix)[2]
      
      #####  we approximated the integral V=\int_t B(t)B(t)'dt using sums of rectangles of width 1
      V <- GramMatrix(basis)
      V_inverse <- solve(V)
      Weights_term2 <- V_inverse%*%t(B_t_matrix)
      
      length_time <- length(spline_seq)
      EYt_true <- matrix(0, nrow = length_time, ncol = 1)
      Time_all <- spline_seq
      
      start <- Sys.time()
      for(i in 1:N){
        
        print(i)
        df_i1 <- filter(df_sim, id==i)
        
        min_time <- min(df_i1$time)
        max_time <- max(df_i1$time)
        df_i1_outcome_max <- df_i1$outcome[nrow(df_i1)]
        df_i1_Prev_outcome_1 <- df_i1$Prev_outcome[1]
        
        temp <- lapply(Time_all,
                       f_data_impute,
                       df_i1 = df_i1,
                       min_time = min_time,
                       max_time = max_time,
                       df_i1_outcome_max = df_i1_outcome_max,
                       df_i1_Prev_outcome_1 = df_i1_Prev_outcome_1)
        df_est <- as.data.frame(do.call(rbind, temp))
        
        Time_means_single  <- matrix(0, nrow = length_time, ncol = 1)
        for(l in 1:length_time){
          temp <- Cond_mean_fn_single(alpha,
                                      X = Xi,
                                      Y = Yi,
                                      x = as.matrix(df_est[l, ]),
                                      beta = SIM_coef,
                                      bandwidth = SIM_band)
          Time_means_single[l] <- temp[[1]]
        }
        
        EYt_true <- EYt_true + Time_means_single / N
        
      }
      end <- Sys.time()
      end - start
      
    }
    
    
    
    Beta_true <- V_inverse %*% (1/2 * (t(B_t_matrix[-length(spline_seq), ]) %*% s_inverse(EYt_true[-length(spline_seq), ]) +
                                         t(B_t_matrix[-1, ]) %*% s_inverse(EYt_true[-1, ])))
    
    month3_orig <- EYt_true[which(spline_seq == 90)]
    month6_orig <- EYt_true[which(spline_seq == 180)]
    month9_orig <- EYt_true[which(spline_seq == 270)]
    month12_orig <- EYt_true[which(spline_seq == 360)]
    
    
    B3 <- evaluate(basis, 90)
    B6 <- evaluate(basis, 180)
    B9 <- evaluate(basis, 270)
    B12 <- evaluate(basis, 360)
    
    month3_link  <- s_fun(B3 %*% Beta_true)
    month6_link  <- s_fun(B6 %*% Beta_true)
    month9_link  <- s_fun(B9 %*% Beta_true)
    month12_link  <- s_fun(B12 %*% Beta_true)
    
    
    return(list(alpha, 
                Beta_true,
                c(month3_orig, month6_orig, month9_orig, month12_orig),
                c(month3_link, month6_link, month9_link, month12_link)
    ))
    
  }
  
  
  
  ### logit
  {
    
    ### treatment group
    
    {
      
      setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_logit")
      
      result_n08 <-  ture_value_ARC_data_single(group = "treat",
                                                type = "logit", # logit, log or identity
                                                alpha = -0.8, 
                                                seed  = 123,
                                                N = 2*10^5)
      # N = 200)
      saveRDS(result_n08, file = 'treat_true_neg08_200000')
      
      result_n04 <-  ture_value_ARC_data_single(group = "treat",
                                                type = "logit", # logit, log or identity
                                                alpha = -0.4, 
                                                seed  = 234,
                                                N = 2*10^5)
      saveRDS(result_n04, file = 'treat_true_neg04_200000')
      
      result_n0 <-  ture_value_ARC_data_single(group = "treat",
                                               type = "logit", # logit, log or identity
                                               alpha = 0,
                                               seed  = 345,
                                               N = 2*10^5)
      saveRDS(result_n0, file = 'treat_true_0_200000')
      
      result_04 <-  ture_value_ARC_data_single(group = "treat",
                                               type = "logit", # logit, log or identity
                                               alpha = 0.4,
                                               seed  = 456,
                                               N = 2*10^5)
      saveRDS(result_04, file = 'treat_true_04_200000')
      
      result_08 <-  ture_value_ARC_data_single(group = "treat", 
                                               type = "logit", # logit, log or identity
                                               alpha = 0.8, 
                                               seed  = 567,
                                               N = 2*10^5)
      saveRDS(result_08, file = 'treat_true_08_200000')
      
      result_12 <-  ture_value_ARC_data_single(group = "treat",
                                               type = "logit", # logit, log or identity
                                               alpha = 1.2, 
                                               seed  = 678,
                                               N = 2*10^5)
      saveRDS(result_12, file = 'treat_true_12_200000')
      
      
    }
    
    
    ### control group
    
    {
      
      setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_logit")
      
      result_n08 <-  ture_value_ARC_data_single(group = "control",
                                                type = "logit", # logit, log or identity
                                                alpha = -0.8, 
                                                seed  = 123,
                                                N = 2*10^5)
      saveRDS(result_n08, file = 'control_true_neg08_200000')
      
      result_n04 <-  ture_value_ARC_data_single(group = "control",
                                                type = "logit", # logit, log or identity
                                                alpha = -0.4, 
                                                seed  = 234,
                                                N = 2*10^5)
      saveRDS(result_n04, file = 'control_true_neg04_200000')
      
      result_n0 <-  ture_value_ARC_data_single(group = "control",
                                               type = "logit", # logit, log or identity
                                               alpha = 0,
                                               seed  = 345,
                                               N = 2*10^5)
      saveRDS(result_n0, file = 'control_true_0_200000')
      
      result_04 <-  ture_value_ARC_data_single(group = "control",
                                               type = "logit", # logit, log or identity
                                               alpha = 0.4,
                                               seed  = 456,
                                               N = 2*10^5)
      saveRDS(result_04, file = 'control_true_04_200000')
      
      result_08 <-  ture_value_ARC_data_single(group = "control", 
                                               type = "logit", # logit, log or identity
                                               alpha = 0.8, 
                                               seed  = 567,
                                               N = 2*10^5)
      saveRDS(result_08, file = 'control_true_08_200000')
      
      result_12 <-  ture_value_ARC_data_single(group = "control",
                                               type = "logit", # logit, log or identity
                                               alpha = 1.2, 
                                               seed  = 678,
                                               N = 2*10^5)
      saveRDS(result_12, file = 'control_true_12_200000')
      
      
    }
    
    
    ##### 
    
    
    #####
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_logit_single")
    
    crl_neg08 <- readRDS(file = 'control_true_neg08_200000') 
    crl_neg04 <- readRDS(file = 'control_true_neg04_200000') 
    crl_0 <- readRDS(file = 'control_true_0_200000') 
    crl_04 <- readRDS(file = 'control_true_04_200000')
    crl_08 <- readRDS(file = 'control_true_08_200000')
    crl_12 <- readRDS(file = 'control_true_12_200000')
    
    trt_neg08 <- readRDS(file = 'treat_true_neg08_200000') 
    trt_neg04 <- readRDS(file = 'treat_true_neg04_200000') 
    trt_0 <- readRDS(file = 'treat_true_0_200000') 
    trt_04 <- readRDS(file = 'treat_true_04_200000')
    trt_08 <- readRDS(file = 'treat_true_08_200000')
    trt_12 <- readRDS(file = 'treat_true_12_200000')
    
    
    # # original   
    trt_neg08[[3]] |> round(3) |> t()
    trt_neg04[[3]] |> round(3) |> t()
    trt_0[[3]] |> round(3) |> t()
    trt_04[[3]] |> round(3) |> t()
    trt_08[[3]] |> round(3) |> t()
    trt_12[[3]] |> round(3) |> t()
    
    crl_neg08[[3]] |> round(3) |> t()
    crl_neg04[[3]] |> round(3) |> t()
    crl_0[[3]] |> round(3) |> t()
    crl_04[[3]] |> round(3) |> t()
    crl_08[[3]] |> round(3) |> t()
    crl_12[[3]] |> round(3) |> t()
    
    # # month <- 3
    # month3_treat <- c(0.171, 0.228, 0.297, 0.375, 0.462, 0.552)
    # month3_control <- c(0.152, 0.202, 0.262, 0.333, 0.412, 0.497)
    # 
    # # month <- 6
    # month6_treat <- c(0.186, 0.245, 0.316, 0.394, 0.481, 0.57)
    # month6_control <- c(0.158, 0.209, 0.271, 0.343, 0.422, 0.506)
    # 
    # # month <- 9
    # month9_treat <- c(0.195, 0.256, 0.328, 0.407, 0.493, 0.58)
    # month9_control <- c(0.162, 0.214, 0.277, 0.348, 0.428, 0.513)
    # 
    # # month <- 12
    # month12_treat <- c(0.2, 0.263, 0.334, 0.414, 0.499, 0.586)
    # month12_control <- c(0.164, 0.216, 0.28, 0.352, 0.433, 0.517)
    
    
    # # link 
    
    trt_neg08[[4]] |> round(3) |> t()
    trt_neg04[[4]] |> round(3) |> t()
    trt_0[[4]] |> round(3) |> t()
    trt_04[[4]] |> round(3) |> t()
    trt_08[[4]] |> round(3) |> t()
    trt_12[[4]] |> round(3) |> t()
    
    crl_neg08[[4]] |> round(3) |> t()
    crl_neg04[[4]] |> round(3) |> t()
    crl_0[[4]] |> round(3) |> t()
    crl_04[[4]] |> round(3) |> t()
    crl_08[[4]] |> round(3) |> t()
    crl_12[[4]] |> round(3) |> t()
    
    # month <- 3
    month3_treat <- c(0.172, 0.229, 0.298, 0.375, 0.462, 0.553)
    month3_control <- c(0.152, 0.203, 0.263, 0.334, 0.413, 0.498)
    
    # month <- 6
    month6_treat <- c(0.187, 0.247, 0.318, 0.396, 0.482, 0.571)
    month6_control <- c(0.158, 0.21, 0.271, 0.343, 0.422, 0.507)
    
    # month <- 9
    month9_treat <- c(0.195, 0.256, 0.328, 0.406, 0.493, 0.58)
    month9_control <- c(0.162, 0.214, 0.277, 0.348, 0.428, 0.513)
    
    # month <- 12
    month12_treat <- c(0.201, 0.263, 0.334, 0.414, 0.5, 0.586)
    month12_control <- c(0.163, 0.217, 0.28, 0.352, 0.433, 0.518)
    
    
    
  }
  
  
  ### log
  {
    
    ### treatment group
    
    {
      
      setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_log")
      
      result_n12 <-  ture_value_ARC_data_single(group = "treat",
                                                type = "log", # logit, log or identity
                                                alpha = -1.2, 
                                                seed  = 678,
                                                N = 2*10^5)
      saveRDS(result_n12, file = 'treat_true_neg12_200000')
      
      result_n08 <-  ture_value_ARC_data_single(group = "treat",
                                                type = "log", # logit, log or identity
                                                alpha = -0.8, 
                                                seed  = 123,
                                                N = 2*10^5)
      saveRDS(result_n08, file = 'treat_true_neg08_200000')
      
      result_n04 <-  ture_value_ARC_data_single(group = "treat",
                                                type = "log", # logit, log or identity
                                                alpha = -0.4, 
                                                seed  = 234,
                                                N = 2*10^5)
      saveRDS(result_n04, file = 'treat_true_neg04_200000')
      
      result_n0 <-  ture_value_ARC_data_single(group = "treat",
                                               type = "log", # logit, log or identity
                                               alpha = 0,
                                               seed  = 345,
                                               N = 2*10^5)
      saveRDS(result_n0, file = 'treat_true_0_200000')
      
      result_04 <-  ture_value_ARC_data_single(group = "treat",
                                               type = "log", # logit, log or identity
                                               alpha = 0.4,
                                               seed  = 456,
                                               N = 2*10^5)
      saveRDS(result_04, file = 'treat_true_04_200000')
      
      result_08 <-  ture_value_ARC_data_single(group = "treat", 
                                               type = "log", # logit, log or identity
                                               alpha = 0.8, 
                                               seed  = 567,
                                               N = 2*10^5)
      saveRDS(result_08, file = 'treat_true_08_200000')
      
      
    }
    
    
    ### control group
    
    {
      
      setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_log")
      
      result_n12 <-  ture_value_ARC_data_single(group = "control",
                                                type = "log", # logit, log or identity
                                                alpha = -1.2, 
                                                seed  = 678,
                                                N = 2*10^5)
      saveRDS(result_n12, file = 'control_true_neg12_200000')
      
      result_n08 <-  ture_value_ARC_data_single(group = "control",
                                                type = "log", # logit, log or identity
                                                alpha = -0.8, 
                                                seed  = 123,
                                                N = 2*10^5)
      saveRDS(result_n08, file = 'control_true_neg08_200000')
      
      result_n04 <-  ture_value_ARC_data_single(group = "control",
                                                type = "log", # logit, log or identity
                                                alpha = -0.4, 
                                                seed  = 234,
                                                N = 2*10^5)
      saveRDS(result_n04, file = 'control_true_neg04_200000')
      
      result_n0 <-  ture_value_ARC_data_single(group = "control",
                                               type = "log", # logit, log or identity
                                               alpha = 0,
                                               seed  = 345,
                                               N = 2*10^5)
      saveRDS(result_n0, file = 'control_true_0_200000')
      
      result_04 <-  ture_value_ARC_data_single(group = "control",
                                               type = "log", # logit, log or identity
                                               alpha = 0.4,
                                               seed  = 456,
                                               N = 2*10^5)
      saveRDS(result_04, file = 'control_true_04_200000')
      
      result_08 <-  ture_value_ARC_data_single(group = "control", 
                                               type = "log", # logit, log or identity
                                               alpha = 0.8, 
                                               seed  = 567,
                                               N = 2*10^5)
      saveRDS(result_08, file = 'control_true_08_200000')
      
      
    }
    
    
    ##### 
    
    
    #####
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_log_single")
    
    
    crl_neg12 <- readRDS(file = 'control_true_neg12_200000')
    crl_neg08 <- readRDS(file = 'control_true_neg08_200000') 
    crl_neg04 <- readRDS(file = 'control_true_neg04_200000') 
    crl_0 <- readRDS(file = 'control_true_0_200000') 
    crl_04 <- readRDS(file = 'control_true_04_200000')
    crl_08 <- readRDS(file = 'control_true_08_200000')
    
    trt_neg12 <- readRDS(file = 'treat_true_neg12_200000')
    trt_neg08 <- readRDS(file = 'treat_true_neg08_200000') 
    trt_neg04 <- readRDS(file = 'treat_true_neg04_200000') 
    trt_0 <- readRDS(file = 'treat_true_0_200000') 
    trt_04 <- readRDS(file = 'treat_true_04_200000')
    trt_08 <- readRDS(file = 'treat_true_08_200000')
    
    
    trt_neg12[[3]] |> round(3) |> t()
    trt_neg08[[3]] |> round(3) |> t()
    trt_neg04[[3]] |> round(3) |> t()
    trt_0[[3]] |> round(3) |> t()
    trt_04[[3]] |> round(3) |> t()
    trt_08[[3]] |> round(3) |> t()
    
    
    crl_neg12[[3]] |> round(3) |> t()
    crl_neg08[[3]] |> round(3) |> t()
    crl_neg04[[3]] |> round(3) |> t()
    crl_0[[3]] |> round(3) |> t()
    crl_04[[3]] |> round(3) |> t()
    crl_08[[3]] |> round(3) |> t()
    
    
    ## Below is the single index model + fixcoef
    # month <- 3
    month3_treat <- c(1.141, 1.434, 1.78, 2.191, 2.703, 3.341)
    month3_control <- c(0.94, 1.323, 1.799, 2.342, 2.92, 3.497)
    
    # month <- 6
    month6_treat <- c(1.113, 1.394, 1.728, 2.13, 2.636, 3.283)
    month6_control <- c(0.921, 1.297, 1.766, 2.304, 2.884, 3.472)
    
    # month <- 9
    month9_treat <- c(1.088, 1.359, 1.686, 2.076, 2.58, 3.23)
    month9_control <- c(0.905, 1.275, 1.737, 2.269, 2.85, 3.443)
    
    # month <- 12
    month12_treat <- c(1.069, 1.334, 1.654, 2.039, 2.541, 3.188)
    month12_control <- c(0.895, 1.26, 1.716, 2.245, 2.828, 3.422)
    
    
  }
  
  
  ### identity
  {
    
    ### treatment group
    
    {
      setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_identity")
      
      
      result_n12 <-  ture_value_ARC_data_single(group = "treat",
                                                type = "identity", # logit, log or identity
                                                alpha = -1.2, 
                                                seed  = 678,
                                                N = 2*10^5)
      saveRDS(result_n12, file = 'treat_true_neg12_200000')
      
      result_n08 <-  ture_value_ARC_data_single(group = "treat",
                                                type = "identity", # logit, log or identity
                                                alpha = -0.8, 
                                                seed  = 123,
                                                N = 2*10^5)
      saveRDS(result_n08, file = 'treat_true_neg08_200000')
      
      result_n04 <-  ture_value_ARC_data_single(group = "treat",
                                                type = "identity", # logit, log or identity
                                                alpha = -0.4, 
                                                seed  = 234,
                                                N = 2*10^5)
      saveRDS(result_n04, file = 'treat_true_neg04_200000')
      
      result_n0 <-  ture_value_ARC_data_single(group = "treat",
                                               type = "identity", # logit, log or identity
                                               alpha = 0,
                                               seed  = 345,
                                               N = 2*10^5)
      saveRDS(result_n0, file = 'treat_true_neg0_200000')
      
      result_04 <-  ture_value_ARC_data_single(group = "treat",
                                               type = "identity", # logit, log or identity
                                               alpha = 0.4,
                                               seed  = 456,
                                               N = 2*10^5)
      saveRDS(result_04, file = 'treat_true_04_200000')
      
      result_08 <-  ture_value_ARC_data_single(group = "treat", 
                                               type = "identity", # logit, log or identity
                                               alpha = 0.8, 
                                               seed  = 567,
                                               N = 2*10^5)
      saveRDS(result_08, file = 'treat_true_08_200000')
      
      
    }
    
    
    ### control group
    
    {
      
      setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_identity")
      
      result_n12 <-  ture_value_ARC_data_single(group = "control",
                                                type = "identity", # logit, log or identity
                                                alpha = -1.2, 
                                                seed  = 678,
                                                N = 2*10^5)
      saveRDS(result_n12, file = 'control_true_neg12_200000')
      
      result_n08 <-  ture_value_ARC_data_single(group = "control",
                                                type = "identity", # logit, log or identity
                                                alpha = -0.8, 
                                                seed  = 123,
                                                N = 2*10^5)
      saveRDS(result_n08, file = 'control_true_neg08_200000')
      
      result_n04 <-  ture_value_ARC_data_single(group = "control",
                                                type = "identity", # logit, log or identity
                                                alpha = -0.4, 
                                                seed  = 234,
                                                N = 2*10^5)
      saveRDS(result_n04, file = 'control_true_neg04_200000')
      
      result_n0 <-  ture_value_ARC_data_single(group = "control",
                                               type = "identity", # logit, log or identity
                                               alpha = 0,
                                               seed  = 345,
                                               N = 2*10^5)
      saveRDS(result_n0, file = 'control_true_neg0_200000')
      
      result_04 <-  ture_value_ARC_data_single(group = "control",
                                               type = "identity", # logit, log or identity
                                               alpha = 0.4,
                                               seed  = 456,
                                               N = 2*10^5)
      saveRDS(result_04, file = 'control_true_04_200000')
      
      result_08 <-  ture_value_ARC_data_single(group = "control", 
                                               type = "identity", # logit, log or identity
                                               alpha = 0.8, 
                                               seed  = 567,
                                               N = 2*10^5)
      saveRDS(result_08, file = 'control_true_08_200000')
      
      
    }
    
    
    ##### 
    
    
    #####
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/true_value/simu_logit")
    
    
    crl_neg12 <- readRDS(file = 'control_true_neg12_200000')
    crl_neg08 <- readRDS(file = 'control_true_neg08_200000') 
    crl_neg04 <- readRDS(file = 'control_true_neg04_200000') 
    crl_0 <- readRDS(file = 'control_true_0_200000') 
    crl_04 <- readRDS(file = 'control_true_04_200000')
    crl_08 <- readRDS(file = 'control_true_08_200000')
    
    trt_neg12 <- readRDS(file = 'treat_true_neg12_200000')
    trt_neg08 <- readRDS(file = 'treat_true_neg08_200000') 
    trt_neg04 <- readRDS(file = 'treat_true_neg04_200000') 
    trt_0 <- readRDS(file = 'treat_true_0_200000') 
    trt_04 <- readRDS(file = 'treat_true_04_200000')
    trt_08 <- readRDS(file = 'treat_true_08_200000')
    
    
    trt_neg12[[2]] |> round(3) |> t()
    trt_neg08[[2]] |> round(3) |> t()
    trt_neg04[[2]] |> round(3) |> t()
    trt_0[[2]] |> round(3) |> t()
    trt_04[[2]] |> round(3) |> t()
    trt_08[[2]] |> round(3) |> t()
    
    
    crl_neg12[[2]] |> round(3) |> t()
    crl_neg08[[2]] |> round(3) |> t()
    crl_neg04[[2]] |> round(3) |> t()
    crl_0[[2]] |> round(3) |> t()
    crl_04[[2]] |> round(3) |> t()
    crl_08[[2]] |> round(3) |> t()
    
    
    ## Below is the single index model - need to change
    # month <- 3
    month3_treat <- c(0.202, 0.247, 0.300, 0.357, 0.419)
    month3_control <- c(0.160, 0.199, 0.243, 0.294, 0.350)
    
    # month <- 6
    month6_treat <- c(0.218, 0.265, 0.320, 0.376, 0.439)
    month6_control <- c(0.170, 0.211, 0.257, 0.310, 0.367)
    
    # month <- 9
    month9_treat <- c(0.228, 0.276, 0.331, 0.389, 0.452)
    month9_control <- c(0.182, 0.225, 0.272, 0.326, 0.385)
    
    # month <- 12
    month12_treat <- c(0.233, 0.282, 0.337, 0.396, 0.458)
    month12_control <- c(0.196, 0.241, 0.290, 0.346, 0.405)
    
    
  }
  
  
}

