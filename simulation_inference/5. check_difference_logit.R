rm(list=ls())
library(dplyr)


### parametric bootstrap - logistic regression model
{
 
   
  # setting
  {
    
    # # month <- 3
    # month3_treat <- c(0.202, 0.247, 0.300, 0.357, 0.419)
    # month3_control <- c(0.160, 0.199, 0.243, 0.294, 0.350)
    # 
    # # month <- 6
    # month6_treat <- c(0.218, 0.265, 0.319, 0.376, 0.439)
    # month6_control <- c(0.171, 0.211, 0.257, 0.310, 0.367)
    # 
    # # month <- 9
    # month9_treat <- c(0.228, 0.276, 0.331, 0.389, 0.452)
    # month9_control <- c(0.182, 0.225, 0.272, 0.326, 0.385)
    # 
    # # month <- 12
    # month12_treat <- c(0.233, 0.282, 0.337, 0.396, 0.458)
    # month12_control <- c(0.196, 0.241, 0.290, 0.346, 0.405)
    
    
    # # single index model +  original
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
    
    # # # single index model + link
    # # month <- 3
    # month3_treat <- c(0.172, 0.229, 0.298, 0.375, 0.462, 0.553)
    # month3_control <- c(0.152, 0.203, 0.263, 0.334, 0.413, 0.498)
    # 
    # # month <- 6
    # month6_treat <- c(0.187, 0.247, 0.318, 0.396, 0.482, 0.571)
    # month6_control <- c(0.158, 0.21, 0.271, 0.343, 0.422, 0.507)
    # 
    # # month <- 9
    # month9_treat <- c(0.195, 0.256, 0.328, 0.406, 0.493, 0.58)
    # month9_control <- c(0.162, 0.214, 0.277, 0.348, 0.428, 0.513)
    # 
    # # month <- 12
    # month12_treat <- c(0.201, 0.263, 0.334, 0.414, 0.5, 0.586)
    # month12_control <- c(0.163, 0.217, 0.28, 0.352, 0.433, 0.518)
    
    # logistic + logit link
    month3_treat <- c(0.17528, 0.23253, 0.29979, 0.37766, 0.46339, 0.55269)
    month6_treat <- c(0.19083, 0.25123, 0.31990, 0.39857, 0.48426, 0.57133)
    month9_treat <- c(0.19837, 0.26035, 0.33039, 0.40918, 0.49422, 0.58054)
    month12_treat <- c(0.20356, 0.26597, 0.33736, 0.41639, 0.50118, 0.58660)

    month3_control <- c(0.13761, 0.18518, 0.24377, 0.31265, 0.39025, 0.47700)
    month6_control <- c(0.14696, 0.19688, 0.25781, 0.32844, 0.40780, 0.49500)
    month9_control <- c(0.15750, 0.20958, 0.27300, 0.34542, 0.42643, 0.51300)
    month12_control <- c(0.17020, 0.22499, 0.29071, 0.36504, 0.44775, 0.53400)
    
    
    Alpha_seq1 <- c(-0.8, -0.4, 0, 0.4, 0.8, 1.2)
    Alpha_seq2 <- c(-0.8, -0.4, 0, 0.4, 0.8, 1.2)
    
    N <- 500
    boot_size <- 200
    
    wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
      CI_left <- mean_orig - Zstar *sqrt(var_orig)
      CI_right <- mean_orig + Zstar *sqrt(var_orig)
      
      return(c(CI_left, CI_right))
    }
    
    atanh <- function(x){
      return(1/2 * log((1 + x)/(1 - x)))
    }
    tanh <- function(x){
      return((exp(x) - exp(-x)) / (exp(x) + exp(-x)))
    }
    
    f_diff <- function(alpha_trt, alpha_crl, month,
                       weight_type = "weight1",
                       true_value){
      
      Alpha_seq1 <- c(-0.8, -0.4, 0, 0.4, 0.8, 1.2)
      Alpha_seq2 <- c(-0.8, -0.4, 0, 0.4, 0.8, 1.2)
      index_trt <- which(Alpha_seq1 == alpha_trt)
      index_crl <- which(Alpha_seq2 == alpha_crl)
      
      index_month <- which(c(1, 3, 6, 9, 12) == month)
      
      # for original data's estimation
      Month_3_trt <- rep(NA, N)
      Var_3_trt <- rep(NA, N) 
      Month_3_crl <- rep(NA, N)
      Var_3_crl <- rep(NA, N) 
      
      Effect <-  rep(NA, N)
      Var_effect <-  rep(NA, N)
      C_month3 <- rep(NA, N)
      
      # for confidence interval
      wald_l_month3_EY <- rep(NA, N)
      wald_r_month3_EY <- rep(NA, N)
      wald_l_month3_Bt <- rep(NA, N)
      wald_r_month3_Bt <- rep(NA, N)
      
      paraboot_l_month3_EY_se <- rep(NA, N)
      paraboot_r_month3_EY_se <- rep(NA, N)
      paraboot_l_month3_Bt_se <- rep(NA, N)
      paraboot_r_month3_Bt_se <- rep(NA, N)
      
      # for coverage
      Covg_wald_EY <- rep(NA, N)
      Covg_wald_Bt <- rep(NA, N)
      Covg_paraboot_EY_se <- rep(NA, N)
      Covg_paraboot_Bt_se <- rep(NA, N)
      
      for(i in 1:N){
        
        Month_3_trt[i] <- orig_treat[[i]][[1]][index_trt, index_month]
        Var_3_trt[i] <- orig_treat[[i]][[2]][index_trt, index_month]
        Month_3_crl[i] <- orig_control[[i]][[1]][index_crl, index_month]
        Var_3_crl[i] <- orig_control[[i]][[2]][index_crl, index_month]
        
        Effect[i] <- Month_3_trt[i] - Month_3_crl[i]
        Var_effect[i] <- Var_3_trt[i] + Var_3_crl[i]
        
        # Wald confidence interval - original
        temp2 <- wald(Effect[i], Var_effect[i])
        wald_l_month3_EY[i]  <- temp2[1]
        wald_r_month3_EY[i]  <- temp2[2]
        
        # Wald confidence interval - trans
        temp3 <- wald(atanh(Effect[i]), Var_effect[i] / ( 1 - Effect[i]^2 )^2 )
        wald_l_month3_Bt[i]  <- tanh(temp3[1])
        wald_r_month3_Bt[i]  <- tanh(temp3[2])
        
        # coverage
        Covg_wald_EY[i] <- ifelse((wald_l_month3_EY[i] < true_value) & (true_value < wald_r_month3_EY[i]), 1, 0)
        Covg_wald_Bt[i] <- ifelse((wald_l_month3_Bt[i] < true_value) & (true_value < wald_r_month3_Bt[i]), 1, 0)
        
        # Load the JN variance
        # if(weight_type == "weight1"){
        #   temp_JN_control <- try(readRDS(paste0("simu_control/para_boot_weight1_400/sim_paraboot_crl_logit_200_", i)), silent = T)
        #   temp_JN_treat <- try(readRDS(paste0("simu_treat/para_boot_weight1_400/sim_paraboot_trt_logit_200_", i)), silent = T)
        # }else{
        #   temp_JN_control <- try(readRDS(paste0("simu_control/para_boot_weight2_400/sim_paraboot_crl_logit_200_", i)), silent = T)
        #   temp_JN_treat <- try(readRDS(paste0("simu_treat/para_boot_weight2_400/sim_paraboot_trt_logit_200_", i)), silent = T)
        # }
        # 
        if(weight_type == "weight1"){
          temp_JN_control <- try(readRDS(paste0("simu_control/para_boot_weight1/sim_paraboot_crl_logit_200_", i)), silent = T)
          temp_JN_treat <- try(readRDS(paste0("simu_treat/para_boot_weight1/sim_paraboot_trt_logit_200_", i)), silent = T)
        }else{
          temp_JN_control <- try(readRDS(paste0("simu_control/para_boot_weight2/sim_paraboot_crl_logit_200_", i)), silent = T)
          temp_JN_treat <- try(readRDS(paste0("simu_treat/para_boot_weight2/sim_paraboot_trt_logit_200_", i)), silent = T)
        }
        
        Est_JN_control <- rep(NA, boot_size)
        Est_JN_treat <- rep(NA, boot_size)
        Est_JN_effect_EY <- rep(NA, boot_size)
        Est_JN_effect_Bt <- rep(NA, boot_size)
        Var_JN_effect_EY <- rep(NA, boot_size)
        
        if(length(temp_JN_control) > 1 & length(temp_JN_treat) > 1){
        
            for(j in 1:boot_size){
              
              Est_JN_control[j] <- temp_JN_control[[j]][[1]][index_crl, index_month]
              Est_JN_treat[j] <- temp_JN_treat[[j]][[1]][index_trt, index_month]
              Est_JN_effect_EY[j] <- Est_JN_treat[j] - Est_JN_control[j]
              Est_JN_effect_Bt[j] <-  atanh(Est_JN_treat[j] - Est_JN_control[j])
              Var_JN_effect_EY[j] <- temp_JN_control[[j]][[2]][index_crl, index_month] + temp_JN_treat[[j]][[2]][index_trt, index_month]
             
            }
            
            # normal CI using bootstrap se
            temp4_EY <- wald(Effect[i], var(Est_JN_effect_EY))
            paraboot_l_month3_EY_se[i]  <- temp4_EY[1]
            paraboot_r_month3_EY_se[i]  <- temp4_EY[2]
            
            temp4_Bt <- wald(atanh(Effect[i]), var(Est_JN_effect_Bt))
            paraboot_l_month3_Bt_se[i]   <- tanh(temp4_Bt[1])
            paraboot_r_month3_Bt_se[i]   <- tanh(temp4_Bt[2])
            
            # Coverage
            Covg_paraboot_EY_se[i] <-  ifelse((paraboot_l_month3_EY_se[i] < true_value) & (true_value < paraboot_r_month3_EY_se[i]), 1, 0)
            Covg_paraboot_Bt_se[i] <-  ifelse((paraboot_l_month3_Bt_se[i] < true_value) & (true_value < paraboot_r_month3_Bt_se[i]), 1, 0)
        
          }
        
      }
      
      result <- data.frame(Emp_mean = mean(Effect), 
                           bias = mean(Effect) - true_value,
                           Wald_EY = mean(Covg_wald_EY, na.rm = T),
                           Wald_Bt = mean(Covg_wald_Bt, na.rm = T),
                           paraboot_se_EY = mean(Covg_paraboot_EY_se, na.rm = T),
                           paraboot_se_Bt = mean(Covg_paraboot_Bt_se, na.rm = T)
      )
      return(result)
    }
    
  }
  
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/logistic_simu_logit_norm1")
  # setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/logistic_simu_logit_fixcoef")
  # setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/single_simu_logit_norm1")
  
  
  # weight 1
  {
    
    # orig_treat <-  readRDS("simu_treat/ARC_simu_res_trt_dfsane_weight1_1e7_200_GL_30")
    # orig_control <-  readRDS("simu_control/ARC_simu_res_crl_dfsane_weight1_1e7_200_GL_30")
    orig_treat <-  readRDS("simu_treat/ARC_simu_res_trt_dfsane_weight1_1e7_200_GL_30_norm1")
    orig_control <-  readRDS("simu_control/ARC_simu_res_crl_dfsane_weight1_1e7_200_GL_30_norm1")
    
    
    # Get all results for one month for Question 1; risk difference
    {
      
      month <- 6
      Result <- c()
      
      for(i in 6:1){
        
        print(i)
        result <- c()
        alpha_trt <- Alpha_seq1[i]
        
        for(j in c(1:6)){ # this is for table in paper 
          
          # for question 1
          alpha_crl <- Alpha_seq2[j]
          
          # for question 1
          index_trt <- which(Alpha_seq1 == alpha_trt)
          index_crl <- which(Alpha_seq2 == alpha_crl)
          
          # for question 1
          true_value <- month6_treat[index_trt] - month6_control[index_crl]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight1",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_EY, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
        
      }
      
      # month 6
      Result |> abs()
      
      
      month <- 12
      Result <- c()
      
      for(i in 6:1){
        print(i)
        result <- c()
        alpha_trt <- Alpha_seq1[i]
        
        for(j in c(1:6)){ # this is for table in paper 
          
          # for question 1
          alpha_crl <- Alpha_seq2[j]
          
          # for question 1
          index_trt <- which(Alpha_seq1 == alpha_trt)
          index_crl <- which(Alpha_seq2 == alpha_crl)
          
          # for question 1
          true_value <- month12_treat[index_trt] - month12_control[index_crl]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight1",
                            true_value)
          
          # print(result1)
          result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_EY, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
      }
      
      # month 12
      Result |> abs()
      
      
    }
    
    
    # Get all results for one month for Question 2; risk difference
    {
      
      month <- 6
      Result <- c()
      
      for(i in 6:1){
        
        print(i)
        result <- c()
        
        for(j in c(1:6)){ # this is for table in paper 
          
          # for question 2
          alpha_trt <- 0
          alpha_crl <- 0
          
          # for question 2
          true_value <- month6_treat[i] - month6_control[j]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight1",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_EY, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
        
      }
      
      # month 6
      Result |> abs()
     
      
      
      month <- 12
      Result <- c()
      
      for(i in 6:1){
        
        print(i)
        result <- c()
        
        for(j in c(1:6)){ # this is for table in paper 
          
          # for question 2
          alpha_trt <- 0
          alpha_crl <- 0
          
          # for question 2
          true_value <- month12_treat[i] - month12_control[j]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight1",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_EY, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
        
      }
      
      # month 12
      Result |> abs()
      
      
      
    }
    
  }
  
  
  # weight 2
  {
    
    
    orig_treat <-  readRDS("simu_treat/ARC_simu_res_trt_dfsane_weight2_1e7_200_GL_30_norm1")
    orig_control <-  readRDS("simu_control/ARC_simu_res_crl_dfsane_weight2_1e7_200_GL_30_norm1")
    
    
    # Get all results for one month for Question 1; risk difference
    {
      
      month <- 6
      Result <- c()
      
      for(i in 6:1){
        print(i)
        result <- c()
        alpha_trt <- Alpha_seq1[i]
        
        for(j in 1:6){ # this is for table in paper 
          
          # for question 1
          alpha_crl <- Alpha_seq2[j]
          
          # for question 1
          index_trt <- which(Alpha_seq1 == alpha_trt)
          index_crl <- which(Alpha_seq2 == alpha_crl)
          
          # for question 1
          true_value <- month6_treat[index_trt] - month6_control[index_crl]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight2",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_Bt, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
        
      }
      
      # month 6
      Result |> abs()
      
      
      month <- 12
      Result <- c()
      
      for(i in 6:1){
        print(i)
        result <- c()
        alpha_trt <- Alpha_seq1[i]
        
        for(j in 1:6){ # this is for table in paper 
          
          # for question 1
          alpha_crl <- Alpha_seq2[j]
          
          # for question 1
          index_trt <- which(Alpha_seq1 == alpha_trt)
          index_crl <- which(Alpha_seq2 == alpha_crl)
          
          # for question 1
          true_value <- month12_treat[index_trt] - month12_control[index_crl]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight2",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_Bt, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
      }
      
      # month 12
      Result |> abs()
      
    }
    
    
    # Get all results for one month for Question 2; risk difference
    {
      
      month <- 6
      Result <- c()
      
      for(i in 6:1){
        
        print(i)
        result <- c()
        
        for(j in c(1:6)){ # this is for table in paper 
          
          # for question 2
          alpha_trt <- 0
          alpha_crl <- 0
          
          # for question 2
          true_value <- month6_treat[i] - month6_control[j]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight2",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_EY, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
        
      }
      
      # month 6
      Result |> abs()
      
      
      
      month <- 12
      Result <- c()
      
      for(i in 6:1){
        
        print(i)
        result <- c()
        
        for(j in c(1:6)){ # this is for table in paper 
          
          # for question 2
          alpha_trt <- 0
          alpha_crl <- 0
          
          # for question 2
          true_value <- month12_treat[i] - month12_control[j]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight2",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_EY, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
        
      }
      
      # month 12
      Result |> abs()
      
      
      
    }
    
  }
  
  
}


### old code ###
{
  
  ### setting
  {
    
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
    
    Alpha_seq1 <- c(-0.6, -0.3, 0, 0.3, 0.6)
    Alpha_seq2 <- c(-0.6, -0.3, 0, 0.3, 0.6)
    
    N_size <- 200
    
    wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
      CI_left <- mean_orig - Zstar *sqrt(var_orig)
      CI_right <- mean_orig + Zstar *sqrt(var_orig)
      
      return(c(CI_left, CI_right))
    }
    
    atanh <- function(x){
      return(1/2 * log((1 + x)/(1 - x)))
    }
    tanh <- function(x){
      return((exp(x) - exp(-x)) / (exp(x) + exp(-x)))
    }
    
    f_diff <- function(alpha_trt, alpha_crl, month,
                       weight_type = "weight1",
                       true_value){
      
      Alpha_seq1 <- c(-0.6, -0.3, 0, 0.3, 0.6)
      Alpha_seq2 <- c(-0.6, -0.3, 0, 0.3, 0.6)
      index_trt <- which(Alpha_seq1 == alpha_trt)
      index_crl <- which(Alpha_seq2 == alpha_crl)
      
      index_month <- which(c(1, 3, 6, 9, 12) == month)
      
      N <- 500
      
      # for original data's estimation
      Month_3_trt <- rep(NA, N)
      Var_3_trt <- rep(NA, N) 
      Month_3_crl <- rep(NA, N)
      Var_3_crl <- rep(NA, N) 
      
      Effect <-  rep(NA, N)
      Var_effect <-  rep(NA, N)
      C_month3 <- rep(NA, N)
      
      # for confidence interval
      wald_l_month3 <- rep(NA, N)
      wald_r_month3 <- rep(NA, N)
      wald_l_month3_trans <- rep(NA, N)
      wald_r_month3_trans <- rep(NA, N)
      
      JN_l_month3 <- rep(NA, N)
      JN_r_month3 <- rep(NA, N)
      JN_l_month3_trans <- rep(NA, N)
      JN_r_month3_trans <- rep(NA, N)
      
      # for coverage
      Covg_wald <- rep(NA, N)
      Covg_wald_trans <- rep(NA, N)
      
      Covg_JN <- rep(NA, N)
      Covg_JN_trans <- rep(NA, N)
      
      for(i in 1:N){
        
        # print(i)
        
        Month_3_trt[i] <- orig_treat[[i]][[1]][index_trt, index_month]
        Var_3_trt[i] <- orig_treat[[i]][[2]][index_trt, index_month]
        
        Month_3_crl[i] <- orig_control[[i]][[1]][index_crl, index_month]
        Var_3_crl[i] <- orig_control[[i]][[2]][index_crl, index_month]
        
        Effect[i] <- Month_3_trt[i] - Month_3_crl[i]
        Var_effect[i] <- Var_3_trt[i] + Var_3_crl[i]
        
        # Wald confidence interval - original
        temp2 <- wald(Effect[i], Var_effect[i])
        wald_l_month3[i]  <- temp2[1]
        wald_r_month3[i]  <- temp2[2]
        
        # Wald confidence interval - trans
        temp3 <- wald(atanh(Effect[i]), Var_effect[i] / ( 1 - Effect[i]^2 )^2 )
        wald_l_month3_trans[i]  <- tanh(temp3[1])
        wald_r_month3_trans[i]  <- tanh(temp3[2])
        
        # coverage
        Covg_wald[i] <- ifelse((wald_l_month3[i] < true_value) & (true_value < wald_r_month3[i]), 1, 0)
        Covg_wald_trans[i] <- ifelse((wald_l_month3_trans[i] < true_value) & (true_value < wald_r_month3_trans[i]), 1, 0)
        
        if(i %in% data_list_effect){
          
          # Load the JN variance
          if(weight_type == "weight1"){
            temp_JN_control <- readRDS(paste0("simu_control/JN_simu_control_weight1_1e7/simu_JN_", i, "_crl"))
            temp_JN_treat <- readRDS(paste0("simu_treat/JN_simu_treat_weight1_1e7/simu_JN_", i, "_trt"))
            
          }else{
            temp_JN_control <- readRDS(paste0("simu_control/JN_simu_control_weight2_1e7/simu_JN_", i, "_crl"))
            temp_JN_treat <- readRDS(paste0("simu_treat/JN_simu_treat_weight2_1e7/simu_JN_", i, "_trt"))
          }
          
          Est_JN_control <- rep(NA, N_size)
          Est_JN_treat <- rep(NA, N_size)
          Est_JN_effect <- rep(NA, N_size)
          Est_JN_effect_trans <- rep(NA, N_size)
          
          for(j in 1:N_size){
            temp_crl <- try(temp_JN_control[[j]][[1]], silent = T)
            temp_trt <- try(temp_JN_treat[[j]][[1]], silent = T)
            
            if(length(temp_crl) > 1 & length(temp_trt) > 1){
              
              Est_JN_control[j] <- temp_JN_control[[j]][[1]][index_crl, index_month]
              Est_JN_treat[j] <- temp_JN_treat[[j]][[1]][index_trt, index_month]
              Est_JN_effect[j] <- Est_JN_treat[j] - Est_JN_control[j]
              Est_JN_effect_trans[j] <-  atanh(Est_JN_treat[j] - Est_JN_control[j])
            }
          }
          
          if(sum(is.na(Est_JN_control)) > 0 & sum(is.na(Est_JN_treat)) > 0){
            
            Covg_JN[i] <- NA
            Covg_JN_trans[i] <- NA
            # Covg_JN_perc[i] <- NA
            # Covg_JN_seudo[i] <- NA
            
          }else{
            
            var_JN_control <- (N_size - 1)/N_size * sum((Est_JN_control - mean(Est_JN_control))^2)
            var_JN_treat <- (N_size - 1)/N_size * sum((Est_JN_treat - mean(Est_JN_treat))^2)
            # var_JN_effect <- (N_size - 1)/N_size *  sum((Est_JN_effect - mean(Est_JN_effect))^2)
            var_JN_effect <- var_JN_control + var_JN_treat
            var_JN_effect_trans <- (N_size - 1)/N_size *  sum((Est_JN_effect_trans - mean(Est_JN_effect_trans))^2)
            
            # JN confidence interval - original
            temp4 <- wald(Effect[i], var_JN_effect)
            JN_l_month3[i]  <- temp4[1]
            JN_r_month3[i]  <- temp4[2]
            
            # JN confidence interval - transformed
            temp5 <- wald(atanh(Effect[i]), var_JN_effect_trans)
            JN_l_month3_trans[i]  <- tanh(temp5[1])
            JN_r_month3_trans[i]  <- tanh(temp5[2])
            
            # Coverage
            Covg_JN[i] <-  ifelse((JN_l_month3[i] < true_value) & (true_value < JN_r_month3[i]), 1, 0)
            Covg_JN_trans[i] <-  ifelse((JN_l_month3_trans[i] < true_value) & (true_value < JN_r_month3_trans[i]), 1, 0)
            
          }
          
        }
        
      }
      
      result <- data.frame(Emp_mean = mean(Effect), 
                           bias = mean(Effect) - true_value,
                           Wald = mean(Covg_wald),
                           Wald_trans = mean(Covg_wald_trans),
                           JN = mean(Covg_JN, na.rm = T),
                           JN_trans = mean(Covg_JN_trans, na.rm = T)
      )
      return(result)
    }
    
    f_ratio <- function(alpha_trt, alpha_crl, month,
                        weight_type = "weight1",
                        true_value){
      
      Alpha_seq1 <- c(-0.6, -0.3, 0, 0.3, 0.6)
      Alpha_seq2 <- c(-0.6, -0.3, 0, 0.3, 0.6)
      index_trt <- which(Alpha_seq1 == alpha_trt)
      index_crl <- which(Alpha_seq2 == alpha_crl)
      
      index_month <- which(c(1, 3, 6, 9, 12) == month)
      
      N <- 500
      
      # for original data's estimation
      Month_3_trt <- rep(NA, N)
      Var_3_trt <- rep(NA, N) 
      Month_3_crl <- rep(NA, N)
      Var_3_crl <- rep(NA, N) 
      
      Effect <-  rep(NA, N)
      Var_effect <-  rep(NA, N)
      C_month3 <- rep(NA, N)
      
      # for confidence interval
      wald_l_month3 <- rep(NA, N)
      wald_r_month3 <- rep(NA, N)
      wald_l_month3_trans <- rep(NA, N)
      wald_r_month3_trans <- rep(NA, N)
      
      JN_l_month3 <- rep(NA, N)
      JN_r_month3 <- rep(NA, N)
      JN_l_month3_trans <- rep(NA, N)
      JN_r_month3_trans <- rep(NA, N)
      
      # for coverage
      Covg_wald <- rep(NA, N)
      Covg_wald_trans <- rep(NA, N)
      
      Covg_JN <- rep(NA, N)
      Covg_JN_trans <- rep(NA, N)
      
      for(i in 1:N){
        
        Month_3_trt[i] <- orig_treat[[i]][[1]][index_trt, index_month]
        Var_3_trt[i] <- orig_treat[[i]][[2]][index_trt, index_month]
        
        Month_3_crl[i] <- orig_control[[i]][[1]][index_crl, index_month]
        Var_3_crl[i] <- orig_control[[i]][[2]][index_crl, index_month]
        
        Effect[i] <- Month_3_trt[i] / Month_3_crl[i]
        Var_effect[i] <- Var_3_trt[i] / (Month_3_crl[i])^2 + (Month_3_trt[i])^2 / (Month_3_crl[i])^4 *Var_3_crl[i]
        
        # Wald confidence interval - original
        temp2 <- wald(Effect[i], Var_effect[i])
        wald_l_month3[i]  <- temp2[1]
        wald_r_month3[i]  <- temp2[2]
        
        # Wald confidence interval - trans
        temp3 <- wald(log(Effect[i]), Var_3_trt[i] / (Month_3_trt[i])^2 + Var_3_crl[i] / (Month_3_crl[i])^2 )
        wald_l_month3_trans[i]  <- exp(temp3[1])
        wald_r_month3_trans[i]  <- exp(temp3[2])
        
        # coverage
        Covg_wald[i] <- ifelse((wald_l_month3[i] < true_value) & (true_value < wald_r_month3[i]), 1, 0)
        Covg_wald_trans[i] <- ifelse((wald_l_month3_trans[i] < true_value) & (true_value < wald_r_month3_trans[i]), 1, 0)
        
        if(i %in% data_list_effect){
          
          # Load the JN variance
          if(weight_type == "weight1"){
            temp_JN_control <- readRDS(paste0("simu_control/JN_simu_control_weight1_1e7/simu_JN_", i, "_crl"))
            temp_JN_treat <- readRDS(paste0("simu_treat/JN_simu_treat_weight1_1e7/simu_JN_", i, "_trt"))
            
          }else{
            temp_JN_control <- readRDS(paste0("simu_control/JN_simu_control_weight2_1e7/simu_JN_", i, "_crl"))
            temp_JN_treat <- readRDS(paste0("simu_treat/JN_simu_treat_weight2_1e7/simu_JN_", i, "_trt"))
          }
          
          Est_JN_control <- rep(NA, N_size)
          Est_JN_treat <- rep(NA, N_size)
          Est_JN_effect <- rep(NA, N_size)
          Est_JN_effect_trans <- rep(NA, N_size)
          
          for(j in 1:N_size){
            temp_crl <- try(temp_JN_control[[j]][[1]], silent = T)
            temp_trt <- try(temp_JN_treat[[j]][[1]], silent = T)
            
            if(length(temp_crl) > 1 & length(temp_trt) > 1){
              
              Est_JN_control[j] <- temp_JN_control[[j]][[1]][index_crl, index_month]
              Est_JN_treat[j] <- temp_JN_treat[[j]][[1]][index_trt, index_month]
              Est_JN_effect[j] <- Est_JN_treat[j]/ Est_JN_control[j]
              Est_JN_effect_trans[j] <-  log(Est_JN_treat[j]/ Est_JN_control[j])
            }
          }
          
          if(sum(is.na(Est_JN_control)) > 0 & sum(is.na(Est_JN_treat)) > 0){
            
            Covg_JN[i] <- NA
            Covg_JN_trans[i] <- NA
            # Covg_JN_perc[i] <- NA
            # Covg_JN_seudo[i] <- NA
            
          }else{
            
            var_JN_effect <- (N_size - 1)/N_size * sum((Est_JN_effect - mean(Est_JN_effect))^2)
            var_JN_effect_trans <- (N_size - 1)/N_size *  sum((Est_JN_effect_trans - mean(Est_JN_effect_trans))^2)
            
            # JN confidence interval - original
            temp4 <- wald(Effect[i], var_JN_effect)
            JN_l_month3[i]  <- temp4[1]
            JN_r_month3[i]  <- temp4[2]
            
            # JN confidence interval - transformed
            temp5 <- wald(log(Effect[i]), var_JN_effect_trans)
            JN_l_month3_trans[i]  <- exp(temp5[1])
            JN_r_month3_trans[i]  <- exp(temp5[2])
            
            # Coverage
            Covg_JN[i] <-  ifelse((JN_l_month3[i] < true_value) & (true_value < JN_r_month3[i]), 1, 0)
            Covg_JN_trans[i] <-  ifelse((JN_l_month3_trans[i] < true_value) & (true_value < JN_r_month3_trans[i]), 1, 0)
            
          }
          
        }
        
      }
      
      result <- data.frame(Emp_mean = mean(Effect), 
                           bias = mean(Effect) - true_value,
                           Wald = mean(Covg_wald),
                           Wald_trans = mean(Covg_wald_trans),
                           JN = mean(Covg_JN, na.rm = T),
                           JN_trans = mean(Covg_JN_trans, na.rm = T)
      )
      return(result)
    }
    
  }
  
  
  ### logit link - weight1
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/binary_simu")
    orig_treat <-  readRDS("simu_treat/ARC_simu_res_trt_dfsane_weight1_1e7")
    orig_control <-  readRDS("simu_control/ARC_simu_res_crl_dfsane_weight1_1e7")
    
    # orig_treat <-  readRDS("simu_treat/ARC_simu_res_trt_dfsane_weight1_coef_1e7_300")
    # orig_control <-  readRDS("simu_control/ARC_simu_res_crl_dfsane_weight1_coef_1e7_300")
    
    data_list_JN_treat <- c() 
    for(k in 1:500){
      temp_data <- try( readRDS(paste0("simu_treat/JN_simu_treat_weight1_1e7/simu_JN_", k, "_trt")) , silent = TRUE)
      if(length(temp_data) > 1){
        data_list_JN_treat <- c(data_list_JN_treat, k)
      }
    }
    
    data_list_JN_control <- c() 
    for(k in 1:500){
      temp_data <- try( readRDS(paste0("simu_control/JN_simu_control_weight1_1e7/simu_JN_", k, "_crl")) , silent = TRUE)
      if(length(temp_data) > 1){
        data_list_JN_control <- c(data_list_JN_control, k)
      }
    }
    
    data_list_effect <- intersect(data_list_JN_treat, data_list_JN_control)
    length(data_list_effect)
    
    
    # Get all results for one month for Question 1; risk difference
    {
      
      month <- 6
      Result <- c()
      
      for(i in 5:1){
        result <- c()
        alpha_trt <- Alpha_seq1[i]
        
        for(j in 1:5){ # this is for table in paper 
          
          # for question 1
          alpha_crl <- Alpha_seq2[j]
          
          # for question 1
          index_trt <- which(Alpha_seq1 == alpha_trt)
          index_crl <- which(Alpha_seq2 == alpha_crl)
          
          # for question 1
          true_value <- month6_treat[index_trt] - month6_control[index_crl]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight1",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$Wald_trans, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
        
      }
      
      Result |> abs()
      
      # rownames(Result) <- c("alpha1_n06", "alpha1_n03", "alpha1_0",
      #                       "alpha1_03", "alpha1_06")
      # colnames(Result) <- c("alpha0_n06_bias", "alpha0_n06_cvg", 
      #                       "alpha0_n03_bias", "alpha0_n03_cvg", 
      #                       "alpha0_0_bias", "alpha0_0_cvg", 
      #                       "alpha0_03_bias", "alpha0_03_cvg", 
      #                       "alpha0_06_bias", "alpha0_06_cvg")
      # Result |> abs()
      
      
      month <- 12
      Result <- c()
      
      for(i in 5:1){
        result <- c()
        alpha_trt <- Alpha_seq1[i]
        
        for(j in 1:5){ # this is for table in paper 
          
          # for question 1
          alpha_crl <- Alpha_seq2[j]
          
          # for question 1
          index_trt <- which(Alpha_seq1 == alpha_trt)
          index_crl <- which(Alpha_seq2 == alpha_crl)
          
          # for question 1
          true_value <- month12_treat[index_trt] - month12_control[index_crl]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight1",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$Wald_trans, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
      }
      
      Result |> abs()
      
    }
    
    
    # Get all results for one month for Question 1; risk ratio 
    #### [the jackkknife version requires control and treatment have same sample size]
    {
      
      month <- 6
      Result <- c()
      
      for(i in 5:1){
        result <- c()
        alpha_trt <- Alpha_seq1[i]
        
        for(j in 1:5){ # this is for table in paper 
          
          # for question 1
          alpha_crl <- Alpha_seq2[j]
          
          # for question 1
          index_trt <- which(Alpha_seq1 == alpha_trt)
          index_crl <- which(Alpha_seq2 == alpha_crl)
          
          # for question 1
          true_value <- month6_treat[index_trt] / month6_control[index_crl]
          # true_value <- month12_treat[index_trt] / month12_control[index_crl]
          
          result1 <- f_ratio(alpha_trt, alpha_crl, month, 
                             weight_type = "weight1",
                             true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$Wald, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
        
      }
      
      Result |> abs()
      
      # rownames(Result) <- c("alpha1_n06", "alpha1_n03", "alpha1_0",
      #                       "alpha1_03", "alpha1_06")
      # colnames(Result) <- c("alpha0_n06_bias", "alpha0_n06_cvg", 
      #                       "alpha0_n03_bias", "alpha0_n03_cvg", 
      #                       "alpha0_0_bias", "alpha0_0_cvg", 
      #                       "alpha0_03_bias", "alpha0_03_cvg", 
      #                       "alpha0_06_bias", "alpha0_06_cvg")
      # Result |> abs()
      
      
      month <- 12
      Result <- c()
      
      for(i in 5:1){
        result <- c()
        alpha_trt <- Alpha_seq1[i]
        
        for(j in 1:5){ # this is for table in paper 
          
          # for question 1
          alpha_crl <- Alpha_seq2[j]
          
          # for question 1
          index_trt <- which(Alpha_seq1 == alpha_trt)
          index_crl <- which(Alpha_seq2 == alpha_crl)
          
          # for question 1
          # true_value <- month3_treat[index_trt] - month3_control[index_crl]
          # true_value <- month6_treat[index_trt] - month6_control[index_crl]
          # true_value <- month9_treat[index_trt] - month9_control[index_crl]
          true_value <- month12_treat[index_trt] - month12_control[index_crl]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            link_type = "log", weight_type = "weight1",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$JN, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
      }
      
      Result |> abs()
      
    }
    
    
    # Get all results for one month for Question 2; risk difference
    {
      
      month <- 6
      Result <- c()
      
      for(i in 5:1){
        result <- c()
        
        for(j in 1:5){
          # for question 2
          alpha_trt <- 0
          alpha_crl <- 0
          
          # for question 2
          true_value <- month6_treat[i] - month6_control[j]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight1",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$Wald_trans, 3))
        }
        Result <- rbind(Result, result)
      }
      
      Result |> abs()
      
      
      
      month <- 12
      Result <- c()
      
      for(i in 5:1){
        result <- c()
        
        for(j in 1:5){
          # for question 2
          alpha_trt <- 0
          alpha_crl <- 0
          
          # for question 2
          true_value <- month12_treat[i] - month12_control[j]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight1",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$Wald_trans, 3))
        }
        Result <- rbind(Result, result)
      }
      
      Result |> abs()
      
    }
    
    
  }
  
  
  ### log link - weight2
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/binary_simu")
    orig_treat <-  readRDS("simu_treat/ARC_simu_res_trt_dfsane_weight2_1e7")
    orig_control <-  readRDS("simu_control/ARC_simu_res_crl_dfsane_weight2_1e7")
    
    data_list_JN_treat <- c() 
    for(k in 1:500){
      temp_data <- try( readRDS(paste0("simu_treat/JN_simu_treat_weight2_1e7/simu_JN_", k, "_trt")) , silent = TRUE)
      if(length(temp_data) > 1){
        data_list_JN_treat <- c(data_list_JN_treat, k)
      }
    }
    
    data_list_JN_control <- c() 
    for(k in 1:500){
      temp_data <- try( readRDS(paste0("simu_control/JN_simu_control_weight2_1e7/simu_JN_", k, "_crl")) , silent = TRUE)
      if(length(temp_data) > 1){
        data_list_JN_control <- c(data_list_JN_control, k)
      }
    }
    
    data_list_effect <- intersect(data_list_JN_treat, data_list_JN_control)
    length(data_list_effect)
    
    
    # Get all results for one month for Question 1
    {
      
      month <- 6
      Result <- c()
      
      for(i in 5:1){
        result <- c() # for question 1
        alpha_trt <- Alpha_seq1[i]
        
        for(j in 1:5){ # this is for table in paper 
          
          # for question 1
          alpha_crl <- Alpha_seq2[j]
          
          # for question 1
          index_trt <- which(Alpha_seq1 == alpha_trt)
          index_crl <- which(Alpha_seq2 == alpha_crl)
          
          # for question 1
          # true_value <- month3_treat[index_trt] - month3_control[index_crl]
          true_value <- month6_treat[index_trt] - month6_control[index_crl]
          # true_value <- month9_treat[index_trt] - month9_control[index_crl]
          # true_value <- month12_treat[index_trt] - month12_control[index_crl]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            weight_type = "weight2",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$JN, 3))
          
        }
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
      }
      
      Result |> abs()
      
      # rownames(Result) <- c("alpha1_n06", "alpha1_n03", "alpha1_0",
      #                       "alpha1_03", "alpha1_06")
      # colnames(Result) <- c("alpha0_n06_bias", "alpha0_n06_cvg", 
      #                       "alpha0_n03_bias", "alpha0_n03_cvg", 
      #                       "alpha0_0_bias", "alpha0_0_cvg", 
      #                       "alpha0_03_bias", "alpha0_03_cvg", 
      #                       "alpha0_06_bias", "alpha0_06_cvg")
      # Result |> abs()
      
      month <- 12
      Result <- c()
      
      for(i in 5:1){
        result <- c()
        # for question 1
        
        alpha_trt <- Alpha_seq1[i]
        
        for(j in 1:5){ # this is for table in paper 
          
          # for question 1
          alpha_crl <- Alpha_seq2[j]
          
          # for question 1
          index_trt <- which(Alpha_seq1 == alpha_trt)
          index_crl <- which(Alpha_seq2 == alpha_crl)
          
          # for question 1
          # true_value <- month3_treat[index_trt] - month3_control[index_crl]
          # true_value <- month6_treat[index_trt] - month6_control[index_crl]
          # true_value <- month9_treat[index_trt] - month9_control[index_crl]
          true_value <- month12_treat[index_trt] - month12_control[index_crl]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            link_type = "log", weight_type = "weight2",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$JN, 3))
          
        }
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
      }
      
      Result |> abs()
      
    }
    
    
    # Get all results for one month for Question 2
    {
      
      month <- 6
      Result <- c()
      
      for(i in 5:1){
        
        result <- c()
        # # for question 1
        # alpha_crl <- Alpha_seq2[2]
        
        for(j in 1:5){
          
          # for question 2
          alpha_trt <- 0
          alpha_crl <- 0
          
          # for question 2
          # true_value <- month3_treat[i] - month3_control[j]
          true_value <- month6_treat[i] - month6_control[j]
          # true_value <- month9_treat[i] - month9_control[j]
          # true_value <- month12_treat[i] - month12_control[j]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month,
                            link_type = "log", weight_type = "weight2",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$JN, 3))
        }
        Result <- rbind(Result, result)
      }
      
      Result |> abs()
      
      # rownames(Result) <- c("alpha1_n06", "alpha1_n03", "alpha1_0",
      #                       "alpha1_03", "alpha1_06")
      # colnames(Result) <- c("alpha0_n06_bias", "alpha0_n06_cvg", 
      #                       "alpha0_n03_bias", "alpha0_n03_cvg", 
      #                       "alpha0_0_bias", "alpha0_0_cvg", 
      #                       "alpha0_03_bias", "alpha0_03_cvg", 
      #                       "alpha0_06_bias", "alpha0_06_cvg")
      # Result |> abs()
      
      month <- 12
      Result <- c()
      
      for(i in 5:1){
        
        result <- c()
        # # for question 1
        # alpha_crl <- Alpha_seq2[2]
        
        for(j in 1:5){
          
          # for question 2
          alpha_trt <- 0
          alpha_crl <- 0
          
          # for question 2
          # true_value <- month3_treat[i] - month3_control[j]
          # true_value <- month6_treat[i] - month6_control[j]
          # true_value <- month9_treat[i] - month9_control[j]
          true_value <- month12_treat[i] - month12_control[j]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month,
                            link_type = "log", weight_type = "weight2",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$JN, 3))
        }
        Result <- rbind(Result, result)
      }
      
      Result |> abs()
      
    }
    
  }
  
}

