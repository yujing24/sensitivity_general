rm(list=ls())
library(dplyr)


### parametric bootstrap - single index model
{
  
  # setting
  {
    
    # month <- 3
    month3_treat <- c(1.601, 1.876, 2.191, 2.563, 3.008)
    month3_control <- c(1.551, 1.930, 2.342, 2.774, 3.211)
    
    # month <- 6
    month6_treat <- c(1.554, 1.821, 2.130, 2.496, 2.945)
    month6_control <- c(1.521, 1.895, 2.304, 2.736, 3.181)
    
    # month <- 9
    month9_treat <- c(1.514, 1.777, 2.076, 2.441, 2.891)
    month9_control <- c(1.495, 1.865, 2.269, 2.702, 3.150)
    
    # month <- 12
    month12_treat <- c(1.486, 1.743, 2.039, 2.402, 2.848)
    month12_control <- c(1.478, 1.842, 2.245, 2.679, 3.127)
    
    Alpha_seq1 <- c(-0.6, -0.3, 0, 0.3, 0.6)
    Alpha_seq2 <- c(-0.6, -0.3, 0, 0.3, 0.6)
    
    N <- 500
    boot_size <- 200
    
    wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
      CI_left <- mean_orig - Zstar *sqrt(var_orig)
      CI_right <- mean_orig + Zstar *sqrt(var_orig)
      
      return(c(CI_left, CI_right))
    }
    
    f_diff <- function(alpha_trt, alpha_crl, month,
                       weight_type = "weight1",
                       true_value){
      
      Alpha_seq1 <- c(-0.6, -0.3, 0, 0.3, 0.6)
      Alpha_seq2 <- c(-0.6, -0.3, 0, 0.3, 0.6)
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
      
      paraboot_l_month3_EY_se <- rep(NA, N)
      paraboot_r_month3_EY_se <- rep(NA, N)
      
      # for coverage
      Covg_wald_EY <- rep(NA, N)
      Covg_paraboot_EY_se <- rep(NA, N)
      
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
        
        # coverage
        Covg_wald_EY[i] <- ifelse((wald_l_month3_EY[i] < true_value) & (true_value < wald_r_month3_EY[i]), 1, 0)
        
        # Load the parametric bootstrap variance
        temp_paraboot_control <- readRDS(paste0("simu_control/para_boot/sim_paraboot_crl_identity_200_", i))
        temp_paraboot_treat <- readRDS(paste0("simu_treat/para_boot/sim_paraboot_trt_identity_200_", i))
        
        Est_paraboot_control <- rep(NA, boot_size)
        Est_paraboot_treat <- rep(NA, boot_size)
        Est_paraboot_effect_EY <- rep(NA, boot_size)
        Var_paraboot_effect_EY <- rep(NA, boot_size)
        
        for(j in 1:boot_size){
          
          Est_paraboot_control[j] <- temp_paraboot_control[[j]][[1]][index_crl, index_month]
          Est_paraboot_treat[j] <- temp_paraboot_treat[[j]][[1]][index_trt, index_month]
          Est_paraboot_effect_EY[j] <- Est_paraboot_treat[j] - Est_paraboot_control[j]
          Var_paraboot_effect_EY[j] <- temp_paraboot_control[[j]][[2]][index_crl, index_month] + temp_paraboot_treat[[j]][[2]][index_trt, index_month]
          
        }
        
        # normal CI using bootstrap se
        temp4_EY <- wald(Effect[i], var(Est_paraboot_effect_EY))
        paraboot_l_month3_EY_se[i]  <- temp4_EY[1]
        paraboot_r_month3_EY_se[i]  <- temp4_EY[2]
  
        # Coverage
        Covg_paraboot_EY_se[i] <-  ifelse((paraboot_l_month3_EY_se[i] < true_value) & (true_value < paraboot_r_month3_EY_se[i]), 1, 0)
        
      }
      
      result <- data.frame(Emp_mean = mean(Effect), 
                           bias = mean(Effect) - true_value,
                           Wald_EY = mean(Covg_wald_EY),
                           paraboot_se_EY = mean(Covg_paraboot_EY_se)
      )
      return(result)
    }
    
  }
  
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/single_simu_identity")
  
  
  # identity
  {
    
    
    orig_treat <-  readRDS("simu_treat/ARC_simu_res_trt_dfsane_1e7_200_GL_30")
    orig_control <-  readRDS("simu_control/ARC_simu_res_crl_dfsane_1e7_200_GL_30")
    
    
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
          
          result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_EY, 3))
          
        }
        
        result <- c(alpha_trt, result)
        Result <- rbind(Result, result)
        
      }
      
      # month 6
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
          
          result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_EY, 3))
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
          
          result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_EY, 3))
        }
        Result <- rbind(Result, result)
      }
      
      Result |> abs()
      
    }
    
  }
  
  
  
}



#### old code ####
{
  
  ### setting
  {
    
    # month <- 3
    month3_treat <- c(1.529, 1.781, 2.094, 2.472, 2.907)
    month3_control <- c(1.441, 1.754, 2.165, 2.675, 3.249)
    
    # month <- 6
    month6_treat <- c(1.479, 1.718, 2.014, 2.372, 2.787)
    month6_control <- c(1.394, 1.691, 2.081, 2.567, 3.125)
    
    # month <- 9
    month9_treat <- c(1.450, 1.681, 1.968, 2.315, 2.721)
    month9_control <- c(1.360, 1.646, 2.021, 2.493, 3.041)
    
    # month <- 12
    month12_treat <- c(1.426, 1.651, 1.930, 2.269, 2.667)
    month12_control <- c(1.332, 1.608, 1.971, 2.431, 2.971)
    
    Alpha_seq1 <- c(-0.6, -0.3, 0, 0.3, 0.6)
    Alpha_seq2 <- c(-0.6, -0.3, 0, 0.3, 0.6)
    
    N_size <- 200
    
    wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
      CI_left <- mean_orig - Zstar *sqrt(var_orig)
      CI_right <- mean_orig + Zstar *sqrt(var_orig)
      
      return(c(CI_left, CI_right))
    }
    
    f_diff <- function(alpha_trt, alpha_crl, month,
                       link_type = "identity", 
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
      
      # for coverage
      Covg_wald <- rep(NA, N)
      Covg_JN <- rep(NA, N)
      
      for(i in 1:N){
        
        Month_3_trt[i] <- orig_treat[[i]][[1]][index_trt, index_month]
        Var_3_trt[i] <- orig_treat[[i]][[2]][index_trt, index_month]
        
        Month_3_crl[i] <- orig_control[[i]][[1]][index_crl, index_month]
        Var_3_crl[i] <- orig_control[[i]][[2]][index_crl, index_month]
        
        Effect[i] <- Month_3_trt[i] - Month_3_crl[i]
        Var_effect[i] <- Var_3_trt[i] + Var_3_crl[i]
        
        # Wald confidence interval
        temp2 <- wald(Effect[i], Var_effect[i])
        wald_l_month3[i]  <- temp2[1]
        wald_r_month3[i]  <- temp2[2]
        # coverage
        Covg_wald[i] <- ifelse((wald_l_month3[i] < true_value) & (true_value < wald_r_month3[i]), 1, 0)
        
        if(i %in% data_list_effect){
          
          # Load the JN variance
          if(link_type == "identity"){
            temp_JN_control <- readRDS(paste0("simu_control/JN_simu_control_1e7/simu_JN_", i, "_crl"))
            temp_JN_treat <- readRDS(paste0("simu_treat/JN_simu_treat_1e7/simu_JN_", i, "_trt"))
          }else{
            if(weight_type == "weight1"){
              temp_JN_control <- readRDS(paste0("simu_control/JN_simu_control_weight1_1e7/simu_JN_", i, "_crl"))
              temp_JN_treat <- readRDS(paste0("simu_treat/JN_simu_treat_weight1_1e7/simu_JN_", i, "_trt"))
              
            }else{
              temp_JN_control <- readRDS(paste0("simu_control/JN_simu_control_weight2_1e7/simu_JN_", i, "_crl"))
              temp_JN_treat <- readRDS(paste0("simu_treat/JN_simu_treat_weight2_1e7/simu_JN_", i, "_trt"))
              
            }
          }
          
          Est_JN_control <- rep(NA, N_size)
          Est_JN_treat <- rep(NA, N_size)
          for(j in 1:N_size){
            Est_JN_control[j] <- temp_JN_control[[j]][[1]][index_crl, index_month]
            Est_JN_treat[j] <- temp_JN_treat[[j]][[1]][index_trt, index_month]
          }
          var_JN_control <- (N_size - 1)/N_size * sum((Est_JN_control - mean(Est_JN_control))^2)
          var_JN_treat <- (N_size - 1)/N_size * sum((Est_JN_treat - mean(Est_JN_treat))^2)
          var_JN_effect <- var_JN_control + var_JN_treat
          
          # JN confidence interval
          temp4 <- wald(Effect[i], var_JN_effect)
          JN_l_month3[i]  <- temp4[1]
          JN_r_month3[i]  <- temp4[2]
          
          # Coverage
          Covg_JN[i] <-  ifelse((JN_l_month3[i] < true_value) & (true_value < JN_r_month3[i]), 1, 0)
          
        }
        
      }
      
      result <- data.frame(Emp_mean = mean(Effect), 
                           bias = mean(Effect) - true_value,
                           Wald = mean(Covg_wald),
                           JN = mean(Covg_JN, na.rm = T)
      )
      return(result)
    }
    
  }
  
  
  ### identity link
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/identity_simu")
    orig_treat <-  readRDS("simu_treat/ARC_simu_res_trt_1e7")
    orig_control <-  readRDS("simu_control/ARC_simu_res_crl_1e7")
    
    data_list_JN_treat <- c() 
    for(k in 1:500){
      temp_data <- try( readRDS(paste0("simu_treat/JN_simu_treat_1e7/simu_JN_", k, "_trt")) , silent = TRUE)
      if(length(temp_data) > 1){
        data_list_JN_treat <- c(data_list_JN_treat, k)
      }
    }
    
    data_list_JN_control <- c() 
    for(k in 1:500){
      temp_data <- try( readRDS(paste0("simu_control/JN_simu_control_1e7/simu_JN_", k, "_crl")) , silent = TRUE)
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
        # for(j in 1:5){ # for data analysis
        result <- c()
        # for question 1
        
        for(j in 1:5){ # this is for table in paper 
          # for(i in 1:5){# for data analysis
          
          # for question 1
          alpha_trt <- Alpha_seq1[i]
          alpha_crl <- Alpha_seq2[j]
          
          # for question 1
          index_trt <- which(Alpha_seq1 == alpha_trt)
          index_crl <- which(Alpha_seq2 == alpha_crl)
          
          # true_value <- month3_treat[index_trt] - month3_control[index_crl]
          true_value <- month6_treat[index_trt] - month6_control[index_crl]
          # true_value <- month9_treat[index_trt] - month9_control[index_crl]
          # true_value <- month12_treat[index_trt] - month12_control[index_crl]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            link_type = "identity", weight_type = "weight1",
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
      
    }
    
    
    # Get all results for one month for Question 2
    {
      
      month <- 6
      Result <- c()
      
      for(i in 5:1){
        # for(j in 1:5){ # for data analysis 
        
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
                            link_type = "identity", weight_type = "weight1",
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
      
    }
    
    
  }
  
  
  ### log link - weight1
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/count_simu")
    orig_treat <-  readRDS("simu_treat/ARC_simu_res_trt_dfsane_weight1_1e7")
    orig_control <-  readRDS("simu_control/ARC_simu_res_crl_dfsane_weight1_1e7")
    
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
    
    
    # Get all results for one month for Question 1
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
          # true_value <- month3_treat[index_trt] - month3_control[index_crl]
          true_value <- month6_treat[index_trt] - month6_control[index_crl]
          # true_value <- month9_treat[index_trt] - month9_control[index_crl]
          # true_value <- month12_treat[index_trt] - month12_control[index_crl]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month, 
                            link_type = "log", weight_type = "weight1",
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
    
    
    # Get all results for one month for Question 2
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
          # true_value <- month3_treat[i] - month3_control[j]
          true_value <- month6_treat[i] - month6_control[j]
          # true_value <- month9_treat[i] - month9_control[j]
          # true_value <- month12_treat[i] - month12_control[j]
          
          result1 <- f_diff(alpha_trt, alpha_crl, month,
                            link_type = "log", weight_type = "weight1",
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
                            link_type = "log", weight_type = "weight1",
                            true_value)
          
          result <- c(result, round(result1$bias, 3), round(result1$JN, 3))
        }
        Result <- rbind(Result, result)
      }
      
      Result |> abs()
      
    }
    
    
  }
  
  
  ### log link - weight2
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/count_simu")
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
                            link_type = "log", weight_type = "weight2",
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
  
  
  ### parametric bootstrap - negative binomial model
  {
    
    # setting
    {
      
      # month <- 3
      month3_treat <- c(1.529, 1.781, 2.094, 2.472, 2.907)
      month3_control <- c(1.441, 1.754, 2.165, 2.675, 3.249)
      
      # month <- 6
      month6_treat <- c(1.479, 1.718, 2.014, 2.372, 2.787)
      month6_control <- c(1.394, 1.691, 2.081, 2.567, 3.125)
      
      # month <- 9
      month9_treat <- c(1.450, 1.681, 1.968, 2.315, 2.721)
      month9_control <- c(1.360, 1.646, 2.021, 2.493, 3.041)
      
      # month <- 12
      month12_treat <- c(1.426, 1.651, 1.930, 2.269, 2.667)
      month12_control <- c(1.332, 1.608, 1.971, 2.431, 2.971)
      
      Alpha_seq1 <- c(-0.6, -0.3, 0, 0.3, 0.6)
      Alpha_seq2 <- c(-0.6, -0.3, 0, 0.3, 0.6)
      
      N <- 500
      boot_size <- 200
      
      wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
        CI_left <- mean_orig - Zstar *sqrt(var_orig)
        CI_right <- mean_orig + Zstar *sqrt(var_orig)
        
        return(c(CI_left, CI_right))
      }
      
      f_diff <- function(alpha_trt, alpha_crl, month,
                         weight_type = "weight1",
                         true_value){
        
        Alpha_seq1 <- c(-0.6, -0.3, 0, 0.3, 0.6)
        Alpha_seq2 <- c(-0.6, -0.3, 0, 0.3, 0.6)
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
        
        # paraboot_l_month3_EY <- rep(NA, N)
        # paraboot_r_month3_EY <- rep(NA, N)
        paraboot_l_month3_EY_se <- rep(NA, N)
        paraboot_r_month3_EY_se <- rep(NA, N)
        paraboot_l_month3_Bt_se <- rep(NA, N)
        paraboot_r_month3_Bt_se <- rep(NA, N)
        
        # for coverage
        Covg_wald_EY <- rep(NA, N)
        Covg_wald_Bt <- rep(NA, N)
        # Covg_paraboot_EY <- rep(NA, N)
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
          temp3 <- wald(log(Month_3_trt[i] / Month_3_crl[i]), Var_3_trt[i] / (Month_3_trt[i]^2) + Var_3_crl[i] / (Month_3_crl[i]^2) )
          wald_l_month3_Bt[i]  <- exp(temp3[1])
          wald_r_month3_Bt[i]  <- exp(temp3[2])
          
          # coverage
          Covg_wald_EY[i] <- ifelse((wald_l_month3_EY[i] < true_value) & (true_value < wald_r_month3_EY[i]), 1, 0)
          Covg_wald_Bt[i] <- ifelse((wald_l_month3_Bt[i] < true_value) & (true_value < wald_r_month3_Bt[i]), 1, 0)
          
          # Load the parametric bootstrap variance
          temp_paraboot_control <- readRDS(paste0("simu_control/para_boot/sim_paraboot_crl_identity_200_", i))
          temp_paraboot_treat <- readRDS(paste0("simu_treat/para_boot/sim_paraboot_trt_identity_200_", i))
          
          Est_paraboot_control <- rep(NA, boot_size)
          Est_paraboot_treat <- rep(NA, boot_size)
          Est_paraboot_effect_EY <- rep(NA, boot_size)
          Est_paraboot_effect_Bt <- rep(NA, boot_size)
          Var_paraboot_effect_EY <- rep(NA, boot_size)
          
          for(j in 1:boot_size){
            
            # temp_crl <- try(temp_paraboot_control[[j]][[1]], silent = T)
            # temp_trt <- try(temp_paraboot_treat[[j]][[1]], silent = T)
            
            # if(length(temp_crl) > 1 & length(temp_trt) > 1){
            
            Est_paraboot_control[j] <- temp_paraboot_control[[j]][[1]][index_crl, index_month]
            Est_paraboot_treat[j] <- temp_paraboot_treat[[j]][[1]][index_trt, index_month]
            Est_paraboot_effect_EY[j] <- Est_paraboot_treat[j] - Est_paraboot_control[j]
            Est_paraboot_effect_Bt[j] <-  log( Est_paraboot_treat[j] / Est_paraboot_control[j] )
            Var_paraboot_effect_EY[j] <- temp_paraboot_control[[j]][[2]][index_crl, index_month] + temp_paraboot_treat[[j]][[2]][index_trt, index_month]
            # }
          }
          
          # if(sum(is.na(Est_JN_control)) > 0 & sum(is.na(Est_JN_treat)) > 0){
          #   
          #   Covg_JN_EY[i] <- NA
          #   Covg_JN_Bt[i] <- NA
          #   # Covg_JN_perc[i] <- NA
          #   # Covg_JN_seudo[i] <- NA
          #   
          # }else{
          #   
          # symmetric boot-t
          # boot_t_EY <- (Effect[i] -  Est_JN_effect_EY) / sqrt(Var_JN_effect_EY)
          # quant_EY <- quantile(abs(boot_t_EY), probs = c(0.95))
          # paraboot_l_month3_EY[i] <- Effect[i] - quant_EY[1] * sqrt(Var_effect[i])
          # paraboot_r_month3_EY[i] <- Effect[i] + quant_EY[1] * sqrt(Var_effect[i])
          
          # normal CI using bootstrap se
          temp4_EY <- wald(Effect[i], var(Est_paraboot_effect_EY))
          paraboot_l_month3_EY_se[i]  <- temp4_EY[1]
          paraboot_r_month3_EY_se[i]  <- temp4_EY[2]
          
          temp4_Bt <- wald(log(Month_3_trt[i] / Month_3_crl[i]), var(Est_paraboot_effect_Bt))
          paraboot_l_month3_Bt_se[i] <- exp(temp4_Bt[1])
          paraboot_r_month3_Bt_se[i] <- exp(temp4_Bt[2])
          
          # Coverage
          # Covg_paraboot_EY[i] <-  ifelse((paraboot_l_month3_EY[i] < true_value) & (true_value < paraboot_r_month3_EY[i]), 1, 0)
          Covg_paraboot_EY_se[i] <-  ifelse((paraboot_l_month3_EY_se[i] < true_value) & (true_value < paraboot_r_month3_EY_se[i]), 1, 0)
          Covg_paraboot_Bt_se[i] <-  ifelse((paraboot_l_month3_Bt_se[i] < true_value) & (true_value < paraboot_r_month3_Bt_se[i]), 1, 0)
          
          # }
          
        }
        
        result <- data.frame(Emp_mean = mean(Effect), 
                             bias = mean(Effect) - true_value,
                             Wald_EY = mean(Covg_wald_EY),
                             Wald_Bt = mean(Covg_wald_Bt),
                             # paraboot_t_EY = mean(Covg_paraboot_EY, na.rm = T),
                             paraboot_se_EY = mean(Covg_paraboot_EY_se),
                             paraboot_se_Bt = mean(Covg_paraboot_Bt_se)
        )
        return(result)
      }
      
    }
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/neg_simu_identity")
    
    
    # weight 1
    {
      
      
      orig_treat <-  readRDS("simu_treat/ARC_simu_res_trt_dfsane_1e7_200_GL_30")
      orig_control <-  readRDS("simu_control/ARC_simu_res_crl_dfsane_1e7_200_GL_30")
      
      
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
            
            result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_EY, 3))
            
          }
          
          result <- c(alpha_trt, result)
          Result <- rbind(Result, result)
          
        }
        
        # month 6
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
            
            result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_EY, 3))
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
            
            result <- c(result, round(result1$bias, 3), round(result1$paraboot_se_EY, 3))
          }
          Result <- rbind(Result, result)
        }
        
        Result |> abs()
        
      }
      
    }
    
    
    
  }
  
}




