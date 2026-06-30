rm(list=ls())
library(dplyr)


#### single index model + parametric bootstrap
{
  
  
  ### setting
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
    
    
    Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
    Month <- c(3, 6, 9, 12)
    N <- 500
    
    wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
      CI_left <- mean_orig - Zstar * sqrt(var_orig)
      CI_right <- mean_orig + Zstar * sqrt(var_orig)
      
      return(c(CI_left, CI_right))
    }
    
    f_point <- function(alpha, month, true_value){
      
      
      Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
      index_alpha <- which(Alpha_seq == alpha)
      index_month <- which(c(1, 3, 6, 9, 12) == month)
      
      # for original data's estimation
      Month_3 <- rep(NA, N)
      Var_3 <-rep(NA, N) 
      
      # confidence interval
      wald_l_month3_EY <- rep(NA, N)
      wald_r_month3_EY <- rep(NA, N)
      
      wald_l_month3_boot_EY <- rep(NA, N)
      wald_r_month3_boot_EY <- rep(NA, N)
      
      wald_l_month3_boot_EY_se <- rep(NA, N)
      wald_r_month3_boot_EY_se <- rep(NA, N)
      
      # for coverage
      Covg_wald_EY <- rep(NA, N)
      Covg_wald_boot_EY <- rep(NA, N)
      Covg_wald_boot_EY_se <- rep(NA, N)
      
      for(i in 1:N){
        
        Month_3[i] <- orig_treat[[i]][[1]][index_alpha, index_month]
        Var_3[i] <- orig_treat[[i]][[2]][index_alpha, index_month]
        
        # load the parametric bootstrap estimation results
        temp0 <- try(temp_paraboot <- readRDS(paste0("para_boot/sim_paraboot_trt_identity_200_", i)) , silent = T)
        
        if(length(temp0) > 1){
          
          ## E[Y(t)] scale
          # Wald confidence interval - on E[Y(t)] scale 
          temp1 <- wald(Month_3[i], Var_3[i])
          wald_l_month3_EY[i]  <- temp1[1]
          wald_r_month3_EY[i]  <- temp1[2]
          
          Est_paraboot_EY <- rep(NA, boot_size)
          Var_paraboot_EY <- rep(NA, boot_size)
          for(j in 1:boot_size){
            Est_paraboot_EY[j] <- temp_paraboot[[j]][[1]][index_alpha, index_month]
            Var_paraboot_EY[j] <- temp_paraboot[[j]][[2]][index_alpha, index_month]
          }
          
          ## boot-t 
          boot_t_EY <- (Month_3[i] -  Est_paraboot_EY) / sqrt(Var_paraboot_EY)
          quant_EY <- quantile(abs(boot_t_EY), probs = c(0.95))
          wald_l_month3_boot_EY[i] <- Month_3[i] - quant_EY[1] * sqrt(Var_3[i])
          wald_r_month3_boot_EY[i] <- Month_3[i] + quant_EY[1] * sqrt(Var_3[i])
          
          # normal CI using bootstrap se
          temp4_EY <- wald(Month_3[i], var(Est_paraboot_EY))
          wald_l_month3_boot_EY_se[i]  <- temp4_EY[1]
          wald_r_month3_boot_EY_se[i]  <- temp4_EY[2]
          
          
        }
        
        # Coverage
        Covg_wald_EY[i] <- ifelse((wald_l_month3_EY[i] < true_value) & (true_value < wald_r_month3_EY[i]), 1, 0)
        Covg_wald_boot_EY[i] <- ifelse((wald_l_month3_boot_EY[i] < true_value) & (true_value < wald_r_month3_boot_EY[i]), 1, 0)
        Covg_wald_boot_EY_se[i] <- ifelse((wald_l_month3_boot_EY_se[i] < true_value) & (true_value < wald_r_month3_boot_EY_se[i]), 1, 0)
        
      }
      
      #print(length(which(is.na(Covg_wald_boot_EY) == F)))
      
      result <- data.frame(Emp_mean = mean(Month_3), 
                           bias = mean(Month_3) - true_value,
                           Wald_EY = mean(Covg_wald_EY, na.rm = T), 
                           Wald_Bt_paraboot_EY = mean(Covg_wald_boot_EY, na.rm = T),
                           Wald_Bt_paraboot_EY_se = mean(Covg_wald_boot_EY_se, na.rm = T)
      )
      
      result
      
      return(result)
    }
    
    f_res_point <- function(alpha = -0.6){
      
      result <- c()
      for(month_i in 1:4){
        
        month <- Month[month_i]
        
        temp <- get(paste0("month", month, "_treat"))
        index_alpha <- which(Alpha_seq == alpha)
        true_value <- temp[index_alpha]
        
        result1 <- f_point(alpha, month, true_value)
        result <- rbind(result, c(alpha, true_value, round(result1, 3)))
        
      }
      colnames(result)[1:2] <- c("alpha", "true_value")
      result
      
    }
    
  }
  
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/single_simu_identity/simu_treat")
  sample_size <- 200
  boot_size <- 200
  N <- 500
  
  
  # identity link
  orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_1e7_200_GL_30")
  f_res_point(alpha = -0.6)
  f_res_point(alpha = -0.3)
  f_res_point(alpha = 0)
  f_res_point(alpha = 0.3)
  f_res_point(alpha = 0.6)
  
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/single_simu_identity/simu_treat")
  miss_index <- c()
  for(i in 1:N){
    temp0 <- try(temp_paraboot <- readRDS(paste0("para_boot/sim_paraboot_trt_identity_200_", i)) , silent = T)
    if(length(temp0) <= 1){
      miss_index <- c(miss_index, i)
    }
  }
  miss_index
  
  
}



#### negative binomial model + parametric bootstrap
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
    
    
    Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
    Month <- c(3, 6, 9, 12)
    N <- 500
    
    wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
      CI_left <- mean_orig - Zstar * sqrt(var_orig)
      CI_right <- mean_orig + Zstar * sqrt(var_orig)
      
      return(c(CI_left, CI_right))
    }
    
    f_point <- function(alpha, month, true_value){
      
      
      Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
      index_alpha <- which(Alpha_seq == alpha)
      index_month <- which(c(1, 3, 6, 9, 12) == month)
      
      # for original data's estimation
      Month_3 <- rep(NA, N)
      Var_3 <-rep(NA, N) 
      
      # confidence interval
      wald_l_month3_EY <- rep(NA, N)
      wald_r_month3_EY <- rep(NA, N)
      wald_l_month3_Bt <- rep(NA, N)
      wald_r_month3_Bt <- rep(NA, N)
      
      wald_l_month3_boot_EY <- rep(NA, N)
      wald_r_month3_boot_EY <- rep(NA, N)
      wald_l_month3_boot_Bt <- rep(NA, N)
      wald_r_month3_boot_Bt <- rep(NA, N)
      
      # wald_l_month3_boot_EY_perct <- rep(NA, N)
      # wald_r_month3_boot_EY_perct <- rep(NA, N)
      # wald_l_month3_boot_Bt_perct <- rep(NA, N)
      # wald_r_month3_boot_Bt_perct <- rep(NA, N)
      # 
      # wald_l_month3_boot_EY_bias <- rep(NA, N)
      # wald_r_month3_boot_EY_bias <- rep(NA, N)
      # wald_l_month3_boot_Bt_bias <- rep(NA, N)
      # wald_r_month3_boot_Bt_bias <- rep(NA, N)
      
      wald_l_month3_boot_EY_se <- rep(NA, N)
      wald_r_month3_boot_EY_se <- rep(NA, N)
      wald_l_month3_boot_Bt_se <- rep(NA, N)
      wald_r_month3_boot_Bt_se <- rep(NA, N)
      
      # for coverage
      Covg_wald_EY <- rep(NA, N)
      Covg_wald_Bt <- rep(NA, N)
      Covg_wald_boot_EY <- rep(NA, N)
      Covg_wald_boot_Bt <- rep(NA, N)
      # Covg_wald_boot_EY_perct <- rep(NA, N)
      # Covg_wald_boot_Bt_perct <- rep(NA, N)
      # Covg_wald_boot_EY_bias <- rep(NA, N)
      # Covg_wald_boot_Bt_bias <- rep(NA, N)
      Covg_wald_boot_EY_se <- rep(NA, N)
      Covg_wald_boot_Bt_se <- rep(NA, N)
      
      for(i in 1:N){
        
        Month_3[i] <- orig_treat[[i]][[1]][index_alpha, index_month]
        Var_3[i] <- orig_treat[[i]][[2]][index_alpha, index_month]
        
        # load the parametric bootstrap estimation results
        temp0 <- try(temp_paraboot <- readRDS(paste0("para_boot/sim_paraboot_trt_identity_200_", i)) , silent = T)
        
        if(length(temp0) > 1){
          
          ## E[Y(t)] scale
          # Wald confidence interval - on E[Y(t)] scale 
          temp1 <- wald(Month_3[i], Var_3[i])
          wald_l_month3_EY[i]  <- temp1[1]
          wald_r_month3_EY[i]  <- temp1[2]
          
          ### B(t)*beta scale
          # Wald confidence interval - trans
          temp3 <- wald(log(Month_3[i]), Var_3[i] / (Month_3[i])^2 )
          wald_l_month3_Bt[i]  <- exp(temp3[1])
          wald_r_month3_Bt[i]  <- exp(temp3[2])
          
          Est_paraboot_EY <- rep(NA, boot_size)
          Est_paraboot_Bt <- rep(NA, boot_size)
          Var_paraboot_EY <- rep(NA, boot_size)
          for(j in 1:boot_size){
            Est_paraboot_EY[j] <- temp_paraboot[[j]][[1]][index_alpha, index_month]
            Est_paraboot_Bt[j] <- log(temp_paraboot[[j]][[1]][index_alpha, index_month])
            Var_paraboot_EY[j] <- temp_paraboot[[j]][[2]][index_alpha, index_month]
          }
          
          ## boot-t 
          boot_t_EY <- (Month_3[i] -  Est_paraboot_EY) / sqrt(Var_paraboot_EY)
          quant_EY <- quantile(abs(boot_t_EY), probs = c(0.95))
          wald_l_month3_boot_EY[i] <- Month_3[i] - quant_EY[1] * sqrt(Var_3[i])
          wald_r_month3_boot_EY[i] <- Month_3[i] + quant_EY[1] * sqrt(Var_3[i])
          
          boot_t_Bt <- (log(Month_3[i]) - Est_paraboot_Bt) / sqrt(Var_paraboot_EY / Est_paraboot_EY^2)
          quant_Bt <- quantile(abs(boot_t_Bt), probs = c(0.95))
          wald_l_month3_boot_Bt[i] <- exp(log(Month_3[i]) - quant_Bt[1] * sqrt(Var_3[i] / (Month_3[i])^2))
          wald_r_month3_boot_Bt[i] <- exp(log(Month_3[i]) + quant_Bt[1] * sqrt(Var_3[i] / (Month_3[i])^2))
          
          
          # boot_t_EY <- (Month_3[i] -  Est_paraboot_EY) / sqrt(Var_paraboot_EY)
          # quant_EY <- quantile(boot_t_EY, probs = c(0.025, 0.975))
          # wald_l_month3_boot_EY[i] <- Month_3[i] + quant_EY[1] * sqrt(Var_3[i])
          # wald_r_month3_boot_EY[i] <- Month_3[i] + quant_EY[2] * sqrt(Var_3[i])
          # 
          # boot_t_Bt <- (logit(Month_3[i]) - Est_paraboot_Bt) / sqrt(Var_paraboot_EY / (Est_paraboot_EY * (1 - Est_paraboot_EY))^2)
          # quant_Bt <- quantile(boot_t_Bt, probs = c(0.025, 0.975))
          # wald_l_month3_boot_Bt[i] <- expit(logit(Month_3[i]) + quant_Bt[1] * sqrt(Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2))
          # wald_r_month3_boot_Bt[i] <- expit(logit(Month_3[i]) + quant_Bt[2] * sqrt(Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2))
          
          
          # # quantile 
          # quant_perct <- quantile(Est_paraboot_EY, probs = c(0.025, 0.975))
          # wald_l_month3_boot_EY_perct[i] <- quant_perct[1]
          # wald_r_month3_boot_EY_perct[i] <- quant_perct[2]
          # 
          # quant_perct_Bt <- quantile(Est_paraboot_Bt, probs = c(0.025, 0.975))
          # wald_l_month3_boot_Bt_perct[i] <- expit(quant_perct_Bt[1])
          # wald_r_month3_boot_Bt_perct[i] <- expit(quant_perct_Bt[2])
          # 
          # # bias-correction CI
          # quant_bias <- quantile(abs(Est_paraboot_EY - Month_3[i]), probs = c(0.95))
          # wald_l_month3_boot_EY_bias[i] <- Month_3[i] - quant_bias[1]
          # wald_r_month3_boot_EY_bias[i] <- Month_3[i] + quant_bias[1]
          # 
          # quant_bias_Bt <- quantile(abs(Est_paraboot_Bt - logit(Month_3[i])), probs = c(0.95))
          # wald_l_month3_boot_Bt_bias[i] <- expit(logit(Month_3[i]) - quant_bias_Bt[1])
          # wald_r_month3_boot_Bt_bias[i] <- expit(logit(Month_3[i]) + quant_bias_Bt[1])
          
          # normal CI using bootstrap se
          temp4_EY <- wald(Month_3[i], var(Est_paraboot_EY))
          wald_l_month3_boot_EY_se[i]  <- temp4_EY[1]
          wald_r_month3_boot_EY_se[i]  <- temp4_EY[2]
          
          temp4_Bt <- wald(log(Month_3[i]), var(Est_paraboot_Bt))
          wald_l_month3_boot_Bt_se[i]  <- exp(temp4_Bt[1])
          wald_r_month3_boot_Bt_se[i]  <- exp(temp4_Bt[2])
          
        }
        
        # Coverage
        Covg_wald_EY[i] <- ifelse((wald_l_month3_EY[i] < true_value) & (true_value < wald_r_month3_EY[i]), 1, 0)
        Covg_wald_Bt[i] <- ifelse((wald_l_month3_Bt[i] < true_value) & (true_value < wald_r_month3_Bt[i]), 1, 0)
        Covg_wald_boot_EY[i] <- ifelse((wald_l_month3_boot_EY[i] < true_value) & (true_value < wald_r_month3_boot_EY[i]), 1, 0)
        Covg_wald_boot_Bt[i] <- ifelse((wald_l_month3_boot_Bt[i] < true_value) & (true_value < wald_r_month3_boot_Bt[i]), 1, 0)
        # Covg_wald_boot_EY_perct[i] <- ifelse((wald_l_month3_boot_EY_perct[i] < true_value) & (true_value < wald_r_month3_boot_EY_perct[i]), 1, 0)
        # Covg_wald_boot_Bt_perct[i] <- ifelse((wald_l_month3_boot_Bt_perct[i] < true_value) & (true_value < wald_r_month3_boot_Bt_perct[i]), 1, 0)
        # Covg_wald_boot_EY_bias[i] <- ifelse((wald_l_month3_boot_EY_bias[i] < true_value) & (true_value < wald_r_month3_boot_EY_bias[i]), 1, 0)
        # Covg_wald_boot_Bt_bias[i] <- ifelse((wald_l_month3_boot_Bt_bias[i] < true_value) & (true_value < wald_r_month3_boot_Bt_bias[i]), 1, 0)
        Covg_wald_boot_EY_se[i] <- ifelse((wald_l_month3_boot_EY_se[i] < true_value) & (true_value < wald_r_month3_boot_EY_se[i]), 1, 0)
        Covg_wald_boot_Bt_se[i] <- ifelse((wald_l_month3_boot_Bt_se[i] < true_value) & (true_value < wald_r_month3_boot_Bt_se[i]), 1, 0)
        
      }
      
      #print(length(which(is.na(Covg_wald_boot_EY) == F)))
      
      result <- data.frame(Emp_mean = mean(Month_3), 
                           bias = mean(Month_3) - true_value,
                           Wald_EY = mean(Covg_wald_EY, na.rm = T), 
                           Wald_Bt = mean(Covg_wald_Bt, na.rm = T), 
                           Wald_Bt_paraboot_EY = mean(Covg_wald_boot_EY, na.rm = T),
                           Wald_Bt_paraboot_Bt = mean(Covg_wald_boot_Bt, na.rm = T),
                           # Wald_Bt_paraboot_EY_perct = mean(Covg_wald_boot_EY_perct, na.rm = T),
                           # Wald_Bt_paraboot_Bt_perct = mean(Covg_wald_boot_Bt_perct, na.rm = T),
                           # Wald_Bt_paraboot_EY_bias = mean(Covg_wald_boot_EY_bias, na.rm = T),
                           # Wald_Bt_paraboot_Bt_bias = mean(Covg_wald_boot_Bt_bias, na.rm = T),
                           Wald_Bt_paraboot_EY_se = mean(Covg_wald_boot_EY_se, na.rm = T),
                           Wald_Bt_paraboot_Bt_se = mean(Covg_wald_boot_Bt_se, na.rm = T)
      )
      
      result
      
      return(result)
    }
    
    f_res_point <- function(alpha = -0.6){
      
      result <- c()
      for(month_i in 1:4){
        
        month <- Month[month_i]
        
        temp <- get(paste0("month", month, "_treat"))
        index_alpha <- which(Alpha_seq == alpha)
        true_value <- temp[index_alpha]
        
        result1 <- f_point(alpha, month, true_value)
        result <- rbind(result, c(alpha, true_value, round(result1, 3)))
        
      }
      colnames(result)[1:2] <- c("alpha", "true_value")
      result
      
    }
    
  }
  
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/identity_simu/simu_treat")
  sample_size <- 200
  boot_size <- 200
  N <- 500
  
  # identity link
  orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_1e7_200_GL_30")
  f_res_point(alpha = -0.6)
  f_res_point(alpha = -0.3)
  f_res_point(alpha = 0)
  f_res_point(alpha = 0.3)
  f_res_point(alpha = 0.6)
  
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/identity_simu/simu_treat")
  miss_index <- c()
  for(i in 1:N){
    temp0 <- try(temp_paraboot <- readRDS(paste0("para_boot/sim_paraboot_trt_identity_200_", i)) , silent = T)
    if(length(temp0) <= 1){
      miss_index <- c(miss_index, i)
    }
  }
  miss_index
  length(miss_index)
  
}


