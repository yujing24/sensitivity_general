rm(list=ls())
library(dplyr)


#### parametric bootstrap
{
  
  
  ### setting
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
    
    
    # single index model +  original
    # month <- 3
    month3_treat <- c(0.171, 0.228, 0.297, 0.375, 0.462, 0.552)
    month3_control <- c(0.152, 0.202, 0.262, 0.333, 0.412, 0.497)

    # month <- 6
    month6_treat <- c(0.186, 0.245, 0.316, 0.394, 0.481, 0.57)
    month6_control <- c(0.158, 0.209, 0.271, 0.343, 0.422, 0.506)

    # month <- 9
    month9_treat <- c(0.195, 0.256, 0.328, 0.407, 0.493, 0.58)
    month9_control <- c(0.162, 0.214, 0.277, 0.348, 0.428, 0.513)

    # month <- 12
    month12_treat <- c(0.2, 0.263, 0.334, 0.414, 0.499, 0.586)
    month12_control <- c(0.164, 0.216, 0.28, 0.352, 0.433, 0.517)

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
    
    # # logistic + logit link
    # month3_treat <- c(0.17528, 0.23253, 0.29979, 0.37766, 0.46339, 0.55269)
    # month6_treat <- c(0.19083, 0.25123, 0.31990, 0.39857, 0.48426, 0.57133)
    # month9_treat <- c(0.19837, 0.26035, 0.33039, 0.40918, 0.49422, 0.58054)
    # month12_treat <- c(0.20356, 0.26597, 0.33736, 0.41639, 0.50118, 0.58660)
    # 
    # month3_control <- c(0.13761, 0.18518, 0.24377, 0.31265, 0.39025, 0.47700)
    # month6_control <- c(0.14696, 0.19688, 0.25781, 0.32844, 0.40780, 0.49500)
    # month9_control <- c(0.15750, 0.20958, 0.27300, 0.34542, 0.42643, 0.51300)
    # month12_control <- c(0.17020, 0.22499, 0.29071, 0.36504, 0.44775, 0.53400)
    
    
    Alpha_seq <- c(-0.8, -0.4, 0, 0.4, 0.8, 1.2)
    Month <- c(3, 6, 9, 12)
    N <- 500
    
    wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
      CI_left <- mean_orig - Zstar * sqrt(var_orig)
      CI_right <- mean_orig + Zstar * sqrt(var_orig)
      
      return(c(CI_left, CI_right))
    }
    
    logit <- function(x){
      return(log( x / (1-x)))
    }
    expit <- function(x){
      return(1 / (1 + exp(-x)))
    }
    
    f_point <- function(alpha, month, true_value,
                        weight_type = "weight2"){
      
      
      Alpha_seq <- c(-0.8, -0.4, 0, 0.4, 0.8, 1.2)
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
      
      wald_l_month3_boot_EY_se <- rep(NA, N)
      wald_r_month3_boot_EY_se <- rep(NA, N)
      wald_l_month3_boot_Bt_se <- rep(NA, N)
      wald_r_month3_boot_Bt_se <- rep(NA, N)
      
      # for coverage
      Covg_wald_EY <- rep(NA, N)
      Covg_wald_Bt <- rep(NA, N)
      Covg_wald_boot_EY <- rep(NA, N)
      Covg_wald_boot_Bt <- rep(NA, N)
      Covg_wald_boot_EY_se <- rep(NA, N)
      Covg_wald_boot_Bt_se <- rep(NA, N)
      
      for(i in 1:N){
        
        Month_3[i] <- orig_treat[[i]][[1]][index_alpha, index_month]
        Var_3[i] <- orig_treat[[i]][[2]][index_alpha, index_month]
        
        # load the parametric bootstrap estimation results
        if(weight_type == "weight1"){
          temp0 <- try(temp_paraboot <- readRDS(paste0("para_boot_weight1/sim_paraboot_trt_logit_200_", i)) , silent = T)
        }else{
          temp0 <- try(temp_paraboot <- readRDS(paste0("para_boot_weight2/sim_paraboot_trt_logit_200_", i)) , silent = T)
        }
        
        if(length(temp0) > 1){
          
          ## E[Y(t)] scale
          # Wald confidence interval - on E[Y(t)] scale 
          temp1 <- wald(Month_3[i], Var_3[i])
          wald_l_month3_EY[i]  <- temp1[1]
          wald_r_month3_EY[i]  <- temp1[2]
          
          ### B(t)*beta scale
          # Wald confidence interval - trans
          temp3 <- wald(logit(Month_3[i]), Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2 )
          wald_l_month3_Bt[i]  <- expit(temp3[1])
          wald_r_month3_Bt[i]  <- expit(temp3[2])
          
          Est_paraboot_EY <- rep(NA, boot_size)
          Est_paraboot_Bt <- rep(NA, boot_size)
          Var_paraboot_EY <- rep(NA, boot_size)
          for(j in 1:boot_size){
            Est_paraboot_EY[j] <- temp_paraboot[[j]][[1]][index_alpha, index_month]
            Est_paraboot_Bt[j] <- logit(temp_paraboot[[j]][[1]][index_alpha, index_month])
            Var_paraboot_EY[j] <- temp_paraboot[[j]][[2]][index_alpha, index_month]
          }
          
          ## boot-t 
          boot_t_EY <- (Month_3[i] -  Est_paraboot_EY) / sqrt(Var_paraboot_EY)
          quant_EY <- quantile(abs(boot_t_EY), probs = c(0.95))
          wald_l_month3_boot_EY[i] <- Month_3[i] - quant_EY[1] * sqrt(Var_3[i])
          wald_r_month3_boot_EY[i] <- Month_3[i] + quant_EY[1] * sqrt(Var_3[i])
          
          boot_t_Bt <- (logit(Month_3[i]) - Est_paraboot_Bt) / sqrt(Var_paraboot_EY / (Est_paraboot_EY * (1 - Est_paraboot_EY))^2)
          quant_Bt <- quantile(abs(boot_t_Bt), probs = c(0.95))
          wald_l_month3_boot_Bt[i] <- expit(logit(Month_3[i]) - quant_Bt[1] * sqrt(Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2))
          wald_r_month3_boot_Bt[i] <- expit(logit(Month_3[i]) + quant_Bt[1] * sqrt(Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2))
          
          # normal CI using bootstrap se
          temp4_EY <- wald(Month_3[i], var(Est_paraboot_EY))
          wald_l_month3_boot_EY_se[i]  <- temp4_EY[1]
          wald_r_month3_boot_EY_se[i]  <- temp4_EY[2]
          
          temp4_Bt <- wald(logit(Month_3[i]), var(Est_paraboot_Bt))
          wald_l_month3_boot_Bt_se[i]  <- expit(temp4_Bt[1])
          wald_r_month3_boot_Bt_se[i]  <- expit(temp4_Bt[2])
          
        }
        
        # Coverage
        Covg_wald_EY[i] <- ifelse((wald_l_month3_EY[i] < true_value) & (true_value < wald_r_month3_EY[i]), 1, 0)
        Covg_wald_Bt[i] <- ifelse((wald_l_month3_Bt[i] < true_value) & (true_value < wald_r_month3_Bt[i]), 1, 0)
        Covg_wald_boot_EY[i] <- ifelse((wald_l_month3_boot_EY[i] < true_value) & (true_value < wald_r_month3_boot_EY[i]), 1, 0)
        Covg_wald_boot_Bt[i] <- ifelse((wald_l_month3_boot_Bt[i] < true_value) & (true_value < wald_r_month3_boot_Bt[i]), 1, 0)
        Covg_wald_boot_EY_se[i] <- ifelse((wald_l_month3_boot_EY_se[i] < true_value) & (true_value < wald_r_month3_boot_EY_se[i]), 1, 0)
        Covg_wald_boot_Bt_se[i] <- ifelse((wald_l_month3_boot_Bt_se[i] < true_value) & (true_value < wald_r_month3_boot_Bt_se[i]), 1, 0)
        
      }
      
      # print(length(which(is.na(Covg_wald_boot_EY) == F)))
      
      result <- data.frame(Emp_mean = mean(Month_3), 
                           bias = mean(Month_3) - true_value,
                           Wald_EY = mean(Covg_wald_EY, na.rm = T),
                           Wald_Bt = mean(Covg_wald_Bt, na.rm = T),
                           Wald_Bt_paraboot_EY = mean(Covg_wald_boot_EY, na.rm = T),
                           Wald_Bt_paraboot_Bt = mean(Covg_wald_boot_Bt, na.rm = T),
                           Wald_Bt_paraboot_EY_se = mean(Covg_wald_boot_EY_se, na.rm = T),
                           Wald_Bt_paraboot_Bt_se = mean(Covg_wald_boot_Bt_se, na.rm = T) )
      
      result
      
      
      return(result)
    }
    
    f_res_point <- function(alpha = -0.6,
                            weight_type){
      
      result <- c()
      for(month_i in 1:4){
        
        month <- Month[month_i]
        
        temp <- get(paste0("month", month, "_treat"))
        index_alpha <- which(Alpha_seq == alpha)
        true_value <- temp[index_alpha]
        
        result1 <- f_point(alpha, month, true_value, weight_type)
        result <- rbind(result, c(alpha, true_value, round(result1, 3)))
        
      }
      colnames(result)[1:2] <- c("alpha", "true_value")
      result
      
    }
    
  }
  
  # logistic + norm1
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/logistic_simu_logit_norm1/simu_treat")
    sample_size <- 200
    boot_size <- 200
    N <- 500
    
    # weight 1
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_1e7_200_GL_30_norm1")
    f_res_point(alpha = -0.8, weight_type = "weight1")
    f_res_point(alpha = -0.4, weight_type = "weight1")
    f_res_point(alpha = 0, weight_type = "weight1")
    f_res_point(alpha = 0.4, weight_type = "weight1")
    f_res_point(alpha = 0.8, weight_type = "weight1")
    f_res_point(alpha = 1.2, weight_type = "weight1")
    
    
    # weight 1
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight2_1e7_200_GL_30_norm1")
    f_res_point(alpha = -0.8, weight_type = "weight2")
    f_res_point(alpha = -0.4, weight_type = "weight2")
    f_res_point(alpha = 0, weight_type = "weight2")
    f_res_point(alpha = 0.4, weight_type = "weight2")
    f_res_point(alpha = 0.8, weight_type = "weight2")
    f_res_point(alpha = 1.2, weight_type = "weight2")
    
    
  }
  
  
  # single index model +  norm1
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/single_simu_logit_norm1/simu_treat")
    sample_size <- 200
    boot_size <- 200
    N <- 500
    
    # weight 1
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_1e7_200_GL_30_norm1")
    f_res_point(alpha = -0.8, weight_type = "weight1")
    f_res_point(alpha = -0.4, weight_type = "weight1")
    f_res_point(alpha = 0, weight_type = "weight1")
    f_res_point(alpha = 0.4, weight_type = "weight1")
    f_res_point(alpha = 0.8, weight_type = "weight1")
    f_res_point(alpha = 1.2, weight_type = "weight1")
    
    
    # weight 2
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight2_1e7_200_GL_30_norm1")
    f_res_point(alpha = -0.8, weight_type = "weight2")
    f_res_point(alpha = -0.4, weight_type = "weight2")
    f_res_point(alpha = 0, weight_type = "weight2")
    f_res_point(alpha = 0.4, weight_type = "weight2")
    f_res_point(alpha = 0.8, weight_type = "weight2")
    f_res_point(alpha = 1.2, weight_type = "weight2")
  }
  
  
  # single index model +  fix_coef
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/single_simu_logit_fixcoef/simu_treat")
    sample_size <- 200
    boot_size <- 200
    N <- 500
    
    # weight 1
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_1e7_200_GL_30")
    f_res_point(alpha = -0.8, weight_type = "weight1")
    f_res_point(alpha = -0.4, weight_type = "weight1")
    f_res_point(alpha = 0, weight_type = "weight1")
    f_res_point(alpha = 0.4, weight_type = "weight1")
    f_res_point(alpha = 0.8, weight_type = "weight1")
    f_res_point(alpha = 1.2, weight_type = "weight1")
    
    
    # weight 2
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight2_1e7_200_GL_30")
    f_res_point(alpha = -0.8, weight_type = "weight2")
    f_res_point(alpha = -0.4, weight_type = "weight2")
    f_res_point(alpha = 0, weight_type = "weight2")
    f_res_point(alpha = 0.4, weight_type = "weight2")
    f_res_point(alpha = 0.8, weight_type = "weight2")
    f_res_point(alpha = 1.2, weight_type = "weight2")
  }
  
  
  # logistic +  fix_coef
  {
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/logistic_simu_logit_fixcoef/simu_treat")
    sample_size <- 200
    boot_size <- 200
    N <- 500
    
    # weight 1
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_1e7_200_GL_30")
    f_res_point(alpha = -0.8, weight_type = "weight1")
    f_res_point(alpha = -0.4, weight_type = "weight1")
    f_res_point(alpha = 0, weight_type = "weight1")
    f_res_point(alpha = 0.4, weight_type = "weight1")
    f_res_point(alpha = 0.8, weight_type = "weight1")
    f_res_point(alpha = 1.2, weight_type = "weight1")
    
    
    # weight 2
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight2_1e7_200_GL_30")
    f_res_point(alpha = -0.8, weight_type = "weight2")
    f_res_point(alpha = -0.4, weight_type = "weight2")
    f_res_point(alpha = 0, weight_type = "weight2")
    f_res_point(alpha = 0.4, weight_type = "weight2")
    f_res_point(alpha = 0.8, weight_type = "weight2")
    f_res_point(alpha = 1.2, weight_type = "weight2")
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
    
    
    Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
    Month <- c(3, 6, 9, 12)
    
    N_size <- 200
    
    wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
      CI_left <- mean_orig - Zstar * sqrt(var_orig)
      CI_right <- mean_orig + Zstar * sqrt(var_orig)
      
      return(c(CI_left, CI_right))
    }
    
    logit <- function(x){
      return(log( x / (1-x)))
    }
    expit <- function(x){
      return(1 / (1 + exp(-x)))
    }
    
    f <- function(alpha, month, true_value, 
                  weight_type = "weight1",
                  N = 500){
      
      Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
      index_alpha <- which(Alpha_seq == alpha)
      index_month <- which(c(1, 3, 6, 9, 12) == month)
      
      # for original data's estimation
      Month_3 <- rep(NA, N)
      Var_3 <-rep(NA, N) 
      
      # for confidence interval
      wald_l_month3 <- rep(NA, N)
      wald_r_month3 <- rep(NA, N)
      wald_l_month3_trans <- rep(NA, N)
      wald_r_month3_trans <- rep(NA, N)
      
      JN_l_month3 <- rep(NA, N)
      JN_r_month3 <- rep(NA, N)
      JN_l_month3_trans <- rep(NA, N)
      JN_r_month3_trans <- rep(NA, N)
      # JN_l_month3_perc <- rep(NA, N)
      # JN_r_month3_perc <- rep(NA, N)
      # JN_l_month3_seudo <- rep(NA, N)
      # JN_r_month3_seudo <- rep(NA, N)
      
      # for coverage
      Covg_wald <- rep(NA, N)
      Covg_wald_trans <- rep(NA, N)
      
      Covg_JN <- rep(NA, N)
      Covg_JN_trans <- rep(NA, N)
      # Covg_JN_perc <- rep(NA, N)
      # Covg_JN_seudo <- rep(NA, N)
      
      for(i in 1:N){
        
        Month_3[i] <- orig_treat[[i]][[1]][index_alpha, index_month]
        Var_3[i] <- orig_treat[[i]][[2]][index_alpha, index_month]
        
        # Wald confidence interval 
        temp2 <- wald(Month_3[i], Var_3[i])
        wald_l_month3[i]  <- temp2[1]
        wald_r_month3[i]  <- temp2[2]
        
        # Wald confidence interval - trans
        temp3 <- wald(logit(Month_3[i]), Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2 )
        wald_l_month3_trans[i]  <- expit(temp3[1])
        wald_r_month3_trans[i]  <- expit(temp3[2])
        
        # Coverage
        Covg_wald[i] <- ifelse((wald_l_month3[i] < true_value) & (true_value < wald_r_month3[i]), 1, 0)
        Covg_wald_trans[i] <- ifelse((wald_l_month3_trans[i] < true_value) & (true_value < wald_r_month3_trans[i]), 1, 0)
        
        if(i %in% data_list_JN){
          
          # load the JN variance
          if(weight_type == "weight1"){
            temp_JN <- readRDS(paste0("JN_simu_treat_weight1_1e7/simu_JN_", i, "_trt")) 
          }else{
            temp_JN <- readRDS(paste0("JN_simu_treat_weight2_1e7/simu_JN_", i, "_trt")) 
          }
          
          Est_JN <- rep(NA, N_size)
          Est_JN_trans <- rep(NA, N_size)
          for(j in 1:N_size){
            temp <- try(temp_JN[[j]][[1]], silent = T)
            if(length(temp) > 1){
              Est_JN[j] <- temp_JN[[j]][[1]][index_alpha, index_month]
              Est_JN_trans[j] <- logit(temp_JN[[j]][[1]][index_alpha, index_month])
            }
          }
          
          if(sum(is.na(Est_JN)) > 0){
            
            Covg_JN[i] <- NA
            Covg_JN_trans[i] <- NA
            # Covg_JN_perc[i] <- NA
            # Covg_JN_seudo[i] <- NA
            
          }else{
            
            var_JN <- (N_size - 1)/N_size * sum((Est_JN - mean(Est_JN))^2)
            var_JN_trans <- (N_size - 1)/N_size * sum((Est_JN_trans - mean(Est_JN_trans))^2)
            
            # JN confidence interval - original scale
            temp4 <- wald(Month_3[i], var_JN)
            JN_l_month3[i] <- temp4[1]
            JN_r_month3[i] <- temp4[2]
            
            # JN confidence interval - transformed scale
            temp5 <- wald(logit(Month_3[i]), var_JN_trans)
            JN_l_month3_trans[i] <- expit(temp5[1])
            JN_r_month3_trans[i] <- expit(temp5[2])
            
            # # JN percentile interval 
            # JN_l_month3_perc[i] <- quantile(Est_JN, probs = 0.025)
            # JN_r_month3_perc[i] <- quantile(Est_JN, probs = 0.975)
            
            # # Pseudovalue jackknife CI
            # Est_seudo <- N_size * Month_3[i] - (N_size - 1) * Est_JN
            # Month_jack_seudo <- mean(Est_seudo)
            # var_jack_seudo <- 1 / (N_size * (N_size - 1)) * sum((Est_seudo - Month_jack_seudo)^2)
            # temp6 <- wald(Month_jack_seudo, var_jack_seudo)
            # JN_l_month3_seudo[i] <- temp6[1]
            # JN_r_month3_seudo[i] <- temp6[2]
            
            # Coverage
            Covg_JN[i] <-  ifelse((JN_l_month3[i] < true_value) & (true_value < JN_r_month3[i]), 1, 0)
            Covg_JN_trans[i] <-  ifelse((JN_l_month3_trans[i] < true_value) & (true_value < JN_r_month3_trans[i]), 1, 0)
            # Covg_JN_perc[i] <-  ifelse((JN_l_month3_perc[i] < true_value) & (true_value < JN_r_month3_perc[i]), 1, 0)
            # Covg_JN_seudo[i] <-  ifelse((JN_l_month3_seudo[i] < true_value) & (true_value < JN_r_month3_seudo[i]), 1, 0)
            
          }
          
        }
        
      }
      
      result <- data.frame(Emp_mean = mean(Month_3), 
                           bias = mean(Month_3) - true_value,
                           Wald = mean(Covg_wald),
                           Wald_trans = mean(Covg_wald_trans),
                           JN = mean(Covg_JN, na.rm = T),
                           JN_trans = mean(Covg_JN_trans, na.rm = T),
                           NA_JN_num = sum(is.na(Covg_JN))
                           # JN_perc = mean(Covg_JN_perc, na.rm = T),
                           # JN_seudo = mean(Covg_JN_seudo, na.rm = T)
      )
      
      return(result)
    }
    
    f_res <- function(alpha = -0.6,
                      weight_type){
      
      result <- c()
      for(month_i in 1:4){
        
        month <- Month[month_i]
        
        temp <- get(paste0("month", month, "_treat"))
        index_alpha <- which(Alpha_seq == alpha)
        true_value <- temp[index_alpha]
        
        result1 <- f(alpha, month, true_value, 
                     weight_type)
        result <- rbind(result, c(true_value, round(result1, 3)))
        
      }
      result
      
    }
    
    f_point <- function(alpha, month, true_value,
                        weight_type = "weight1",
                        N = 500){
      
      Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
      index_alpha <- which(Alpha_seq == alpha)
      index_month <- which(c(1, 3, 6, 9, 12) == month)
      
      # for original data's estimation
      Month_3 <- rep(NA, N)
      Var_3 <-rep(NA, N) 
      
      # for confidence interval
      wald_l_month3 <- rep(NA, N)
      wald_r_month3 <- rep(NA, N)
      wald_l_month3_trans <- rep(NA, N)
      wald_r_month3_trans <- rep(NA, N)
      
      # for coverage
      Covg_wald <- rep(NA, N)
      Covg_wald_trans <- rep(NA, N)
      
      for(i in 1:N){
        
        Month_3[i] <- orig_treat[[i]][[1]][index_alpha, index_month]
        Var_3[i] <- orig_treat[[i]][[2]][index_alpha, index_month]
        
        # Wald confidence interval 
        temp2 <- wald(Month_3[i], Var_3[i])
        wald_l_month3[i]  <- temp2[1]
        wald_r_month3[i]  <- temp2[2]
        
        # Wald confidence interval - trans
        temp3 <- wald(logit(Month_3[i]), Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2 )
        wald_l_month3_trans[i]  <- expit(temp3[1])
        wald_r_month3_trans[i]  <- expit(temp3[2])
        
        # Coverage
        Covg_wald[i] <- ifelse((wald_l_month3[i] < true_value) & (true_value < wald_r_month3[i]), 1, 0)
        Covg_wald_trans[i] <- ifelse((wald_l_month3_trans[i] < true_value) & (true_value < wald_r_month3_trans[i]), 1, 0)
      }
      
      result <- data.frame(Emp_mean = mean(Month_3), 
                           bias = mean(Month_3) - true_value,
                           Wald = mean(Covg_wald), 
                           Wald_trans = mean(Covg_wald_trans))
      return(result)
    }
    
    f_res_point <- function(alpha = -0.6,
                            weight_type){
      
      result <- c()
      for(month_i in 1:4){
        
        month <- Month[month_i]
        
        temp <- get(paste0("month", month, "_treat"))
        index_alpha <- which(Alpha_seq == alpha)
        true_value <- temp[index_alpha]
        
        result1 <- f_point(alpha, month, true_value, weight_type)
        result <- rbind(result, c(true_value, round(result1, 3)))
        
      }
      result
      
    }
    
  }
  
  
  # logit link - weight1
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/binary_simu/simu_treat")
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_1e7")
    # orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_1e9")
    # orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_norm1") # not well
    # orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_band_1e7")
    # orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_band_1e9") # similar as fix coef
    
    
    # no jackknife, only point estimation
    f_res_point(alpha = -0.6, weight_type = "weight1")
    f_res_point(alpha = -0.3, weight_type = "weight1")
    f_res_point(alpha = 0, weight_type = "weight1")
    f_res_point(alpha = 0.3, weight_type = "weight1")
    f_res_point(alpha = 0.6, weight_type = "weight1")
    
    
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_coef_1e7_400")
    f_res_point(alpha = -0.6, weight_type = "weight1")
    f_res_point(alpha = -0.3, weight_type = "weight1")
    f_res_point(alpha = 0, weight_type = "weight1")
    f_res_point(alpha = 0.3, weight_type = "weight1")
    f_res_point(alpha = 0.6, weight_type = "weight1")
    
    
    data_list <- c(1:500)
    data_list_JN <- c() 
    for(k in 1:500){
      temp_data <- try( readRDS(paste0("JN_simu_treat_weight1_norm1/simu_JN_", k, "_trt")) , silent = TRUE)
      if(length(temp_data) > 1){
        data_list_JN <- c(data_list_JN, k)
      }
    }
    length(data_list_JN)
    setdiff(data_list, data_list_JN)
    
    # including jackknife 
    f_res(alpha = -0.6, weight_type = "weight1")
    f_res(alpha = -0.3, weight_type = "weight1")
    f_res(alpha = 0, weight_type = "weight1")
    f_res(alpha = 0.3, weight_type = "weight1")
    f_res(alpha = 0.6, weight_type = "weight1")
    
  }
  
  
  # logit link - weight2
  {
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/binary_simu/simu_treat")
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight2_1e7")
    
    
    data_list <- c(1:500)
    data_list_JN <- c() 
    for(k in 1:500){
      temp_data <- try( readRDS(paste0("JN_simu_treat_weight2_1e7/simu_JN_", k, "_trt")) , silent = TRUE)
      if(length(temp_data) > 1){
        data_list_JN <- c(data_list_JN, k)
      }
    }
    length(data_list_JN)
    setdiff(data_list, data_list_JN)
    
    
    # including jackknife 
    f_res(alpha = -0.6, weight_type = "weight2")
    f_res(alpha = -0.3, weight_type = "weight2")
    f_res(alpha = 0, weight_type = "weight2")
    f_res(alpha = 0.3, weight_type = "weight2")
    f_res(alpha = 0.6, weight_type = "weight2")
    
  }
  
  
  #### Edgeworth expansion for skewness-corrected CI
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
      
      
      Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
      Month <- c(3, 6, 9, 12)
      N <- 500
      
      wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
        CI_left <- mean_orig - Zstar * sqrt(var_orig)
        CI_right <- mean_orig + Zstar * sqrt(var_orig)
        
        return(c(CI_left, CI_right))
      }
      
      wald_skew <- function(mean_orig, var_orig, Zstar = qnorm(0.975), skew){
        
        # Z_star1 <- Zstar +  skew / 6 / sqrt(sample_size) * ( 2 * Zstar^2 + 1)
        # Z_star2 <- Zstar - skew / 6 / sqrt(sample_size) * ( 2 * Zstar^2 + 1)
        
        # another method from one paper
        # CI_left <- mean_orig - sqrt(var_orig) * (Zstar + skew / 6 / sqrt(sample_size) * (2 * Zstar^2 + 1))
        # CI_right <- mean_orig + sqrt(var_orig) * (Zstar - skew / 6 / sqrt(sample_size) * (2 * Zstar^2 + 1))
        
        # Edgeworth expansion
        CI_left <- mean_orig - sqrt(var_orig) * (Zstar + skew / 6 / sqrt(sample_size) * (Zstar^2 - 1))
        CI_right <- mean_orig - sqrt(var_orig) * (-Zstar + skew / 6 / sqrt(sample_size) * (Zstar^2 + 1))
        
        
        return(c(CI_left, CI_right))
      }
      wald_skew_peter <- function(mean_orig, var_orig, Zstar = qnorm(0.975), skew){
        
        # Z_star1 <- Zstar +  skew / 6 / sqrt(sample_size) * ( 2 * Zstar^2 + 1)
        # Z_star2 <- Zstar - skew / 6 / sqrt(sample_size) * ( 2 * Zstar^2 + 1)
        
        # another method from one paper - Peter Hall's paper
        CI_left <- mean_orig - sqrt(var_orig) * (Zstar + skew / 6 / sqrt(sample_size) * (2 * Zstar^2 + 1))
        CI_right <- mean_orig + sqrt(var_orig) * (Zstar - skew / 6 / sqrt(sample_size) * (2 * Zstar^2 + 1))
        
        # # Edgeworth expansion - Chatgpt
        # CI_left <- mean_orig - sqrt(var_orig) * (Zstar + skew / 6 / sqrt(sample_size) * (Zstar^2 - 1))
        # CI_right <- mean_orig - sqrt(var_orig) * (-Zstar + skew / 6 / sqrt(sample_size) * (Zstar^2 + 1))
        
        return(c(CI_left, CI_right))
      }
      
      logit <- function(x){
        return(log( x / (1-x)))
      }
      expit <- function(x){
        return(1 / (1 + exp(-x)))
      }
      
      
      f_point <- function(alpha, month, true_value,
                          weight_type = "weight1"){
        
        Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
        index_alpha <- which(Alpha_seq == alpha)
        index_month <- which(c(1, 3, 6, 9, 12) == month)
        
        # for original data's estimation
        Month_3 <- rep(NA, N)
        Var_3 <-rep(NA, N) 
        Skew_Bt <- rep(NA, N)
        Skew_EY <- rep(NA, N)
        
        # confidence interval
        wald_l_month3_EY <- rep(NA, N)
        wald_r_month3_EY <- rep(NA, N)
        wald_l_month3_Bt <- rep(NA, N)
        wald_r_month3_Bt <- rep(NA, N)
        
        # skewness-corrected CI - chatgpt
        wald_l_month3_EY_skew <- rep(NA, N)
        wald_r_month3_EY_skew <- rep(NA, N)
        wald_l_month3_Bt_skew <- rep(NA, N)
        wald_r_month3_Bt_skew <- rep(NA, N)
        
        # skewness-corrected CI - Peter
        wald_l_month3_EY_skew_peter <- rep(NA, N)
        wald_r_month3_EY_skew_peter <- rep(NA, N)
        wald_l_month3_Bt_skew_peter <- rep(NA, N)
        wald_r_month3_Bt_skew_peter <- rep(NA, N)
        
        # for coverage
        Covg_wald_EY <- rep(NA, N)
        Covg_wald_Bt <- rep(NA, N)
        Covg_wald_EY_skew <- rep(NA, N)
        Covg_wald_Bt_skew <- rep(NA, N)
        Covg_wald_EY_skew_peter <- rep(NA, N)
        Covg_wald_Bt_skew_peter <- rep(NA, N)
        
        
        for(i in 1:N){
          
          Month_3[i] <- orig_treat[[i]][[1]][index_alpha, index_month]
          Var_3[i] <- orig_treat[[i]][[2]][index_alpha, index_month]
          Skew_Bt[i] <- orig_treat[[i]][[3]][index_alpha, index_month] # on B(t) * beta scale
          Skew_EY[i] <- orig_treat[[i]][[4]][index_alpha, index_month] # on s(B(t) * beta) scale
          
          ## E[Y(t)] scale
          # Wald confidence interval - on E[Y(t)] scale 
          temp1 <- wald(Month_3[i], Var_3[i])
          wald_l_month3_EY[i]  <- temp1[1]
          wald_r_month3_EY[i]  <- temp1[2]
          
          # Wald confidence interval - on E[Y(t)] scale - skewness-corrected
          temp2 <- wald_skew(Month_3[i], Var_3[i], Zstar = qnorm(0.975), skew = Skew_EY[i])
          wald_l_month3_EY_skew[i]  <- temp2[1]
          wald_r_month3_EY_skew[i]  <- temp2[2]
          
          
          # Wald confidence interval - on E[Y(t)] scale - skewness-corrected (peter)
          temp2_peter <- wald_skew_peter(Month_3[i], Var_3[i], Zstar = qnorm(0.975), skew = Skew_EY[i])
          wald_l_month3_EY_skew_peter[i]  <- temp2_peter[1]
          wald_r_month3_EY_skew_peter[i]  <- temp2_peter[2]
          
          
          ## B(t)*beta scale (transfered)- skewness-corrected
          temp3 <- wald(logit(Month_3[i]), Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2 )
          wald_l_month3_Bt[i]  <- expit(temp3[1])
          wald_r_month3_Bt[i]  <- expit(temp3[2])
          
          # Wald confidence interval - on B(t)*beta scale - skewness-corrected
          temp4 <- wald_skew(logit(Month_3[i]), Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2, Zstar = qnorm(0.975), skew = Skew_Bt[i])
          wald_l_month3_Bt_skew[i]  <- expit(temp4[1])
          wald_r_month3_Bt_skew[i]  <- expit(temp4[2])
          
          # Wald confidence interval - on B(t)*beta scale - skewness-corrected (peter)
          temp4_peter <- wald_skew_peter(logit(Month_3[i]), Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2, Zstar = qnorm(0.975), skew = Skew_Bt[i])
          wald_l_month3_Bt_skew_peter[i]  <- expit(temp4_peter[1])
          wald_r_month3_Bt_skew_peter[i]  <- expit(temp4_peter[2])
          
          
          # Coverage
          Covg_wald_EY[i] <- ifelse((wald_l_month3_EY[i] < true_value) & (true_value < wald_r_month3_EY[i]), 1, 0)
          Covg_wald_EY_skew[i] <- ifelse((wald_l_month3_EY_skew[i] < true_value) & (true_value < wald_r_month3_EY_skew[i]), 1, 0)
          Covg_wald_EY_skew_peter[i] <- ifelse((wald_l_month3_EY_skew_peter[i] < true_value) & (true_value < wald_r_month3_EY_skew_peter[i]), 1, 0)
          
          Covg_wald_Bt[i] <- ifelse((wald_l_month3_Bt[i] < true_value) & (true_value < wald_r_month3_Bt[i]), 1, 0)
          Covg_wald_Bt_skew[i] <- ifelse((wald_l_month3_Bt_skew[i] < true_value) & (true_value < wald_r_month3_Bt_skew[i]), 1, 0)
          Covg_wald_Bt_skew_peter[i] <- ifelse((wald_l_month3_Bt_skew_peter[i] < true_value) & (true_value < wald_r_month3_Bt_skew_peter[i]), 1, 0)
          
        }
        
        result <- data.frame(Emp_mean = mean(Month_3), 
                             bias = mean(Month_3) - true_value,
                             Wald_EY = mean(Covg_wald_EY), 
                             Wald_EY_skew = mean(Covg_wald_EY_skew),
                             Wald_EY_skew_peter = mean(Covg_wald_EY_skew_peter),
                             Wald_Bt = mean(Covg_wald_Bt), 
                             Wald_Bt_skew = mean(Covg_wald_Bt_skew),
                             Wald_Bt_skew_peter = mean(Covg_wald_Bt_skew_peter)
        )
        return(result)
      }
      
      f_res_point <- function(alpha = -0.6,
                              weight_type){
        
        result <- c()
        for(month_i in 1:4){
          
          month <- Month[month_i]
          
          temp <- get(paste0("month", month, "_treat"))
          index_alpha <- which(Alpha_seq == alpha)
          true_value <- temp[index_alpha]
          
          result1 <- f_point(alpha, month, true_value, weight_type)
          result <- rbind(result, c(alpha, true_value, round(result1, 3)))
          
        }
        colnames(result)[1:2] <- c("alpha", "true_value")
        result
        
      }
      
    }
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/binary_simu/simu_treat")
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_1e7_400_skew")
    sample_size <- 400
    
    
    # no jackknife
    f_res_point(alpha = -0.6, weight_type = "weight1")
    f_res_point(alpha = -0.3, weight_type = "weight1")
    f_res_point(alpha = 0, weight_type = "weight1")
    f_res_point(alpha = 0.3, weight_type = "weight1")
    f_res_point(alpha = 0.6, weight_type = "weight1")
    
    
  }
  
  
  #### multiplier bootstrap
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
      
      
      Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
      Month <- c(3, 6, 9, 12)
      N <- 500
      
      wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
        CI_left <- mean_orig - Zstar * sqrt(var_orig)
        CI_right <- mean_orig + Zstar * sqrt(var_orig)
        
        return(c(CI_left, CI_right))
      }
      
      logit <- function(x){
        return(log( x / (1-x)))
      }
      expit <- function(x){
        return(1 / (1 + exp(-x)))
      }
      
      f_point <- function(alpha, month, true_value,
                          weight_type = "weight1"){
        
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
        wald_l_month3_boot_EY2 <- rep(NA, N)
        wald_r_month3_boot_EY2 <- rep(NA, N)
        
        wald_l_month3_boot_Bt <- rep(NA, N)
        wald_r_month3_boot_Bt <- rep(NA, N)
        wald_l_month3_boot_Bt2 <- rep(NA, N)
        wald_r_month3_boot_Bt2 <- rep(NA, N)
        
        # for coverage
        Covg_wald_EY <- rep(NA, N)
        Covg_wald_Bt <- rep(NA, N)
        Covg_wald_boot_EY <- rep(NA, N)
        Covg_wald_boot_EY2 <- rep(NA, N)
        Covg_wald_boot_Bt <- rep(NA, N)
        Covg_wald_boot_Bt2 <- rep(NA, N)
        
        for(i in 1:N){
          
          Month_3[i] <- orig_treat[[i]][[1]][index_alpha, index_month]
          Var_3[i] <- orig_treat[[i]][[2]][index_alpha, index_month]
          Boot_temp_EY <- orig_treat[[i]][[3]][[index_alpha]][,  index_month] # on E[Y(t)] scale
          Boot_temp_Bt <- logit(orig_treat[[i]][[3]][[index_alpha]][,  index_month]) # on B(t) * beta scale
          Var_boot1_EY <- orig_treat[[i]][[4]][[index_alpha]][,  index_month] 
          Var_boot2_EY <- orig_treat[[i]][[5]][[index_alpha]][,  index_month] 
          
          ## E[Y(t)] scale
          # Wald confidence interval - on E[Y(t)] scale 
          temp1 <- wald(Month_3[i], Var_3[i])
          wald_l_month3_EY[i]  <- temp1[1]
          wald_r_month3_EY[i]  <- temp1[2]
          
          ### B(t)*beta scale
          # Wald confidence interval - trans
          temp3 <- wald(logit(Month_3[i]), Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2 )
          wald_l_month3_Bt[i]  <- expit(temp3[1])
          wald_r_month3_Bt[i]  <- expit(temp3[2])
          
          ## multipler bootstrap - EY
          
          ## quantile - less than others
          # quant_EY <- quantile(Boot_temp_EY - Month_3[i], probs = c(0.025, 0.975))
          # wald_l_month3_boot_EY[i] <- Month_3[i] - quant_EY[2]
          # wald_r_month3_boot_EY[i] <- Month_3[i] - quant_EY[1]
          # 
          # ## sd -- similar as others
          # temp4_EY <- wald(Month_3[i], var(Boot_temp_EY))
          # wald_l_month3_boot_EY[i]  <- temp4_EY[1]
          # wald_r_month3_boot_EY[i]  <- temp4_EY[2]
          
          ## boot-t - less than others 
          boot_t_EY <- (Month_3[i] - Boot_temp_EY) / sqrt(Var_boot1_EY)
          quant_EY <- quantile(boot_t_EY, probs = c(0.025, 0.975))
          wald_l_month3_boot_EY[i] <- Month_3[i] + quant_EY[1] * sqrt(Var_3[i])
          wald_r_month3_boot_EY[i] <- Month_3[i] + quant_EY[2] * sqrt(Var_3[i])
          
          boot_t_EY2 <- (Month_3[i] - Boot_temp_EY) / sqrt(Var_boot2_EY)
          quant_EY2 <- quantile(boot_t_EY2, probs = c(0.025, 0.975))
          wald_l_month3_boot_EY2[i] <- Month_3[i] + quant_EY2[1] * sqrt(Var_3[i])
          wald_r_month3_boot_EY2[i] <- Month_3[i] + quant_EY2[2] * sqrt(Var_3[i])
          
          
          # quant_Bt <- quantile(Boot_temp_Bt - logit(Month_3[i]), probs = c(0.025, 0.975))
          # wald_l_month3_boot_Bt[i] <- expit(log(Month_3[i]) - quant_Bt[2])
          # wald_r_month3_boot_Bt[i] <- expit(log(Month_3[i]) - quant_Bt[1])
          # 
          # temp4_Bt <- wald(logit(Month_3[i]), var(Boot_temp_Bt))
          # wald_l_month3_boot_Bt[i]  <- expit(temp4_Bt[1])
          # wald_r_month3_boot_Bt[i]  <- expit(temp4_Bt[2])
          
          boot_t_Bt <- (logit(Month_3[i]) - Boot_temp_Bt) / sqrt(Var_boot1_EY / (Boot_temp_EY * (1 - Boot_temp_EY))^2)
          quant_Bt <- quantile(boot_t_Bt, probs = c(0.025, 0.975), na.rm = T)
          wald_l_month3_boot_Bt[i] <- expit(logit(Month_3[i]) + quant_Bt[1] * sqrt(Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2))
          wald_r_month3_boot_Bt[i] <- expit(logit(Month_3[i]) + quant_Bt[2] * sqrt(Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2))
          
          boot_t_Bt2 <- (logit(Month_3[i]) - Boot_temp_Bt) / sqrt(Var_boot2_EY / (Boot_temp_EY * (1 - Boot_temp_EY))^2)
          quant_Bt2 <- quantile(boot_t_Bt2, probs = c(0.025, 0.975), na.rm = T)
          wald_l_month3_boot_Bt2[i] <- expit(logit(Month_3[i]) + quant_Bt2[1] * sqrt(Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2))
          wald_r_month3_boot_Bt2[i] <- expit(logit(Month_3[i]) + quant_Bt2[2] * sqrt(Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2))
          
          
          # boot_t_Bt2 <- (logit(Month_3[i]) - Boot_temp_Bt) / sqrt(Var_boot2_EY / (Boot_temp_EY * (1 - Boot_temp_EY))^2)
          # quant_Bt2 <- quantile(boot_t_Bt2, probs = c(0.025, 0.975))
          # wald_l_month3_boot_Bt2[i] <- expit(logit(Month_3[i]) + quant_Bt2[1] * sqrt(Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2))
          # wald_r_month3_boot_Bt2[i] <- expit(logit(Month_3[i]) + quant_Bt2[2] * sqrt(Var_3[i] / (Month_3[i] * (1 - Month_3[i]))^2))
          # 
          
          # Coverage
          Covg_wald_EY[i] <- ifelse((wald_l_month3_EY[i] < true_value) & (true_value < wald_r_month3_EY[i]), 1, 0)
          Covg_wald_Bt[i] <- ifelse((wald_l_month3_Bt[i] < true_value) & (true_value < wald_r_month3_Bt[i]), 1, 0)
          Covg_wald_boot_EY[i] <- ifelse((wald_l_month3_boot_EY[i] < true_value) & (true_value < wald_r_month3_boot_EY[i]), 1, 0)
          Covg_wald_boot_EY2[i] <- ifelse((wald_l_month3_boot_EY2[i] < true_value) & (true_value < wald_r_month3_boot_EY2[i]), 1, 0)
          Covg_wald_boot_Bt[i] <- ifelse((wald_l_month3_boot_Bt[i] < true_value) & (true_value < wald_r_month3_boot_Bt[i]), 1, 0)
          Covg_wald_boot_Bt2[i] <- ifelse((wald_l_month3_boot_Bt2[i] < true_value) & (true_value < wald_r_month3_boot_Bt2[i]), 1, 0)
          
        }
        
        result <- data.frame(Emp_mean = mean(Month_3), 
                             bias = mean(Month_3) - true_value,
                             Wald_EY = mean(Covg_wald_EY), 
                             Wald_Bt = mean(Covg_wald_Bt), 
                             Wald_Bt_boot_EY = mean(Covg_wald_boot_EY),
                             Wald_Bt_boot_EY2 = mean(Covg_wald_boot_EY2),
                             Wald_Bt_boot_Bt = mean(Covg_wald_boot_Bt),
                             Wald_Bt_boot_Bt2 = mean(Covg_wald_boot_Bt2)
        )
        result
        
        return(result)
      }
      
      f_res_point <- function(alpha = -0.6,
                              weight_type){
        
        result <- c()
        for(month_i in 1:4){
          
          month <- Month[month_i]
          
          temp <- get(paste0("month", month, "_treat"))
          index_alpha <- which(Alpha_seq == alpha)
          true_value <- temp[index_alpha]
          
          result1 <- f_point(alpha, month, true_value, weight_type)
          result <- rbind(result, c(alpha, true_value, round(result1, 3)))
          
        }
        colnames(result)[1:2] <- c("alpha", "true_value")
        result
        
      }
      
    }
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/binary_simu/simu_treat")
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_1e7_200_boot_rade")
    sample_size <- 200
    boot_size <- 200
    
    # no jackknife
    norm_200_neg06 <- f_res_point(alpha = -0.6, weight_type = "weight1")
    norm_200_neg03 <- f_res_point(alpha = -0.3, weight_type = "weight1")
    norm_200_00 <- f_res_point(alpha = 0, weight_type = "weight1")
    norm_200_03 <- f_res_point(alpha = 0.3, weight_type = "weight1")
    norm_200_06 <- f_res_point(alpha = 0.6, weight_type = "weight1")
    
    
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_1e7_400_boot_norm")
    sample_size <- 400
    boot_size <- 200
    
    # no jackknife
    norm_400_neg06 <- f_res_point(alpha = -0.6, weight_type = "weight1")
    norm_400_neg03 <- f_res_point(alpha = -0.3, weight_type = "weight1")
    norm_400_00 <- f_res_point(alpha = 0, weight_type = "weight1")
    norm_400_03 <- f_res_point(alpha = 0.3, weight_type = "weight1")
    norm_400_06 <- f_res_point(alpha = 0.6, weight_type = "weight1")
    
    
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_1e7_200_boot_rade")
    sample_size <- 200
    boot_size <- 200
    
    # no jackknife
    rade_200_neg06 <- f_res_point(alpha = -0.6, weight_type = "weight1")
    rade_200_neg03 <- f_res_point(alpha = -0.3, weight_type = "weight1")
    rade_200_00 <- f_res_point(alpha = 0, weight_type = "weight1")
    rade_200_03 <- f_res_point(alpha = 0.3, weight_type = "weight1")
    rade_200_06 <- f_res_point(alpha = 0.6, weight_type = "weight1")
    
    
    
    orig_treat <-  readRDS("ARC_simu_res_trt_dfsane_weight1_1e7_400_boot_rade")
    sample_size <- 400
    boot_size <- 200
    
    # no jackknife
    rade_400_neg06 <- f_res_point(alpha = -0.6, weight_type = "weight1")
    rade_400_neg03 <- f_res_point(alpha = -0.3, weight_type = "weight1")
    rade_400_00 <- f_res_point(alpha = 0, weight_type = "weight1")
    rade_400_03 <- f_res_point(alpha = 0.3, weight_type = "weight1")
    rade_400_06 <- f_res_point(alpha = 0.6, weight_type = "weight1")
    
    
    
  }
  
}