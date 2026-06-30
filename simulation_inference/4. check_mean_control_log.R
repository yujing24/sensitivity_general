rm(list=ls())
library(dplyr)


#### parametric bootstrap - single index model
{
  
  ### setting
  {
    
    # # month <- 3
    # month3_treat <- c(1.601, 1.876, 2.191, 2.563, 3.008)
    # month3_control <- c(1.551, 1.930, 2.342, 2.774, 3.211)
    # 
    # # month <- 6
    # month6_treat <- c(1.554, 1.821, 2.130, 2.496, 2.945)
    # month6_control <- c(1.521, 1.895, 2.304, 2.736, 3.181)
    # 
    # # month <- 9
    # month9_treat <- c(1.514, 1.777, 2.076, 2.441, 2.891)
    # month9_control <- c(1.495, 1.865, 2.269, 2.702, 3.150)
    # 
    # # month <- 12
    # month12_treat <- c(1.486, 1.743, 2.039, 2.402, 2.848)
    # month12_control <- c(1.478, 1.842, 2.245, 2.679, 3.127)
    
    # ## Below is the single index model + fixcoef
    # # month <- 3
    # month3_treat <- c(1.141, 1.434, 1.78, 2.191, 2.703, 3.341)
    # month3_control <- c(0.94, 1.323, 1.799, 2.342, 2.92, 3.497)
    # 
    # # month <- 6
    # month6_treat <- c(1.113, 1.394, 1.728, 2.13, 2.636, 3.283)
    # month6_control <- c(0.921, 1.297, 1.766, 2.304, 2.884, 3.472)
    # 
    # # month <- 9
    # month9_treat <- c(1.088, 1.359, 1.686, 2.076, 2.58, 3.23)
    # month9_control <- c(0.905, 1.275, 1.737, 2.269, 2.85, 3.443)
    # 
    # # month <- 12
    # month12_treat <- c(1.069, 1.334, 1.654, 2.039, 2.541, 3.188)
    # month12_control <- c(0.895, 1.26, 1.716, 2.245, 2.828, 3.422)
    # 
    # negative binomial + log link
    month3_treat <- c(1.16734, 1.39203, 1.69211, 2.09621, 2.61513, 3.22004)
    month6_treat <- c(1.13348, 1.34842, 1.63375, 2.01638, 2.50784, 3.09132)
    month9_treat <-  c(1.11376, 1.32203, 1.59900, 1.96872, 2.44638, 3.01779)
    month12_treat <- c(1.09788, 1.30106, 1.57155, 1.93054, 2.39640, 2.95888)
    
    month3_control <- c(1.02777, 1.27702, 1.63975, 2.16540, 2.86166, 3.63791)
    month6_control <- c(1.00005, 1.23851, 1.58236, 2.08165, 2.74729, 3.51121)
    month9_control <- c(0.98010, 1.21039, 1.54180, 2.02261, 2.67015, 3.42586)
    month12_control <- c(0.96400, 1.18720, 1.50757, 1.97351, 2.60517, 3.35527)
    
    
    
    Alpha_seq <- c(-1.2, -0.8, -0.4, 0, 0.4, 0.8)
    Month <- c(3, 6, 9, 12)
    N <- 500
    
    wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
      CI_left <- mean_orig - Zstar * sqrt(var_orig)
      CI_right <- mean_orig + Zstar * sqrt(var_orig)
      
      return(c(CI_left, CI_right))
    }
    
    f_point <- function(alpha, month, true_value,
                        weight_type = "weight1"){
      
      
      Alpha_seq <- c(-1.2, -0.8, -0.4, 0, 0.4, 0.8)
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
        if(weight_type == 'weight1'){
          temp0 <- try(temp_paraboot <- readRDS(paste0("para_boot_weight1/sim_paraboot_crl_log_200_", i)) , silent = T)
        }else{
          temp0 <- try(temp_paraboot <- readRDS(paste0("para_boot_weight2/sim_paraboot_crl_log_200_", i)) , silent = T)
        }
        
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
            Est_paraboot_Bt[j] <- log(Est_paraboot_EY[j])
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
          
          # normal CI using bootstrap se
          temp4_EY <- wald(Month_3[i], var(Est_paraboot_EY))
          wald_l_month3_boot_EY_se[i]  <- temp4_EY[1]
          wald_r_month3_boot_EY_se[i]  <- temp4_EY[2]
          
          temp4_Bt <- wald(log(Month_3[i]), var(log(Est_paraboot_EY)))
          wald_l_month3_boot_Bt_se[i]  <- exp(temp4_Bt[1])
          wald_r_month3_boot_Bt_se[i]  <- exp(temp4_Bt[2])
          
        }
        
        # Coverage
        Covg_wald_EY[i] <- ifelse((wald_l_month3_EY[i] < true_value) & (true_value < wald_r_month3_EY[i]), 1, 0)
        Covg_wald_Bt[i] <- ifelse((wald_l_month3_Bt[i] < true_value) & (true_value < wald_r_month3_Bt[i]), 1, 0)
        Covg_wald_boot_EY[i] <- ifelse((wald_l_month3_boot_EY[i] < true_value) & (true_value < wald_r_month3_boot_EY[i]), 1, 0)
        Covg_wald_boot_Bt[i] <- ifelse((wald_l_month3_boot_Bt[i] < true_value) & (true_value < wald_r_month3_boot_Bt[i]), 1, 0)
        Covg_wald_boot_EY_se[i] <- ifelse((wald_l_month3_boot_EY_se[i] < true_value) & (true_value < wald_r_month3_boot_EY_se[i]), 1, 0)
        Covg_wald_boot_Bt_se[i] <- ifelse((wald_l_month3_boot_Bt_se[i] < true_value) & (true_value < wald_r_month3_boot_Bt_se[i]), 1, 0)
        
      }
      
      result <- data.frame(Emp_mean = mean(Month_3), 
                           bias = mean(Month_3) - true_value,
                           Wald_EY = mean(Covg_wald_EY, na.rm = T), 
                           Wald_Bt = mean(Covg_wald_Bt, na.rm = T), 
                           Wald_Bt_paraboot_EY = mean(Covg_wald_boot_EY, na.rm = T),
                           Wald_Bt_paraboot_Bt = mean(Covg_wald_boot_Bt, na.rm = T),
                           Wald_Bt_paraboot_EY_se = mean(Covg_wald_boot_EY_se, na.rm = T),
                           Wald_Bt_paraboot_Bt_se = mean(Covg_wald_boot_Bt_se, na.rm = T)
      )
      
      result
      
      
      return(result)
    }
    
    f_res_point <- function(alpha = -0.6,
                            weight_type){
      
      result <- c()
      for(month_i in 1:4){
        
        month <- Month[month_i]
        
        temp <- get(paste0("month", month, "_control"))
        index_alpha <- which(Alpha_seq == alpha)
        true_value <- temp[index_alpha]
        
        result1 <- f_point(alpha, month, true_value, weight_type)
        result <- rbind(result, c(alpha, true_value, round(result1, 3)))
        
      }
      colnames(result)[1:2] <- c("alpha", "true_value")
      result
      
    }
    
  }
  
  
  # negative + norm1
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/neg_simu_log_norm1/simu_control")
    sample_size <- 200
    boot_size <- 200
    N <- 500
    
    # weight 1
    orig_treat <-  readRDS("ARC_simu_res_crl_dfsane_weight1_1e7_200_GL_30_norm1")
    f_res_point(alpha = -1.2, weight_type = "weight1")
    f_res_point(alpha = -0.8, weight_type = "weight1")
    f_res_point(alpha = -0.4, weight_type = "weight1")
    f_res_point(alpha = 0, weight_type = "weight1")
    f_res_point(alpha = 0.4, weight_type = "weight1")
    f_res_point(alpha = 0.8, weight_type = "weight1")
    
    
    # weight 2
    orig_treat <-  readRDS("ARC_simu_res_crl_dfsane_weight2_1e7_200_GL_30")
    f_res_point(alpha = -0.6, weight_type = "weight2")
    f_res_point(alpha = -0.3, weight_type = "weight2")
    f_res_point(alpha = 0, weight_type = "weight2")
    f_res_point(alpha = 0.3, weight_type = "weight2")
    f_res_point(alpha = 0.6, weight_type = "weight2")
    
  }
  
}


### old ### 
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
    
    N_size <- 200
    
    wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
      CI_left <- mean_orig - Zstar * sqrt(var_orig)
      CI_right <- mean_orig + Zstar * sqrt(var_orig)
      
      return(c(CI_left, CI_right))
    }
    
    f <- function(alpha, month, true_value, 
                  link_type = "identity", 
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
      
      # for coverage
      Covg_wald <- rep(NA, N)
      Covg_wald_trans <- rep(NA, N)
      Covg_JN <- rep(NA, N)
      Covg_JN_trans <- rep(NA, N)
      
      for(i in 1:N){
        
        Month_3[i] <- orig_treat[[i]][[1]][index_alpha, index_month]
        Var_3[i] <- orig_treat[[i]][[2]][index_alpha, index_month]
        
        # Wald confidence interval 
        temp2 <- wald(Month_3[i], Var_3[i])
        wald_l_month3[i]  <- temp2[1]
        wald_r_month3[i]  <- temp2[2]
        
        # Wald confidence interval - trans
        temp3 <- wald(log(Month_3[i]), Var_3[i] / (Month_3[i])^2 )
        wald_l_month3_trans[i]  <- exp(temp3[1])
        wald_r_month3_trans[i]  <- exp(temp3[2])
        
        # Coverage
        Covg_wald[i] <- ifelse((wald_l_month3[i] < true_value) & (true_value < wald_r_month3[i]), 1, 0)
        Covg_wald_trans[i] <- ifelse((wald_l_month3_trans[i] < true_value) & (true_value < wald_r_month3_trans[i]), 1, 0)
        
        if(i %in% data_list_JN){
          
          # load the JN variance
          if(link_type == "identity"){
            temp_JN <- readRDS(paste0("JN_simu_control_1e7/simu_JN_", i, "_crl"))
          }else{
            if(weight_type == "weight1"){
              temp_JN <- readRDS(paste0("JN_simu_control_weight1_1e7/simu_JN_", i, "_crl")) 
            }else{
              temp_JN <- readRDS(paste0("JN_simu_control_weight2_1e7/simu_JN_", i, "_crl")) 
            }
          }
          
          Est_JN <- rep(NA, N_size)
          Est_JN_trans <- rep(NA, N_size)
          for(j in 1:N_size){
            temp <- try(temp_JN[[j]][[1]], silent = T)
            if(length(temp) > 1){
              Est_JN[j] <- temp_JN[[j]][[1]][index_alpha, index_month]
              Est_JN_trans[j] <- log(temp_JN[[j]][[1]][index_alpha, index_month])
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
            temp5 <- wald(log(Month_3[i]), var_JN_trans)
            JN_l_month3_trans[i] <- exp(temp5[1])
            JN_r_month3_trans[i] <- exp(temp5[2])
            
            # Coverage
            Covg_JN[i] <-  ifelse((JN_l_month3[i] < true_value) & (true_value < JN_r_month3[i]), 1, 0)
            Covg_JN_trans[i] <-  ifelse((JN_l_month3_trans[i] < true_value) & (true_value < JN_r_month3_trans[i]), 1, 0)
            
          }
          
        }
        
      }
      
      result <- data.frame(Emp_mean = mean(Month_3), 
                           bias = mean(Month_3) - true_value,
                           Wald = mean(Covg_wald),
                           Wald_trans = mean(Covg_wald_trans),
                           JN = mean(Covg_JN, na.rm = T), 
                           JN_trans = mean(Covg_JN_trans, na.rm = T)
      )
      return(result)
    }
    
    f_res <- function(alpha = -0.6,
                      link_type,
                      weight_type){
      
      result <- c()
      for(month_i in 1:4){
        
        month <- Month[month_i]
        
        temp <- get(paste0("month", month, "_control"))
        index_alpha <- which(Alpha_seq == alpha)
        true_value <- temp[index_alpha]
        
        result1 <- f(alpha, month, true_value, 
                     link_type, weight_type)
        result <- rbind(result, c(true_value, round(result1, 3)))
        
      }
      result
      
    }
    
  }
  
  
  ### log link - weight1
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/count_simu/simu_control")
    orig_treat <-  readRDS("ARC_simu_res_crl_dfsane_weight1_1e7")
    orig_treat <-  readRDS("ARC_simu_res_crl_dfsane_weight1_1e7_400")
    
    data_list <- c(1:500)
    data_list_JN <- c() 
    for(k in 1:500){
      temp_data <- try(readRDS(paste0("JN_simu_control_weight1_1e7/simu_JN_", k, "_crl")) , silent = TRUE)
      if(length(temp_data) > 1){
        data_list_JN <- c(data_list_JN, k)
      }
    }
    length(data_list_JN)
    setdiff(data_list, data_list_JN)
    
    
    f_res(alpha = -0.6, link_type = "log", weight_type = "weight1") 
    f_res(alpha = -0.3, link_type = "log", weight_type = "weight1") 
    f_res(alpha = 0, link_type = "log", weight_type = "weight1") 
    f_res(alpha = 0.3, link_type = "log", weight_type = "weight1") 
    f_res(alpha = 0.6, link_type = "log", weight_type = "weight1") 
    
    
  }
  
  
  ### log link - weight2
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/count_simu/simu_control")
    orig_treat <-  readRDS("ARC_simu_res_crl_dfsane_weight2_1e7")
    
    data_list <- c(1:500)
    data_list_JN <- c() 
    for(k in 1:500){
      temp_data <- try(readRDS(paste0("JN_simu_control_weight2_1e7/simu_JN_", k, "_crl")) , silent = TRUE)
      if(length(temp_data) > 1){
        data_list_JN <- c(data_list_JN, k)
      }
    }
    length(data_list_JN)
    setdiff(data_list, data_list_JN)
    
    
    f_res(alpha = -0.6, link_type = "log", weight_type = "weight2")
    f_res(alpha = -0.3, link_type = "log", weight_type = "weight2")
    f_res(alpha = 0, link_type = "log", weight_type = "weight2")
    f_res(alpha = 0.3, link_type = "log", weight_type = "weight2")
    f_res(alpha = 0.6, link_type = "log", weight_type = "weight2")
    
  }
  
  
  ### identity link
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/identity_simu/simu_control")
    orig_treat <-  readRDS("ARC_simu_res_crl_1e7")
    
    data_list <- c(1:500)
    data_list_JN <- c() 
    for(k in 1:500){
      temp_data <- try( readRDS(paste0("JN_simu_control_1e7/simu_JN_", k, "_crl")) , silent = TRUE)
      if(length(temp_data) > 1){
        data_list_JN <- c(data_list_JN, k)
      }
    }
    length(data_list_JN)
    setdiff(data_list, data_list_JN)
    
    
    f_res(alpha = -0.6, link_type = "identity")
  }
  
  
  #### Edgeworth expansion for skewness-corrected CI
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
          
          
          ### B(t)*beta scale
          ## B(t)*beta scale (transfered)- skewness-corrected
          # Wald confidence interval - trans
          temp3 <- wald(log(Month_3[i]), Var_3[i] / (Month_3[i])^2)
          wald_l_month3_Bt[i]  <- exp(temp3[1])
          wald_r_month3_Bt[i]  <- exp(temp3[2])
          
          # Wald confidence interval - on B(t)*beta scale - skewness-corrected
          temp4 <- wald_skew(log(Month_3[i]), Var_3[i] / (Month_3[i])^2, Zstar = qnorm(0.975), skew = Skew_Bt[i])
          wald_l_month3_Bt_skew[i]  <- exp(temp4[1])
          wald_r_month3_Bt_skew[i]  <- exp(temp4[2])
          
          # Wald confidence interval - on B(t)*beta scale - skewness-corrected (peter)
          temp4_peter <- wald_skew_peter(log(Month_3[i]), Var_3[i] / (Month_3[i])^2, Zstar = qnorm(0.975), skew = Skew_Bt[i])
          wald_l_month3_Bt_skew_peter[i]  <- exp(temp4_peter[1])
          wald_r_month3_Bt_skew_peter[i]  <- exp(temp4_peter[2])
          
          
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
          
          temp <- get(paste0("month", month, "_control"))
          index_alpha <- which(Alpha_seq == alpha)
          true_value <- temp[index_alpha]
          
          result1 <- f_point(alpha, month, true_value, weight_type)
          result <- rbind(result, c(alpha, true_value, round(result1, 3)))
          
        }
        colnames(result)[1:2] <- c("alpha", "true_value")
        result
        
      }
      
    }
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/count_simu/simu_control")
    orig_treat <-  readRDS("ARC_simu_res_crl_dfsane_weight1_1e7_400_skew")
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
          Boot_temp_EY <- orig_treat[[i]][[3]][[index_alpha]][,  index_month] # on B(t) * beta scale
          Boot_temp_Bt <- log(orig_treat[[i]][[3]][[index_alpha]][,  index_month]) # on B(t) * beta scale
          Var_boot1_EY <- orig_treat[[i]][[4]][[index_alpha]][,  index_month] 
          Var_boot2_EY <- orig_treat[[i]][[5]][[index_alpha]][,  index_month] 
          
          ## E[Y(t)] scale
          # Wald confidence interval - on E[Y(t)] scale 
          temp1 <- wald(Month_3[i], Var_3[i])
          wald_l_month3_EY[i]  <- temp1[1]
          wald_r_month3_EY[i]  <- temp1[2]
          
          ### B(t)*beta scale
          # Wald confidence interval - trans
          temp3 <- wald(log(Month_3[i]), Var_3[i] / (Month_3[i])^2)
          wald_l_month3_Bt[i]  <- exp(temp3[1])
          wald_r_month3_Bt[i]  <- exp(temp3[2])
          
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
          
          
          # quant_Bt <- quantile(Boot_temp_Bt - log(Month_3[i]), probs = c(0.025, 0.975))
          # wald_l_month3_boot_Bt[i] <- exp(log(Month_3[i]) - quant_Bt[2])
          # wald_r_month3_boot_Bt[i] <- exp(log(Month_3[i]) - quant_Bt[1])
          
          # temp4_Bt <- wald(log(Month_3[i]), var(Boot_temp_Bt))
          # wald_l_month3_boot_Bt[i]  <- exp(temp4_Bt[1])
          # wald_r_month3_boot_Bt[i]  <- exp(temp4_Bt[2])
          
          boot_t_Bt <- (log(Month_3[i]) - Boot_temp_Bt) / sqrt(Var_boot1_EY / (Boot_temp_EY)^2)
          quant_Bt <- quantile(boot_t_Bt, probs = c(0.025, 0.975))
          wald_l_month3_boot_Bt[i] <- exp(log(Month_3[i]) + quant_Bt[1] * sqrt(Var_3[i] / (Month_3[i])^2))
          wald_r_month3_boot_Bt[i] <- exp(log(Month_3[i]) + quant_Bt[2] * sqrt(Var_3[i] / (Month_3[i])^2))
          
          boot_t_Bt2 <- (log(Month_3[i]) - Boot_temp_Bt) / sqrt(Var_boot2_EY / (Boot_temp_EY)^2)
          quant_Bt2 <- quantile(boot_t_Bt2, probs = c(0.025, 0.975))
          wald_l_month3_boot_Bt2[i] <- exp(log(Month_3[i]) + quant_Bt2[1] * sqrt(Var_3[i] / (Month_3[i])^2))
          wald_r_month3_boot_Bt2[i] <- exp(log(Month_3[i]) + quant_Bt2[2] * sqrt(Var_3[i] / (Month_3[i])^2))
          
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
          
          temp <- get(paste0("month", month, "_control"))
          index_alpha <- which(Alpha_seq == alpha)
          true_value <- temp[index_alpha]
          
          result1 <- f_point(alpha, month, true_value, weight_type)
          result <- rbind(result, c(alpha, true_value, round(result1, 3)))
          
        }
        colnames(result)[1:2] <- c("alpha", "true_value")
        result
        
      }
      
    }
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/count_simu/simu_control")
    orig_treat <-  readRDS("ARC_simu_res_crl_dfsane_weight1_1e7_200_boot_norm")
    sample_size <- 200
    boot_size <- 200
    N <- 500
    
    # no jackknife
    norm_200_neg06 <- f_res_point(alpha = -0.6, weight_type = "weight1")
    norm_200_neg03 <- f_res_point(alpha = -0.3, weight_type = "weight1")
    norm_200_00 <- f_res_point(alpha = 0, weight_type = "weight1")
    norm_200_03 <- f_res_point(alpha = 0.3, weight_type = "weight1")
    norm_200_06 <- f_res_point(alpha = 0.6, weight_type = "weight1")
    
    
    orig_treat <-  readRDS("ARC_simu_res_crl_dfsane_weight1_1e7_400_boot_norm")
    sample_size <- 400
    boot_size <- 200
    
    # no jackknife
    norm_400_neg06 <- f_res_point(alpha = -0.6, weight_type = "weight1")
    norm_400_neg03 <- f_res_point(alpha = -0.3, weight_type = "weight1")
    norm_400_00 <- f_res_point(alpha = 0, weight_type = "weight1")
    norm_400_03 <- f_res_point(alpha = 0.3, weight_type = "weight1")
    norm_400_06 <- f_res_point(alpha = 0.6, weight_type = "weight1")
    
    
    orig_treat <-  readRDS("ARC_simu_res_crl_dfsane_weight1_1e7_200_boot_rade")
    sample_size <- 200
    boot_size <- 200
    
    # no jackknife
    rade_200_neg06 <- f_res_point(alpha = -0.6, weight_type = "weight1")
    rade_200_neg03 <- f_res_point(alpha = -0.3, weight_type = "weight1")
    rade_200_00 <- f_res_point(alpha = 0, weight_type = "weight1")
    rade_200_03 <- f_res_point(alpha = 0.3, weight_type = "weight1")
    rade_200_06 <- f_res_point(alpha = 0.6, weight_type = "weight1")
    
    
    
    orig_treat <-  readRDS("ARC_simu_res_crl_dfsane_weight1_1e7_400_boot_rade")
    sample_size <- 400
    boot_size <- 200
    
    # no jackknife
    rade_400_neg06 <- f_res_point(alpha = -0.6, weight_type = "weight1")
    rade_400_neg03 <- f_res_point(alpha = -0.3, weight_type = "weight1")
    rade_400_00 <- f_res_point(alpha = 0, weight_type = "weight1")
    rade_400_03 <- f_res_point(alpha = 0.3, weight_type = "weight1")
    rade_400_06 <- f_res_point(alpha = 0.6, weight_type = "weight1")
    
  }
  
  
  #### parametric bootstrap - negative binomial model
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
          if(weight_type == 'weight1'){
            temp0 <- try(temp_paraboot <- readRDS(paste0("para_boot_weight1/sim_paraboot_crl_count_200_", i)) , silent = T)
          }else{
            temp0 <- try(temp_paraboot <- readRDS(paste0("para_boot_weight2/sim_paraboot_crl_count_200_", i)) , silent = T)
          }
          
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
      
      f_res_point <- function(alpha = -0.6,
                              weight_type){
        
        result <- c()
        for(month_i in 1:4){
          
          month <- Month[month_i]
          
          temp <- get(paste0("month", month, "_control"))
          index_alpha <- which(Alpha_seq == alpha)
          true_value <- temp[index_alpha]
          
          result1 <- f_point(alpha, month, true_value, weight_type)
          result <- rbind(result, c(alpha, true_value, round(result1, 3)))
          
        }
        colnames(result)[1:2] <- c("alpha", "true_value")
        result
        
      }
      
    }
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/single_simu_log/simu_control")
    sample_size <- 200
    boot_size <- 200
    N <- 500
    
    # weight 1
    orig_treat <-  readRDS("ARC_simu_res_crl_dfsane_weight1_1e7_200_GL_30")
    f_res_point(alpha = -0.6, weight_type = "weight1")
    f_res_point(alpha = -0.3, weight_type = "weight1")
    f_res_point(alpha = 0, weight_type = "weight1")
    f_res_point(alpha = 0.3, weight_type = "weight1")
    f_res_point(alpha = 0.6, weight_type = "weight1")
    
    
    # weight 2
    orig_treat <-  readRDS("ARC_simu_res_crl_dfsane_weight2_1e7_200_GL_30")
    f_res_point(alpha = -0.6, weight_type = "weight2")
    f_res_point(alpha = -0.3, weight_type = "weight2")
    f_res_point(alpha = 0, weight_type = "weight2")
    f_res_point(alpha = 0.3, weight_type = "weight2")
    f_res_point(alpha = 0.6, weight_type = "weight2")
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/count_simu/simu_control")
    miss_index <- c()
    for(i in 1:N){
      temp0 <- try(temp_paraboot <- readRDS(paste0("para_boot_weight1/sim_paraboot_crl_count_200_", i)) , silent = T)
      if(length(temp0) <= 1){
        miss_index <- c(miss_index, i)
      }
    }
    miss_index
    
    
  }
  
  
}


