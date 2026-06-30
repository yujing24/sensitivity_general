rm(list=ls())
library(dplyr)
# library(SensIAT)
library(orthogonalsplinebasis)
library(MASS)
library(ggplot2)
library(monotone) # unimonotone 



##### Point estimation  - negative binomial model simulated data
{
  
  
  # below is the true value
  
  alpha_crl_l_true <- -0.7492188
  alpha_crl_u_true <- 0.4492188
  alpha_trt_l_true <- -0.9585937
  alpha_trt_u_true <- 0.6254688
  
  
  #### One-step multiplier bootstrap - log link
  {
    
    ## Treat group
    {
      
      link_type = "log"
      sample_size = 200
      group = "treat"
      IF_type = "weight2"
      
      if(group == "treat"){
        grp = "trt"
      }else{
        grp = "crl"
      }
      Boot_num <- 100
      
      
      simu_size <- 500
      
      Alpha_u_trt_point_orig <- rep(NA, simu_size)
      Alpha_l_trt_point_orig <- rep(NA, simu_size) 
      Alpha_u_trt_point_corr_adap <- rep(NA, simu_size)
      Alpha_l_trt_point_corr_adap <- rep(NA, simu_size) 
      
      Alpha_u_trt_point_corr_NA_count_adap <- rep(NA, simu_size)
      Alpha_l_trt_point_corr_NA_count_adap <- rep(NA, simu_size)
      
      Convergence_all <- rep(NA, simu_size)
      
      setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
      
      for(i in 1:simu_size){
         
        print(i)
        data_orig <- try(readRDS(file = paste0('orig_point_est/', 
                                               link_type, "_", sample_size,
                                               "/",  
                                               group, "_", IF_type, 
                                               "/select_", grp, "_", sample_size, "_", link_type, "_", IF_type, "_", i)), silent = T)
        data_bias_adap <-  try(readRDS(file = paste0("bias_correct_", Boot_num,
                                                     "/",
                                                     link_type, "_", sample_size,
                                                     "/", 
                                                     group, "_", IF_type, 
                                                     "/", 
                                                     "adaptive_", grp, "_", sample_size, "_", link_type, "_", IF_type,"_", i, "_corr")), silent = T)
        
        
        Alpha_l_trt_point_orig[i] <- data_orig[[1]][1]
        Alpha_u_trt_point_orig[i] <- data_orig[[1]][2]
        
        Boot_est_u_adap <- rep(NA, Boot_num)
        Boot_est_l_adap <- rep(NA, Boot_num)
        
        # for onestep
        for(j in 1:Boot_num){
          
          temp <- data_bias_adap[[j]]
          Boot_est_l_adap[j] <- temp[[1]][1]
          Boot_est_u_adap[j] <- temp[[1]][2]
          
        }
        
        Alpha_u_trt_point_corr_adap[i] <- 2 * Alpha_u_trt_point_orig[i] - mean(Boot_est_u_adap, na.rm = T)
        Alpha_l_trt_point_corr_adap[i] <- 2 * Alpha_l_trt_point_orig[i] - mean(Boot_est_l_adap, na.rm = T)
        Alpha_u_trt_point_corr_NA_count_adap[i] <- length(which(is.na(Boot_est_u_adap) == T))
        Alpha_l_trt_point_corr_NA_count_adap[i] <- length(which(is.na(Boot_est_l_adap) == T))
        
        Convergence <- data_orig[[5]]
        Convergence_all[i] <- length(which(Convergence > 0))
        
      }
      
      max(Convergence_all)
      max(Alpha_u_trt_point_corr_NA_count_adap) 
      max(Alpha_l_trt_point_corr_NA_count_adap) 
      Alpha_u_trt_point_corr_NA_count_adap[which(Alpha_u_trt_point_corr_NA_count_adap > 0)]
      Alpha_l_trt_point_corr_NA_count_adap[which(Alpha_l_trt_point_corr_NA_count_adap > 0)]
      
      saveRDS(Alpha_u_trt_point_orig, paste0('value/orig/Alpha_u_', grp, '_orig_log_', IF_type, '_', sample_size))
      saveRDS(Alpha_l_trt_point_orig, paste0('value/orig/Alpha_l_', grp, '_orig_log_', IF_type, '_', sample_size))
      
      saveRDS(Alpha_u_trt_point_corr_adap, paste0('value/bias_corr/Alpha_u_', grp, '_corr_log_', IF_type, '_', sample_size))
      saveRDS(Alpha_l_trt_point_corr_adap, paste0('value/bias_corr/Alpha_l_', grp, '_corr_log_', IF_type, '_', sample_size))
      
      mean(Alpha_u_trt_point_orig, na.rm = T) -  alpha_trt_u_true
      mean(Alpha_u_trt_point_corr_adap, na.rm = T) -  alpha_trt_u_true
      (mean(Alpha_u_trt_point_orig, na.rm = T) -  alpha_trt_u_true) /  alpha_trt_u_true
      (mean(Alpha_u_trt_point_corr_adap, na.rm = T) -  alpha_trt_u_true) /  alpha_trt_u_true
      
      mean(Alpha_l_trt_point_orig, na.rm = T) -  alpha_trt_l_true
      mean(Alpha_l_trt_point_corr_adap, na.rm = T) -  alpha_trt_l_true
      (mean(Alpha_l_trt_point_orig, na.rm = T) -  alpha_trt_l_true) /  alpha_trt_l_true
      (mean(Alpha_l_trt_point_corr_adap, na.rm = T) -  alpha_trt_l_true) /  alpha_trt_l_true
      
      
    }
    
    
    ## Control group
    {
      
      link_type = "log"
      sample_size = 200
      group = "control"
      IF_type = "weight2"
      
      if(group == "treat"){
        grp = "trt"
      }else{
        grp = "crl"
      }
      Boot_num <- 100
      
      
      simu_size <- 500
      
      Alpha_u_trt_point_orig <- rep(NA, simu_size)
      Alpha_l_trt_point_orig <- rep(NA, simu_size) 
      Alpha_u_trt_point_corr_adap <- rep(NA, simu_size)
      Alpha_l_trt_point_corr_adap <- rep(NA, simu_size) 
      
      Alpha_u_trt_point_corr_NA_count_adap <- rep(NA, simu_size)
      Alpha_l_trt_point_corr_NA_count_adap <- rep(NA, simu_size)
      
      Convergence_all <- rep(NA, simu_size)
      
      setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
      
      for(i in 1:simu_size){
        
        print(i)
        data_orig <- try(readRDS(file = paste0('orig_point_est/', 
                                               link_type, "_", sample_size,
                                               "/",  
                                               group, "_", IF_type, 
                                               "/select_", grp, "_", sample_size, "_", link_type, "_", IF_type, "_", i)), silent = T)
        data_bias_adap <-  try(readRDS(file = paste0("bias_correct_", Boot_num,
                                                     "/",
                                                     link_type, "_", sample_size,
                                                     "/", 
                                                     group, "_", IF_type, 
                                                     "/", 
                                                     "adaptive_", grp, "_", sample_size, "_", link_type, "_", IF_type,"_", i, "_corr")), silent = T)
        
        
        Alpha_l_trt_point_orig[i] <- data_orig[[1]][1]
        Alpha_u_trt_point_orig[i] <- data_orig[[1]][2]
        
        Boot_est_u_adap <- rep(NA, Boot_num)
        Boot_est_l_adap <- rep(NA, Boot_num)
        
        # for onestep
        for(j in 1:Boot_num){
          
          temp <- data_bias_adap[[j]]
          Boot_est_l_adap[j] <- temp[[1]][1]
          Boot_est_u_adap[j] <- temp[[1]][2]
          
        }
        
        Alpha_u_trt_point_corr_adap[i] <- 2 * Alpha_u_trt_point_orig[i] - mean(Boot_est_u_adap, na.rm = T)
        Alpha_l_trt_point_corr_adap[i] <- 2 * Alpha_l_trt_point_orig[i] - mean(Boot_est_l_adap, na.rm = T)
        Alpha_u_trt_point_corr_NA_count_adap[i] <- length(which(is.na(Boot_est_u_adap) == T))
        Alpha_l_trt_point_corr_NA_count_adap[i] <- length(which(is.na(Boot_est_l_adap) == T))
        
        Convergence <- data_orig[[5]]
        Convergence_all[i] <- length(which(Convergence > 0))
        
      }
      
      
      max(Convergence_all)
      max(Alpha_u_trt_point_corr_NA_count_adap) 
      max(Alpha_l_trt_point_corr_NA_count_adap) 
      Alpha_u_trt_point_corr_NA_count_adap[which(Alpha_u_trt_point_corr_NA_count_adap > 0)]
      Alpha_l_trt_point_corr_NA_count_adap[which(Alpha_l_trt_point_corr_NA_count_adap > 0)]
      
      saveRDS(Alpha_u_trt_point_orig, paste0('value/orig/Alpha_u_', grp, '_orig_log_', IF_type, '_', sample_size))
      saveRDS(Alpha_l_trt_point_orig, paste0('value/orig/Alpha_l_', grp, '_orig_log_', IF_type, '_', sample_size))
      
      saveRDS(Alpha_u_trt_point_corr_adap, paste0('value/bias_corr/Alpha_u_', grp, '_corr_log_', IF_type, '_', sample_size))
      saveRDS(Alpha_l_trt_point_corr_adap, paste0('value/bias_corr/Alpha_l_', grp, '_corr_log_', IF_type, '_', sample_size))
      
      mean(Alpha_u_trt_point_orig, na.rm = T) -  alpha_crl_u_true
      mean(Alpha_u_trt_point_corr_adap, na.rm = T) -  alpha_crl_u_true
      (mean(Alpha_u_trt_point_orig, na.rm = T) -  alpha_crl_u_true) /  alpha_crl_u_true
      (mean(Alpha_u_trt_point_corr_adap, na.rm = T) -  alpha_crl_u_true) /  alpha_crl_u_true
      
      mean(Alpha_l_trt_point_orig, na.rm = T) -  alpha_crl_l_true
      mean(Alpha_l_trt_point_corr_adap, na.rm = T) -  alpha_crl_l_true
      (mean(Alpha_l_trt_point_orig, na.rm = T) -  alpha_crl_l_true) /  alpha_crl_l_true
      (mean(Alpha_l_trt_point_corr_adap, na.rm = T) -  alpha_crl_l_true) /  alpha_crl_l_true
      
      
    }
    
    
  }

  
  #### One-step multiplier bootstrap - identity link
  {
    
    
    ## Treat group
    {
      
      link_type = "identity"
      sample_size = 200
      group = "treat"
      
      if(group == "treat"){
        grp = "trt"
      }else{
        grp = "crl"
      }
      Boot_num <- 100
      
      
      simu_size <- 500
      
      Alpha_u_trt_point_orig <- rep(NA, simu_size)
      Alpha_l_trt_point_orig <- rep(NA, simu_size) 
      Alpha_u_trt_point_corr_adap <- rep(NA, simu_size)
      Alpha_l_trt_point_corr_adap <- rep(NA, simu_size) 
      
      Alpha_u_trt_point_corr_NA_count_adap <- rep(NA, simu_size)
      Alpha_l_trt_point_corr_NA_count_adap <- rep(NA, simu_size)
      
      Convergence_all <- rep(NA, simu_size)
      
      setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
      
      for(i in 1:simu_size){
        
        print(i)
        data_orig <- try(readRDS(file = paste0('orig_point_est/', 
                                               link_type, "_", sample_size,
                                               "/",  
                                               group, 
                                               "/select_", grp, "_", sample_size, "_", link_type, "_", i)), silent = T)
        data_bias_adap <-  try(readRDS(file = paste0("bias_correct_", Boot_num,
                                                     "/",
                                                     link_type, "_", sample_size,
                                                     "/", 
                                                     group, 
                                                     "/", 
                                                     "adaptive_", grp, "_", sample_size, "_", link_type, "_", i, "_corr")), silent = T)
        
        
        Alpha_l_trt_point_orig[i] <- data_orig[[1]][1]
        Alpha_u_trt_point_orig[i] <- data_orig[[1]][2]
        
        Boot_est_u_adap <- rep(NA, Boot_num)
        Boot_est_l_adap <- rep(NA, Boot_num)
        
        # for onestep
        for(j in 1:Boot_num){
          
          temp <- data_bias_adap[[j]]
          Boot_est_l_adap[j] <- temp[[1]][1]
          Boot_est_u_adap[j] <- temp[[1]][2]
          
        }
        
        Alpha_u_trt_point_corr_adap[i] <- 2 * Alpha_u_trt_point_orig[i] - mean(Boot_est_u_adap, na.rm = T)
        Alpha_l_trt_point_corr_adap[i] <- 2 * Alpha_l_trt_point_orig[i] - mean(Boot_est_l_adap, na.rm = T)
        Alpha_u_trt_point_corr_NA_count_adap[i] <- length(which(is.na(Boot_est_u_adap) == T))
        Alpha_l_trt_point_corr_NA_count_adap[i] <- length(which(is.na(Boot_est_l_adap) == T))
        
        Convergence <- data_orig[[5]]
        Convergence_all[i] <- length(which(Convergence > 0))
        
      }
      
      max(Convergence_all)
      max(Alpha_u_trt_point_corr_NA_count_adap) 
      max(Alpha_l_trt_point_corr_NA_count_adap) 
      Alpha_u_trt_point_corr_NA_count_adap[which(Alpha_u_trt_point_corr_NA_count_adap > 0)]
      Alpha_l_trt_point_corr_NA_count_adap[which(Alpha_l_trt_point_corr_NA_count_adap > 0)]
      
      saveRDS(Alpha_u_trt_point_orig, paste0('value/orig/Alpha_u_', grp, '_orig_identity_', sample_size))
      saveRDS(Alpha_l_trt_point_orig, paste0('value/orig/Alpha_l_', grp, '_orig_identity_', sample_size))
      
      saveRDS(Alpha_u_trt_point_corr_adap, paste0('value/bias_corr/Alpha_u_', grp, '_corr_identity_', sample_size))
      saveRDS(Alpha_l_trt_point_corr_adap, paste0('value/bias_corr/Alpha_l_', grp, '_corr_identity_', sample_size))
      
      mean(Alpha_u_trt_point_orig, na.rm = T) -  alpha_trt_u_true
      mean(Alpha_u_trt_point_corr_adap, na.rm = T) -  alpha_trt_u_true
      (mean(Alpha_u_trt_point_orig, na.rm = T) -  alpha_trt_u_true) /  alpha_trt_u_true
      (mean(Alpha_u_trt_point_corr_adap, na.rm = T) -  alpha_trt_u_true) /  alpha_trt_u_true
      
      mean(Alpha_l_trt_point_orig, na.rm = T) -  alpha_trt_l_true
      mean(Alpha_l_trt_point_corr_adap, na.rm = T) -  alpha_trt_l_true
      (mean(Alpha_l_trt_point_orig, na.rm = T) -  alpha_trt_l_true) /  alpha_trt_l_true
      (mean(Alpha_l_trt_point_corr_adap, na.rm = T) -  alpha_trt_l_true) /  alpha_trt_l_true
      
      
    }
    
    
    ## Control group
    {
      
      link_type = "identity"
      sample_size = 200
      group = "control"
      
      if(group == "treat"){
        grp = "trt"
      }else{
        grp = "crl"
      }
      Boot_num <- 100
      
      
      simu_size <- 500
      
      Alpha_u_trt_point_orig <- rep(NA, simu_size)
      Alpha_l_trt_point_orig <- rep(NA, simu_size) 
      Alpha_u_trt_point_corr_adap <- rep(NA, simu_size)
      Alpha_l_trt_point_corr_adap <- rep(NA, simu_size) 
      
      Alpha_u_trt_point_corr_NA_count_adap <- rep(NA, simu_size)
      Alpha_l_trt_point_corr_NA_count_adap <- rep(NA, simu_size)
      
      Convergence_all <- rep(NA, simu_size)
      
      setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
      
      for(i in 1:simu_size){
        
        print(i)
        data_orig <- try(readRDS(file = paste0('orig_point_est/', 
                                               link_type, "_", sample_size,
                                               "/",  
                                               group, 
                                               "/select_", grp, "_", sample_size, "_", link_type,  "_", i)), silent = T)
        data_bias_adap <-  try(readRDS(file = paste0("bias_correct_", Boot_num,
                                                     "/",
                                                     link_type, "_", sample_size,
                                                     "/", 
                                                     group, 
                                                     "/", 
                                                     "adaptive_", grp, "_", sample_size, "_", link_type, "_", i, "_corr")), silent = T)
        
        
        Alpha_l_trt_point_orig[i] <- data_orig[[1]][1]
        Alpha_u_trt_point_orig[i] <- data_orig[[1]][2]
        
        Boot_est_u_adap <- rep(NA, Boot_num)
        Boot_est_l_adap <- rep(NA, Boot_num)
        
        # for onestep
        for(j in 1:Boot_num){
          
          temp <- data_bias_adap[[j]]
          Boot_est_l_adap[j] <- temp[[1]][1]
          Boot_est_u_adap[j] <- temp[[1]][2]
          
        }
        
        Alpha_u_trt_point_corr_adap[i] <- 2 * Alpha_u_trt_point_orig[i] - mean(Boot_est_u_adap, na.rm = T)
        Alpha_l_trt_point_corr_adap[i] <- 2 * Alpha_l_trt_point_orig[i] - mean(Boot_est_l_adap, na.rm = T)
        Alpha_u_trt_point_corr_NA_count_adap[i] <- length(which(is.na(Boot_est_u_adap) == T))
        Alpha_l_trt_point_corr_NA_count_adap[i] <- length(which(is.na(Boot_est_l_adap) == T))
        
        Convergence <- data_orig[[5]]
        Convergence_all[i] <- length(which(Convergence > 0))
        
      }
      
      
      max(Convergence_all)
      max(Alpha_u_trt_point_corr_NA_count_adap) 
      max(Alpha_l_trt_point_corr_NA_count_adap) 
      Alpha_u_trt_point_corr_NA_count_adap[which(Alpha_u_trt_point_corr_NA_count_adap > 0)]
      Alpha_l_trt_point_corr_NA_count_adap[which(Alpha_l_trt_point_corr_NA_count_adap > 0)]
      
      saveRDS(Alpha_u_trt_point_orig, paste0('value/orig/Alpha_u_', grp, '_orig_identity_',  sample_size))
      saveRDS(Alpha_l_trt_point_orig, paste0('value/orig/Alpha_l_', grp, '_orig_identity_', sample_size))
      
      saveRDS(Alpha_u_trt_point_corr_adap, paste0('value/bias_corr/Alpha_u_', grp, '_corr_identity_', sample_size))
      saveRDS(Alpha_l_trt_point_corr_adap, paste0('value/bias_corr/Alpha_l_', grp, '_corr_identity_', sample_size))
      
      mean(Alpha_u_trt_point_orig, na.rm = T) -  alpha_crl_u_true
      mean(Alpha_u_trt_point_corr_adap, na.rm = T) -  alpha_crl_u_true
      (mean(Alpha_u_trt_point_orig, na.rm = T) -  alpha_crl_u_true) /  alpha_crl_u_true
      (mean(Alpha_u_trt_point_corr_adap, na.rm = T) -  alpha_crl_u_true) /  alpha_crl_u_true
      
      mean(Alpha_l_trt_point_orig, na.rm = T) -  alpha_crl_l_true
      mean(Alpha_l_trt_point_corr_adap, na.rm = T) -  alpha_crl_l_true
      (mean(Alpha_l_trt_point_orig, na.rm = T) -  alpha_crl_l_true) /  alpha_crl_l_true
      (mean(Alpha_l_trt_point_corr_adap, na.rm = T) -  alpha_crl_l_true) /  alpha_crl_l_true
      
      
    }
    
    
  }
  

  ### plot for identity link (one-step)
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
    
    ## alpha_u
    Alpha_u_crl_point_orig <- readRDS("value/orig/Alpha_u_crl_orig_identity_200")
    Alpha_u_trt_point_orig <- readRDS("value/orig/Alpha_u_trt_orig_identity_200")
    
    Alpha_u_crl_point_corr <- readRDS("value/bias_corr/Alpha_u_crl_corr_identity_200")
    Alpha_u_trt_point_corr <- readRDS("value/bias_corr/Alpha_u_trt_corr_identity_200")
    
    
    ## alpha_l
    Alpha_l_crl_point_orig <- readRDS("value/orig/Alpha_l_crl_orig_identity_200")
    Alpha_l_trt_point_orig <- readRDS("value/orig/Alpha_l_trt_orig_identity_200")
    
    Alpha_l_crl_point_corr <- readRDS("value/bias_corr/Alpha_l_crl_corr_identity_200")
    Alpha_l_trt_point_corr <- readRDS("value/bias_corr/Alpha_l_trt_corr_identity_200")
    
    
    M <- 500
    data_plot_all <- data.frame(bias = c(Alpha_u_crl_point_orig - alpha_crl_u_true,
                                         Alpha_u_crl_point_corr - alpha_crl_u_true,
                                         Alpha_u_trt_point_orig - alpha_trt_u_true,
                                         Alpha_u_trt_point_corr - alpha_trt_u_true,
                                         Alpha_l_crl_point_orig - alpha_crl_l_true,
                                         Alpha_l_crl_point_corr - alpha_crl_l_true,
                                         Alpha_l_trt_point_orig - alpha_trt_l_true,
                                         Alpha_l_trt_point_corr - alpha_trt_l_true),
                                relative_bias =  c((Alpha_u_crl_point_orig - alpha_crl_u_true) / alpha_crl_u_true,
                                                   (Alpha_u_crl_point_corr - alpha_crl_u_true) / alpha_crl_u_true,
                                                   (Alpha_u_trt_point_orig - alpha_trt_u_true) / alpha_trt_u_true,
                                                   (Alpha_u_trt_point_corr - alpha_trt_u_true) / alpha_trt_u_true,
                                                   (Alpha_l_crl_point_orig - alpha_crl_l_true) / alpha_crl_l_true,
                                                   (Alpha_l_crl_point_corr - alpha_crl_l_true) / alpha_crl_l_true,
                                                   (Alpha_l_trt_point_orig - alpha_trt_l_true) / alpha_trt_l_true,
                                                   (Alpha_l_trt_point_corr - alpha_trt_l_true) / alpha_trt_l_true),
                                Type = c(rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M)),
                                # treatment = c(rep("Control Group (alpha_u)", 2 * M),
                                #               rep("Treatment Group (alpha_u)", 2 * M),
                                #               rep("Control Group (alpha_l)", 2 * M),
                                #               rep("Treatment Group (alpha_l)", 2 * M)),
                                treatment = c(rep("Control Group", 2 * M),
                                              rep("Treatment Group", 2 * M),
                                              rep("Control Group", 2 * M),
                                              rep("Treatment Group", 2 * M)),
                                Alpha_type = c(rep("Alpha_u", 4 * M),
                                               rep("Alpha_l", 4 * M)))
    
    data_plot_all$Type <- factor(data_plot_all$Type, 
                                 levels = c("Original",
                                            "Bias Correction"))
    
    data_plot_all$treatment <- as.factor(data_plot_all$treatment)
    data_plot_all$Alpha_type <- as.factor(data_plot_all$Alpha_type)
    
    
    plot3 <- ggplot(data = data_plot_all, aes(x = treatment, 
                                              y = relative_bias, 
                                              fill = Type)) +
      geom_boxplot(alpha = 0.3, linetype = "dashed", outlier.shape = NA) +
      stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
      stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
      stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) + 
      scale_fill_manual(values=c("white", 
                                 "grey"),
                                 # "grey50"),
                        breaks=c("Original",
                                 "Bias Correction")) +
      scale_y_continuous(breaks = seq(-1.25, 0.75, by=0.25)) +
      ylim(-1.25, 0.75) +
      expand_limits(y=c(-1.25, 0.75)) +
      geom_hline(yintercept = 0, linewidth = 0.5) + 
      facet_grid(. ~ Alpha_type) +
      theme_light()+
      xlab("Group")+
      ylab(expression(paste("Estimation Relative Bias of Alpha")))+
      theme(
        # legend.position="none",
        panel.grid.major =element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, color = "black"), 
        strip.background = element_rect(
          color="white", fill="white", linewidth = 1.5, linetype="solid")
      )+
      labs(fill='Estimation Type') +
      ggtitle("Relative Bias in Alpha_u and Alpha_l Before and After Bias Correction (Identity Link)")
    
    plot3
    
    ggsave(plot = plot3, 
           width = 9, 
           height = 5, 
           dpi = 500, 
           filename = "plots/bias_correct_identity.pdf")
    
    
  }
  
    
  ### plot for log link (weight1)
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
    
    ## alpha_u
    Alpha_u_crl_point_orig <- readRDS("value/orig/Alpha_u_crl_orig_log_weight1_200")
    Alpha_u_trt_point_orig <- readRDS("value/orig/Alpha_u_trt_orig_log_weight1_200")
    
    Alpha_u_crl_point_corr <- readRDS("value/bias_corr/Alpha_u_crl_corr_log_weight1_200")
    Alpha_u_trt_point_corr <- readRDS("value/bias_corr/Alpha_u_trt_corr_log_weight1_200")
    
    
    ## alpha_l
    Alpha_l_crl_point_orig <- readRDS("value/orig/Alpha_l_crl_orig_log_weight1_200")
    Alpha_l_trt_point_orig <- readRDS("value/orig/Alpha_l_trt_orig_log_weight1_200")
    
    Alpha_l_crl_point_corr <- readRDS("value/bias_corr/Alpha_l_crl_corr_log_weight1_200")
    Alpha_l_trt_point_corr <- readRDS("value/bias_corr/Alpha_l_trt_corr_log_weight1_200")
    
    
    M <- 500
    data_plot_all <- data.frame(bias = c(Alpha_u_crl_point_orig - alpha_crl_u_true,
                                         Alpha_u_crl_point_corr - alpha_crl_u_true,
                                         Alpha_u_trt_point_orig - alpha_trt_u_true,
                                         Alpha_u_trt_point_corr - alpha_trt_u_true,
                                         Alpha_l_crl_point_orig - alpha_crl_l_true,
                                         Alpha_l_crl_point_corr - alpha_crl_l_true,
                                         Alpha_l_trt_point_orig - alpha_trt_l_true,
                                         Alpha_l_trt_point_corr - alpha_trt_l_true),
                                relative_bias =  c((Alpha_u_crl_point_orig - alpha_crl_u_true) / alpha_crl_u_true,
                                                   (Alpha_u_crl_point_corr - alpha_crl_u_true) / alpha_crl_u_true,
                                                   (Alpha_u_trt_point_orig - alpha_trt_u_true) / alpha_trt_u_true,
                                                   (Alpha_u_trt_point_corr - alpha_trt_u_true) / alpha_trt_u_true,
                                                   (Alpha_l_crl_point_orig - alpha_crl_l_true) / alpha_crl_l_true,
                                                   (Alpha_l_crl_point_corr - alpha_crl_l_true) / alpha_crl_l_true,
                                                   (Alpha_l_trt_point_orig - alpha_trt_l_true) / alpha_trt_l_true,
                                                   (Alpha_l_trt_point_corr - alpha_trt_l_true) / alpha_trt_l_true),
                                Type = c(rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M)),
                                # treatment = c(rep("Control Group (alpha_u)", 2 * M),
                                #               rep("Treatment Group (alpha_u)", 2 * M),
                                #               rep("Control Group (alpha_l)", 2 * M),
                                #               rep("Treatment Group (alpha_l)", 2 * M)),
                                treatment = c(rep("Control Group", 2 * M),
                                              rep("Treatment Group", 2 * M),
                                              rep("Control Group", 2 * M),
                                              rep("Treatment Group", 2 * M)),
                                Alpha_type = c(rep("Alpha_u", 4 * M),
                                               rep("Alpha_l", 4 * M)))
    
    data_plot_all$Type <- factor(data_plot_all$Type, 
                                 levels = c("Original",
                                            "Bias Correction"))
    
    data_plot_all$treatment <- as.factor(data_plot_all$treatment)
    data_plot_all$Alpha_type <- as.factor(data_plot_all$Alpha_type)
    
    
    plot1 <- ggplot(data = data_plot_all, aes(x = treatment, 
                                              y = relative_bias, 
                                              fill = Type)) +
      geom_boxplot(alpha = 0.3, linetype = "dashed", outlier.shape = NA) +
      stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
      stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
      stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) + 
      scale_fill_manual(values=c("white", 
                                 "grey"),
                        # "grey50"),
                        breaks=c("Original",
                                 "Bias Correction")) +
      scale_y_continuous(breaks = seq(-1.25, 0.75, by=0.25)) +
      ylim(-1.25, 0.75) +
      expand_limits(y=c(-1.25, 0.75)) +
      geom_hline(yintercept = 0, linewidth = 0.5) + 
      facet_grid(. ~ Alpha_type) +
      theme_light()+
      xlab("Group")+
      ylab(expression(paste("Estimation Relative Bias of Alpha")))+
      theme(
        # legend.position="none",
        panel.grid.major =element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, color = "black"), 
        strip.background = element_rect(
          color="white", fill="white", linewidth = 1.5, linetype="solid")
      )+
      labs(fill='Estimation Type') +
      ggtitle("Relative Bias in Alpha_u and Alpha_l Before and After Bias Correction (Log Link - Weight1)")
    
    plot1
    
    ggsave(plot = plot1, 
           width = 9, 
           height = 5, 
           dpi = 500, 
           filename = "plots/bias_correct_log_weight1.pdf")
    
    
  }
  
  
  ### plot for log link (weight2)
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
    
    ## alpha_u
    Alpha_u_crl_point_orig <- readRDS("value/orig/Alpha_u_crl_orig_log_weight2_200")
    Alpha_u_trt_point_orig <- readRDS("value/orig/Alpha_u_trt_orig_log_weight2_200")
    
    Alpha_u_crl_point_corr <- readRDS("value/bias_corr/Alpha_u_crl_corr_log_weight2_200")
    Alpha_u_trt_point_corr <- readRDS("value/bias_corr/Alpha_u_trt_corr_log_weight2_200")
    
    
    ## alpha_l
    Alpha_l_crl_point_orig <- readRDS("value/orig/Alpha_l_crl_orig_log_weight2_200")
    Alpha_l_trt_point_orig <- readRDS("value/orig/Alpha_l_trt_orig_log_weight2_200")
    
    Alpha_l_crl_point_corr <- readRDS("value/bias_corr/Alpha_l_crl_corr_log_weight2_200")
    Alpha_l_trt_point_corr <- readRDS("value/bias_corr/Alpha_l_trt_corr_log_weight2_200")
    
    
    M <- 500
    data_plot_all <- data.frame(bias = c(Alpha_u_crl_point_orig - alpha_crl_u_true,
                                         Alpha_u_crl_point_corr - alpha_crl_u_true,
                                         Alpha_u_trt_point_orig - alpha_trt_u_true,
                                         Alpha_u_trt_point_corr - alpha_trt_u_true,
                                         Alpha_l_crl_point_orig - alpha_crl_l_true,
                                         Alpha_l_crl_point_corr - alpha_crl_l_true,
                                         Alpha_l_trt_point_orig - alpha_trt_l_true,
                                         Alpha_l_trt_point_corr - alpha_trt_l_true),
                                relative_bias =  c((Alpha_u_crl_point_orig - alpha_crl_u_true) / alpha_crl_u_true,
                                                   (Alpha_u_crl_point_corr - alpha_crl_u_true) / alpha_crl_u_true,
                                                   (Alpha_u_trt_point_orig - alpha_trt_u_true) / alpha_trt_u_true,
                                                   (Alpha_u_trt_point_corr - alpha_trt_u_true) / alpha_trt_u_true,
                                                   (Alpha_l_crl_point_orig - alpha_crl_l_true) / alpha_crl_l_true,
                                                   (Alpha_l_crl_point_corr - alpha_crl_l_true) / alpha_crl_l_true,
                                                   (Alpha_l_trt_point_orig - alpha_trt_l_true) / alpha_trt_l_true,
                                                   (Alpha_l_trt_point_corr - alpha_trt_l_true) / alpha_trt_l_true),
                                Type = c(rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M)),
                                # treatment = c(rep("Control Group (alpha_u)", 2 * M),
                                #               rep("Treatment Group (alpha_u)", 2 * M),
                                #               rep("Control Group (alpha_l)", 2 * M),
                                #               rep("Treatment Group (alpha_l)", 2 * M)),
                                treatment = c(rep("Control Group", 2 * M),
                                              rep("Treatment Group", 2 * M),
                                              rep("Control Group", 2 * M),
                                              rep("Treatment Group", 2 * M)),
                                Alpha_type = c(rep("Alpha_u", 4 * M),
                                               rep("Alpha_l", 4 * M)))
    
    data_plot_all$Type <- factor(data_plot_all$Type, 
                                 levels = c("Original",
                                            "Bias Correction"))
    
    data_plot_all$treatment <- as.factor(data_plot_all$treatment)
    data_plot_all$Alpha_type <- as.factor(data_plot_all$Alpha_type)
    
    
    plot2 <- ggplot(data = data_plot_all, aes(x = treatment, 
                                              y = relative_bias, 
                                              fill = Type)) +
      geom_boxplot(alpha = 0.3, linetype = "dashed", outlier.shape = NA) +
      stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
      stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
      stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) + 
      scale_fill_manual(values=c("white", 
                                 "grey"),
                        # "grey50"),
                        breaks=c("Original",
                                 "Bias Correction")) +
      scale_y_continuous(breaks = seq(-1.25, 0.75, by=0.25)) +
      ylim(-1.25, 0.75) +
      expand_limits(y=c(-1.25, 0.75)) +
      geom_hline(yintercept = 0, linewidth = 0.5) + 
      facet_grid(. ~ Alpha_type) +
      theme_light()+
      xlab("Group")+
      ylab(expression(paste("Estimation Relative Bias of Alpha")))+
      theme(
        # legend.position="none",
        panel.grid.major =element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, color = "black"), 
        strip.background = element_rect(
          color="white", fill="white", linewidth = 1.5, linetype="solid")
      )+
      labs(fill='Estimation Type') +
      ggtitle("Relative Bias in Alpha_u and Alpha_l Before and After Bias Correction (Log Link - Weight2)")
    
    plot2
    
    ggsave(plot = plot2, 
           width = 9, 
           height = 5, 
           dpi = 500, 
           filename = "plots/bias_correct_log_weight2.pdf")
    
    
  }
  
  
  
}
  


##### Point estimation  - logistic regression model simulated data
{
  
  
  # below is the true value
  alpha_crl_l_true <- -0.2632813
  alpha_crl_u_true <-  0.9632812
  alpha_trt_l_true <- -0.5685547
  alpha_trt_u_true <-  0.7760938
  
  
  #### One-step multiplier bootstrap - logit link
  {
    
    
    ## Treat group
    {
      
      link_type = "logit"
      sample_size = 200
      group = "treat"
      IF_type = "weight1"
      
      if(group == "treat"){
        grp = "trt"
      }else{
        grp = "crl"
      }
      Boot_num <- 100
      
      
      simu_size <- 500
      
      Alpha_u_trt_point_orig <- rep(NA, simu_size)
      Alpha_l_trt_point_orig <- rep(NA, simu_size) 
      Alpha_u_trt_point_corr_adap <- rep(NA, simu_size)
      Alpha_l_trt_point_corr_adap <- rep(NA, simu_size) 
      
      Alpha_u_trt_point_corr_NA_count_adap <- rep(NA, simu_size)
      Alpha_l_trt_point_corr_NA_count_adap <- rep(NA, simu_size)
      
      Convergence_all <- rep(NA, simu_size)
      
      setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
      
      for(i in 1:simu_size){
        
        print(i)
        data_orig <- try(readRDS(file = paste0('orig_point_est/', 
                                               link_type, "_", sample_size,
                                               "/",  
                                               group, "_", IF_type, 
                                               "/select_", grp, "_", sample_size, "_", link_type, "_", IF_type, "_", i)), silent = T)
        data_bias_adap <-  try(readRDS(file = paste0("bias_correct_", Boot_num,
                                                     "/",
                                                     link_type, "_", sample_size,
                                                     "/", 
                                                     group, "_", IF_type, 
                                                     "/", 
                                                     "adaptive_", grp, "_", sample_size, "_", link_type, "_", IF_type,"_", i, "_corr")), silent = T)
        
        
        Alpha_l_trt_point_orig[i] <- data_orig[[1]][1]
        Alpha_u_trt_point_orig[i] <- data_orig[[1]][2]
        
        Boot_est_u_adap <- rep(NA, Boot_num)
        Boot_est_l_adap <- rep(NA, Boot_num)
        
        # for onestep
        for(j in 1:Boot_num){
          
          temp <- data_bias_adap[[j]]
          Boot_est_l_adap[j] <- temp[[1]][1]
          Boot_est_u_adap[j] <- temp[[1]][2]
          
        }
        
        Alpha_u_trt_point_corr_adap[i] <- 2 * Alpha_u_trt_point_orig[i] - mean(Boot_est_u_adap, na.rm = T)
        Alpha_l_trt_point_corr_adap[i] <- 2 * Alpha_l_trt_point_orig[i] - mean(Boot_est_l_adap, na.rm = T)
        Alpha_u_trt_point_corr_NA_count_adap[i] <- length(which(is.na(Boot_est_u_adap) == T))
        Alpha_l_trt_point_corr_NA_count_adap[i] <- length(which(is.na(Boot_est_l_adap) == T))
        
        Convergence <- data_orig[[5]]
        Convergence_all[i] <- length(which(Convergence > 0))
        
      }
      
      max(Convergence_all)
      max(Alpha_u_trt_point_corr_NA_count_adap) 
      max(Alpha_l_trt_point_corr_NA_count_adap) 
      Alpha_u_trt_point_corr_NA_count_adap[which(Alpha_u_trt_point_corr_NA_count_adap > 0)]
      Alpha_l_trt_point_corr_NA_count_adap[which(Alpha_l_trt_point_corr_NA_count_adap > 0)]
      
      saveRDS(Alpha_u_trt_point_orig, paste0('value/orig/Alpha_u_', grp, '_orig_logit_', IF_type, '_', sample_size))
      saveRDS(Alpha_l_trt_point_orig, paste0('value/orig/Alpha_l_', grp, '_orig_logit_', IF_type, '_', sample_size))
      
      saveRDS(Alpha_u_trt_point_corr_adap, paste0('value/bias_corr/Alpha_u_', grp, '_corr_logit_', IF_type, '_', sample_size))
      saveRDS(Alpha_l_trt_point_corr_adap, paste0('value/bias_corr/Alpha_l_', grp, '_corr_logit_', IF_type, '_', sample_size))
      
      mean(Alpha_u_trt_point_orig, na.rm = T) -  alpha_trt_u_true
      mean(Alpha_u_trt_point_corr_adap, na.rm = T) -  alpha_trt_u_true
      (mean(Alpha_u_trt_point_orig, na.rm = T) -  alpha_trt_u_true) /  alpha_trt_u_true
      (mean(Alpha_u_trt_point_corr_adap, na.rm = T) -  alpha_trt_u_true) /  alpha_trt_u_true
      
      mean(Alpha_l_trt_point_orig, na.rm = T) -  alpha_trt_l_true
      mean(Alpha_l_trt_point_corr_adap, na.rm = T) -  alpha_trt_l_true
      (mean(Alpha_l_trt_point_orig, na.rm = T) -  alpha_trt_l_true) /  alpha_trt_l_true
      (mean(Alpha_l_trt_point_corr_adap, na.rm = T) -  alpha_trt_l_true) /  alpha_trt_l_true
      
      
    }
    
    
    ## Control group
    {
      
      link_type = "logit"
      sample_size = 200
      group = "control"
      IF_type = "weight1"
      
      if(group == "treat"){
        grp = "trt"
      }else{
        grp = "crl"
      }
      Boot_num <- 100
      
      
      simu_size <- 500
      
      Alpha_u_trt_point_orig <- rep(NA, simu_size)
      Alpha_l_trt_point_orig <- rep(NA, simu_size) 
      Alpha_u_trt_point_corr_adap <- rep(NA, simu_size)
      Alpha_l_trt_point_corr_adap <- rep(NA, simu_size) 
      
      Alpha_u_trt_point_corr_NA_count_adap <- rep(NA, simu_size)
      Alpha_l_trt_point_corr_NA_count_adap <- rep(NA, simu_size)
      
      Convergence_all <- rep(NA, simu_size)
      
      setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
      
      for(i in 1:simu_size){
        
        print(i)
        data_orig <- try(readRDS(file = paste0('orig_point_est/', 
                                               link_type, "_", sample_size,
                                               "/",  
                                               group, "_", IF_type, 
                                               "/select_", grp, "_", sample_size, "_", link_type, "_", IF_type, "_", i)), silent = T)
        data_bias_adap <-  try(readRDS(file = paste0("bias_correct_", Boot_num,
                                                     "/",
                                                     link_type, "_", sample_size,
                                                     "/", 
                                                     group, "_", IF_type, 
                                                     "/", 
                                                     "adaptive_", grp, "_", sample_size, "_", link_type, "_", IF_type,"_", i, "_corr")), silent = T)
        
        
        Alpha_l_trt_point_orig[i] <- data_orig[[1]][1]
        Alpha_u_trt_point_orig[i] <- data_orig[[1]][2]
        
        Boot_est_u_adap <- rep(NA, Boot_num)
        Boot_est_l_adap <- rep(NA, Boot_num)
        
        # for onestep
        for(j in 1:Boot_num){
          
          temp <- data_bias_adap[[j]]
          Boot_est_l_adap[j] <- temp[[1]][1]
          Boot_est_u_adap[j] <- temp[[1]][2]
          
        }
        
        Alpha_u_trt_point_corr_adap[i] <- 2 * Alpha_u_trt_point_orig[i] - mean(Boot_est_u_adap, na.rm = T)
        Alpha_l_trt_point_corr_adap[i] <- 2 * Alpha_l_trt_point_orig[i] - mean(Boot_est_l_adap, na.rm = T)
        Alpha_u_trt_point_corr_NA_count_adap[i] <- length(which(is.na(Boot_est_u_adap) == T))
        Alpha_l_trt_point_corr_NA_count_adap[i] <- length(which(is.na(Boot_est_l_adap) == T))
        
        Convergence <- data_orig[[5]]
        Convergence_all[i] <- length(which(Convergence > 0))
        
      }
      
      
      max(Convergence_all)
      max(Alpha_u_trt_point_corr_NA_count_adap) 
      max(Alpha_l_trt_point_corr_NA_count_adap) 
      Alpha_u_trt_point_corr_NA_count_adap[which(Alpha_u_trt_point_corr_NA_count_adap > 0)]
      Alpha_l_trt_point_corr_NA_count_adap[which(Alpha_l_trt_point_corr_NA_count_adap > 0)]
      
      saveRDS(Alpha_u_trt_point_orig, paste0('value/orig/Alpha_u_', grp, '_orig_logit_', IF_type, '_', sample_size))
      saveRDS(Alpha_l_trt_point_orig, paste0('value/orig/Alpha_l_', grp, '_orig_logit_', IF_type, '_', sample_size))
      
      saveRDS(Alpha_u_trt_point_corr_adap, paste0('value/bias_corr/Alpha_u_', grp, '_corr_logit_', IF_type, '_', sample_size))
      saveRDS(Alpha_l_trt_point_corr_adap, paste0('value/bias_corr/Alpha_l_', grp, '_corr_logit_', IF_type, '_', sample_size))
      
      mean(Alpha_u_trt_point_orig, na.rm = T) -  alpha_crl_u_true
      mean(Alpha_u_trt_point_corr_adap, na.rm = T) -  alpha_crl_u_true
      (mean(Alpha_u_trt_point_orig, na.rm = T) -  alpha_crl_u_true) /  alpha_crl_u_true
      (mean(Alpha_u_trt_point_corr_adap, na.rm = T) -  alpha_crl_u_true) /  alpha_crl_u_true
      
      mean(Alpha_l_trt_point_orig, na.rm = T) -  alpha_crl_l_true
      mean(Alpha_l_trt_point_corr_adap, na.rm = T) -  alpha_crl_l_true
      (mean(Alpha_l_trt_point_orig, na.rm = T) -  alpha_crl_l_true) /  alpha_crl_l_true
      (mean(Alpha_l_trt_point_corr_adap, na.rm = T) -  alpha_crl_l_true) /  alpha_crl_l_true
      
      
    }
    
    
  }
  
  
  ### plot for logit link (weight2)
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
    
    ## alpha_u
    Alpha_u_crl_point_orig <- readRDS("value/orig/Alpha_u_crl_orig_logit_weight2_200")
    Alpha_u_trt_point_orig <- readRDS("value/orig/Alpha_u_trt_orig_logit_weight2_200")
    
    Alpha_u_crl_point_corr <- readRDS("value/bias_corr/Alpha_u_crl_corr_logit_weight2_200")
    Alpha_u_trt_point_corr <- readRDS("value/bias_corr/Alpha_u_trt_corr_logit_weight2_200")
    
    
    ## alpha_l
    Alpha_l_crl_point_orig <- readRDS("value/orig/Alpha_l_crl_orig_logit_weight2_200")
    Alpha_l_trt_point_orig <- readRDS("value/orig/Alpha_l_trt_orig_logit_weight2_200")
    
    Alpha_l_crl_point_corr <- readRDS("value/bias_corr/Alpha_l_crl_corr_logit_weight2_200")
    Alpha_l_trt_point_corr <- readRDS("value/bias_corr/Alpha_l_trt_corr_logit_weight2_200")
    
    
    M <- 500
    data_plot_all <- data.frame(bias = c(Alpha_u_crl_point_orig - alpha_crl_u_true,
                                         Alpha_u_crl_point_corr - alpha_crl_u_true,
                                         Alpha_u_trt_point_orig - alpha_trt_u_true,
                                         Alpha_u_trt_point_corr - alpha_trt_u_true,
                                         Alpha_l_crl_point_orig - alpha_crl_l_true,
                                         Alpha_l_crl_point_corr - alpha_crl_l_true,
                                         Alpha_l_trt_point_orig - alpha_trt_l_true,
                                         Alpha_l_trt_point_corr - alpha_trt_l_true),
                                relative_bias =  c((Alpha_u_crl_point_orig - alpha_crl_u_true) / alpha_crl_u_true,
                                                   (Alpha_u_crl_point_corr - alpha_crl_u_true) / alpha_crl_u_true,
                                                   (Alpha_u_trt_point_orig - alpha_trt_u_true) / alpha_trt_u_true,
                                                   (Alpha_u_trt_point_corr - alpha_trt_u_true) / alpha_trt_u_true,
                                                   (Alpha_l_crl_point_orig - alpha_crl_l_true) / alpha_crl_l_true,
                                                   (Alpha_l_crl_point_corr - alpha_crl_l_true) / alpha_crl_l_true,
                                                   (Alpha_l_trt_point_orig - alpha_trt_l_true) / alpha_trt_l_true,
                                                   (Alpha_l_trt_point_corr - alpha_trt_l_true) / alpha_trt_l_true),
                                Type = c(rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M)),
                                # treatment = c(rep("Control Group (alpha_u)", 2 * M),
                                #               rep("Treatment Group (alpha_u)", 2 * M),
                                #               rep("Control Group (alpha_l)", 2 * M),
                                #               rep("Treatment Group (alpha_l)", 2 * M)),
                                treatment = c(rep("Control Group", 2 * M),
                                              rep("Treatment Group", 2 * M),
                                              rep("Control Group", 2 * M),
                                              rep("Treatment Group", 2 * M)),
                                Alpha_type = c(rep("Alpha_u", 4 * M),
                                               rep("Alpha_l", 4 * M)))
    
    data_plot_all$Type <- factor(data_plot_all$Type, 
                                 levels = c("Original",
                                            "Bias Correction"))
    
    data_plot_all$treatment <- as.factor(data_plot_all$treatment)
    data_plot_all$Alpha_type <- as.factor(data_plot_all$Alpha_type)
    
    
    plot4 <- ggplot(data = data_plot_all, aes(x = treatment, 
                                              y = relative_bias, 
                                              fill = Type)) +
      geom_boxplot(alpha = 0.3, linetype = "dashed", outlier.shape = NA) +
      stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
      stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
      stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) + 
      scale_fill_manual(values=c("white", 
                                 "grey"),
                        # "grey50"),
                        breaks=c("Original",
                                 "Bias Correction")) +
      scale_y_continuous(breaks = seq(-1.25, 0.75, by=0.25)) +
      ylim(-1.25, 0.75) +
      expand_limits(y=c(-1.25, 0.75)) +
      geom_hline(yintercept = 0, linewidth = 0.5) + 
      facet_grid(. ~ Alpha_type) +
      theme_light()+
      xlab("Group")+
      ylab(expression(paste("Estimation Relative Bias of Alpha")))+
      theme(
        # legend.position="none",
        panel.grid.major =element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, color = "black"), 
        strip.background = element_rect(
          color="white", fill="white", linewidth = 1.5, linetype="solid")
      )+
      labs(fill='Estimation Type') +
      ggtitle("Relative Bias in Alpha_u and Alpha_l Before and After Bias Correction (Logit Link - Weight2)")
    
    plot4
    
    ggsave(plot = plot4, 
           width = 9, 
           height = 5, 
           dpi = 500, 
           filename = "plots/bias_correct_logit_weight2.pdf")
    
    
  }
  
  
  ### plot for logit link (weight1)
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
    
    ## alpha_u
    Alpha_u_crl_point_orig <- readRDS("value/orig/Alpha_u_crl_orig_logit_weight1_200")
    Alpha_u_trt_point_orig <- readRDS("value/orig/Alpha_u_trt_orig_logit_weight1_200")
    
    Alpha_u_crl_point_corr <- readRDS("value/bias_corr/Alpha_u_crl_corr_logit_weight1_200")
    Alpha_u_trt_point_corr <- readRDS("value/bias_corr/Alpha_u_trt_corr_logit_weight1_200")
    
    
    ## alpha_l
    Alpha_l_crl_point_orig <- readRDS("value/orig/Alpha_l_crl_orig_logit_weight1_200")
    Alpha_l_trt_point_orig <- readRDS("value/orig/Alpha_l_trt_orig_logit_weight1_200")
    
    Alpha_l_crl_point_corr <- readRDS("value/bias_corr/Alpha_l_crl_corr_logit_weight1_200")
    Alpha_l_trt_point_corr <- readRDS("value/bias_corr/Alpha_l_trt_corr_logit_weight1_200")
    
    
    M <- 500
    data_plot_all <- data.frame(bias = c(Alpha_u_crl_point_orig - alpha_crl_u_true,
                                         Alpha_u_crl_point_corr - alpha_crl_u_true,
                                         Alpha_u_trt_point_orig - alpha_trt_u_true,
                                         Alpha_u_trt_point_corr - alpha_trt_u_true,
                                         Alpha_l_crl_point_orig - alpha_crl_l_true,
                                         Alpha_l_crl_point_corr - alpha_crl_l_true,
                                         Alpha_l_trt_point_orig - alpha_trt_l_true,
                                         Alpha_l_trt_point_corr - alpha_trt_l_true),
                                relative_bias =  c((Alpha_u_crl_point_orig - alpha_crl_u_true) / alpha_crl_u_true,
                                                   (Alpha_u_crl_point_corr - alpha_crl_u_true) / alpha_crl_u_true,
                                                   (Alpha_u_trt_point_orig - alpha_trt_u_true) / alpha_trt_u_true,
                                                   (Alpha_u_trt_point_corr - alpha_trt_u_true) / alpha_trt_u_true,
                                                   (Alpha_l_crl_point_orig - alpha_crl_l_true) / alpha_crl_l_true,
                                                   (Alpha_l_crl_point_corr - alpha_crl_l_true) / alpha_crl_l_true,
                                                   (Alpha_l_trt_point_orig - alpha_trt_l_true) / alpha_trt_l_true,
                                                   (Alpha_l_trt_point_corr - alpha_trt_l_true) / alpha_trt_l_true),
                                Type = c(rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M),
                                         rep("Original", M), rep("Bias Correction", M)),
                                # treatment = c(rep("Control Group (alpha_u)", 2 * M),
                                #               rep("Treatment Group (alpha_u)", 2 * M),
                                #               rep("Control Group (alpha_l)", 2 * M),
                                #               rep("Treatment Group (alpha_l)", 2 * M)),
                                treatment = c(rep("Control Group", 2 * M),
                                              rep("Treatment Group", 2 * M),
                                              rep("Control Group", 2 * M),
                                              rep("Treatment Group", 2 * M)),
                                Alpha_type = c(rep("Alpha_u", 4 * M),
                                               rep("Alpha_l", 4 * M)))
    
    data_plot_all$Type <- factor(data_plot_all$Type, 
                                 levels = c("Original",
                                            "Bias Correction"))
    
    data_plot_all$treatment <- as.factor(data_plot_all$treatment)
    data_plot_all$Alpha_type <- as.factor(data_plot_all$Alpha_type)
    
    
    plot5 <- ggplot(data = data_plot_all, aes(x = treatment, 
                                              y = relative_bias, 
                                              fill = Type)) +
      geom_boxplot(alpha = 0.3, linetype = "dashed", outlier.shape = NA) +
      stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
      stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
      stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) + 
      scale_fill_manual(values=c("white", 
                                 "grey"),
                        # "grey50"),
                        breaks=c("Original",
                                 "Bias Correction")) +
      scale_y_continuous(breaks = seq(-1.25, 0.75, by=0.25)) +
      ylim(-1.25, 0.75) +
      expand_limits(y=c(-1.25, 0.75)) +
      geom_hline(yintercept = 0, linewidth = 0.5) + 
      facet_grid(. ~ Alpha_type) +
      theme_light()+
      xlab("Group")+
      ylab(expression(paste("Estimation Relative Bias of Alpha")))+
      theme(
        # legend.position="none",
        panel.grid.major =element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, color = "black"), 
        strip.background = element_rect(
          color="white", fill="white", linewidth = 1.5, linetype="solid")
      )+
      labs(fill='Estimation Type') +
      ggtitle("Relative Bias in Alpha_u and Alpha_l Before and After Bias Correction (Logit Link - Weight1)")
    
    plot5
    
    ggsave(plot = plot5, 
           width = 9, 
           height = 5, 
           dpi = 500, 
           filename = "plots/bias_correct_logit_weight1.pdf")
    
    
  }
  
  
}



### Check coverage  - negative binomial model simulated data
{
  
  
  # below is the true value
  
  alpha_crl_l_true <- -0.7492188
  alpha_crl_u_true <- 0.4492188
  alpha_trt_l_true <- -0.9585937
  alpha_trt_u_true <- 0.6254688
  
  
  ## Log link - treatment group, weight1, weight2 done
  {
    
    
    link_type = "log"
    sample_size = 200
    group = "treat"
    IF_type = "weight2"
    
    if(group == "treat"){
      grp = "trt"
    }else{
      grp = "crl"
    }
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
    
    M <- 500
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_corr1 <- rep(NA, M)
    Cov_perc_alpha_u_trt_corr1 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_corr1 <-  rep(NA, M)
    Cov_alpha_CI_corr1 <- rep(NA, M)
    
    CI_perc_upper_alpha_u_trt_orig1 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_orig1 <- rep(NA, M)
    Cov_perc_alpha_u_trt_orig1 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_orig1 <-  rep(NA, M)
    Cov_alpha_CI_orig1 <- rep(NA, M)
    
    # Wald CI
    CI_perc_upper_alpha_u_trt_corr2 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_corr2 <- rep(NA, M)
    Cov_perc_alpha_u_trt_corr2 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_corr2 <-  rep(NA, M)
    Cov_alpha_CI_corr2 <- rep(NA, M)
    
    CI_perc_upper_alpha_u_trt_orig2 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_orig2 <- rep(NA, M)
    Cov_perc_alpha_u_trt_orig2 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_orig2 <-  rep(NA, M)
    Cov_alpha_CI_orig2 <- rep(NA, M)
    
    NA_alpha_u_avg <- rep(NA, M)
    NA_alpha_l_avg <- rep(NA, M)
    
    B_outer <- 400
    B_inner <- 50
    
    # Regarding these simulated data in log, the original optimization process cannot converge
    # Thus, only the outer bootstrap has the results, the inner bootstrap all failed.
    
    for(i in c(1:M)){
        
        print(i)
        
        Orig <- try(readRDS(file = paste0('orig_point_est/', 
                                             link_type, "_", sample_size,
                                             "/",  
                                             group, "_", IF_type, 
                                             "/select_", grp, "_", sample_size, "_", link_type, "_", IF_type, "_", i)), silent = T)
        Double <-  try(readRDS(file = paste0('double_boot_400_50/',
                                             link_type, "_", sample_size,
                                             "/", 
                                             group, "_", IF_type, 
                                             "/", 
                                             "adaptive_", grp, "_", sample_size, "_", link_type, "_", IF_type,"_", i, "_double")), silent = T)
      
      
        if(length(Double) > 1){
          
          alpha_u_boot_outer_orig <- rep(NA, B_outer)
          alpha_l_boot_outer_orig <- rep(NA, B_outer)
          alpha_u_boot_outer_corr <- rep(NA, B_outer)
          alpha_l_boot_outer_corr <- rep(NA, B_outer)
          NA_count_u_i <-  rep(NA, B_outer)
          NA_count_l_i <-  rep(NA, B_outer)
          
          for(j in 1:B_outer){
            
            temp <-  Double[[j]]
            alpha_l_boot_outer_orig[j] <- temp[[1]][1]
            alpha_u_boot_outer_orig[j] <- temp[[1]][2]
            alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[2]], na.rm = T)
            alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[3]], na.rm = T)
            
            NA_count_u_i[j] <- temp[[4]][1]
            NA_count_l_i[j] <- temp[[4]][2]
            
          }
          
          Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
          Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
          
          # Percentile CI
          CI_perc_upper_alpha_u_trt_corr1[i] <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
          CI_perc_low_alpha_l_trt_corr1[i] <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
          
          CI_perc_upper_alpha_u_trt_orig1[i] <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
          CI_perc_low_alpha_l_trt_orig1[i] <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
          
          Cov_perc_alpha_u_trt_corr1[i] <- ifelse(alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_corr1[i], 1, 0)
          Cov_perc_alpha_l_trt_corr1[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_corr1[i], 1, 0)
          Cov_perc_alpha_u_trt_orig1[i] <- ifelse(alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_orig1[i], 1, 0)
          Cov_perc_alpha_l_trt_orig1[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_orig1[i], 1, 0)
          
          Cov_alpha_CI_corr1[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_corr1[i] & alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_corr1[i], 1, 0)
          Cov_alpha_CI_orig1[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_orig1[i] & alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_orig1[i], 1, 0)
          
          # Wald CI
          CI_perc_low_alpha_l_trt_corr2[i] <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
          CI_perc_upper_alpha_u_trt_corr2[i] <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
          
          CI_perc_low_alpha_l_trt_orig2[i] <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
          CI_perc_upper_alpha_u_trt_orig2[i] <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
          
          Cov_perc_alpha_u_trt_corr2[i] <- ifelse(alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_corr2[i], 1, 0)
          Cov_perc_alpha_l_trt_corr2[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_corr2[i], 1, 0)
          Cov_perc_alpha_u_trt_orig2[i] <- ifelse(alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_orig2[i], 1, 0)
          Cov_perc_alpha_l_trt_orig2[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_orig2[i], 1, 0)
          
          Cov_alpha_CI_corr2[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_corr2[i] & alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_corr2[i], 1, 0)
          Cov_alpha_CI_orig2[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_orig2[i] & alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_orig2[i], 1, 0)
          
          NA_alpha_u_avg[i] <- mean(NA_count_u_i, na.rm = T)
          NA_alpha_l_avg[i] <- mean(NA_count_l_i, na.rm = T)
          
        }
        
      # }
      
    }
    
    # check NA values
    which(is.na(Cov_alpha_CI_corr1) == T)
    which(is.na(Cov_alpha_CI_corr1) == T) |> length()
    
    
    # Percentile CI
    c(mean(Cov_perc_alpha_l_trt_orig1, na.rm = T), 
      mean(Cov_perc_alpha_u_trt_orig1, na.rm = T), 
      mean(Cov_alpha_CI_orig1, na.rm = T))
    c( mean(Cov_perc_alpha_l_trt_corr1, na.rm = T),
       mean(Cov_perc_alpha_u_trt_corr1, na.rm = T),
       mean(Cov_alpha_CI_corr1, na.rm = T))
    
    
    # Wald CI
    c(mean(Cov_perc_alpha_l_trt_orig2, na.rm = T),
      mean(Cov_perc_alpha_u_trt_orig2, na.rm = T),
      mean(Cov_alpha_CI_orig2, na.rm = T))
    c(mean(Cov_perc_alpha_l_trt_corr2, na.rm = T),
      mean(Cov_perc_alpha_u_trt_corr2, na.rm = T),
      mean(Cov_alpha_CI_corr2, na.rm = T))
    
    
    NA_alpha_u_avg
    NA_alpha_l_avg
    
    
  }
  
  
  ## Log link - control group, weight1, weight2 done
  {
    
    
    link_type = "log"
    sample_size = 200
    group = "control"
    IF_type = "weight2"
    
    if(group == "treat"){
      grp = "trt"
    }else{
      grp = "crl"
    }
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
    
    M <- 500
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_corr1 <- rep(NA, M)
    Cov_perc_alpha_u_trt_corr1 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_corr1 <-  rep(NA, M)
    Cov_alpha_CI_corr1 <- rep(NA, M)
    
    CI_perc_upper_alpha_u_trt_orig1 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_orig1 <- rep(NA, M)
    Cov_perc_alpha_u_trt_orig1 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_orig1 <-  rep(NA, M)
    Cov_alpha_CI_orig1 <- rep(NA, M)
    
    # Wald CI
    CI_perc_upper_alpha_u_trt_corr2 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_corr2 <- rep(NA, M)
    Cov_perc_alpha_u_trt_corr2 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_corr2 <-  rep(NA, M)
    Cov_alpha_CI_corr2 <- rep(NA, M)
    
    CI_perc_upper_alpha_u_trt_orig2 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_orig2 <- rep(NA, M)
    Cov_perc_alpha_u_trt_orig2 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_orig2 <-  rep(NA, M)
    Cov_alpha_CI_orig2 <- rep(NA, M)
    
    NA_alpha_u_avg <- rep(NA, M)
    NA_alpha_l_avg <- rep(NA, M)
    
    B_outer <- 400
    B_inner <- 50
    
    # Regarding these simulated data in log, the original optimization process cannot converge
    # Thus, only the outer bootstrap has the results, the inner bootstrap all failed.
    
    for(i in c(1:M)){
      
      print(i)
      
      Orig <- try(readRDS(file = paste0('orig_point_est/', 
                                        link_type, "_", sample_size,
                                        "/",  
                                        group, "_", IF_type, 
                                        "/select_", grp, "_", sample_size, "_", link_type, "_", IF_type, "_", i)), silent = T)
      Double <-  try(readRDS(file = paste0('double_boot_400_50/',
                                           link_type, "_", sample_size,
                                           "/", 
                                           group, "_", IF_type, 
                                           "/", 
                                           "adaptive_", grp, "_", sample_size, "_", link_type, "_", IF_type,"_", i, "_double")), silent = T)
      
      
      if(length(Double) > 1){
        
        alpha_u_boot_outer_orig <- rep(NA, B_outer)
        alpha_l_boot_outer_orig <- rep(NA, B_outer)
        alpha_u_boot_outer_corr <- rep(NA, B_outer)
        alpha_l_boot_outer_corr <- rep(NA, B_outer)
        NA_count_u_i <-  rep(NA, B_outer)
        NA_count_l_i <-  rep(NA, B_outer)
        
        for(j in 1:B_outer){
          
          temp <-  Double[[j]]
          alpha_l_boot_outer_orig[j] <- temp[[1]][1]
          alpha_u_boot_outer_orig[j] <- temp[[1]][2]
          alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[2]], na.rm = T)
          alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[3]], na.rm = T)
          
          NA_count_u_i[j] <- temp[[4]][1]
          NA_count_l_i[j] <- temp[[4]][2]
          
        }
        
        Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
        Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
        
        # Percentile CI
        CI_perc_upper_alpha_u_trt_corr1[i] <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
        CI_perc_low_alpha_l_trt_corr1[i] <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
        
        CI_perc_upper_alpha_u_trt_orig1[i] <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
        CI_perc_low_alpha_l_trt_orig1[i] <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
        
        Cov_perc_alpha_u_trt_corr1[i] <- ifelse(alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_corr1[i], 1, 0)
        Cov_perc_alpha_l_trt_corr1[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_corr1[i], 1, 0)
        Cov_perc_alpha_u_trt_orig1[i] <- ifelse(alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_orig1[i], 1, 0)
        Cov_perc_alpha_l_trt_orig1[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_orig1[i], 1, 0)
        
        Cov_alpha_CI_corr1[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_corr1[i] & alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_corr1[i], 1, 0)
        Cov_alpha_CI_orig1[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_orig1[i] & alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_orig1[i], 1, 0)
        
        # Wald CI
        CI_perc_low_alpha_l_trt_corr2[i] <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
        CI_perc_upper_alpha_u_trt_corr2[i] <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
        
        CI_perc_low_alpha_l_trt_orig2[i] <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
        CI_perc_upper_alpha_u_trt_orig2[i] <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
        
        Cov_perc_alpha_u_trt_corr2[i] <- ifelse(alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_corr2[i], 1, 0)
        Cov_perc_alpha_l_trt_corr2[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_corr2[i], 1, 0)
        Cov_perc_alpha_u_trt_orig2[i] <- ifelse(alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_orig2[i], 1, 0)
        Cov_perc_alpha_l_trt_orig2[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_orig2[i], 1, 0)
        
        Cov_alpha_CI_corr2[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_corr2[i] & alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_corr2[i], 1, 0)
        Cov_alpha_CI_orig2[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_orig2[i] & alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_orig2[i], 1, 0)
        
        NA_alpha_u_avg[i] <- mean(NA_count_u_i, na.rm = T)
        NA_alpha_l_avg[i] <- mean(NA_count_l_i, na.rm = T)
        
      }
      
      # }
      
    }
    
    # check NA values
    which(is.na(Cov_alpha_CI_corr1) == T)
    which(is.na(Cov_alpha_CI_corr1) == T) |> length()
     
    
    # Percentile CI
    c(mean(Cov_perc_alpha_l_trt_orig1, na.rm = T), 
      mean(Cov_perc_alpha_u_trt_orig1, na.rm = T), 
      mean(Cov_alpha_CI_orig1, na.rm = T))
    c( mean(Cov_perc_alpha_l_trt_corr1, na.rm = T),
       mean(Cov_perc_alpha_u_trt_corr1, na.rm = T),
       mean(Cov_alpha_CI_corr1, na.rm = T))
    
    
    # Wald CI
    c(mean(Cov_perc_alpha_l_trt_orig2, na.rm = T),
      mean(Cov_perc_alpha_u_trt_orig2, na.rm = T),
      mean(Cov_alpha_CI_orig2, na.rm = T))
    c(mean(Cov_perc_alpha_l_trt_corr2, na.rm = T),
      mean(Cov_perc_alpha_u_trt_corr2, na.rm = T),
      mean(Cov_alpha_CI_corr2, na.rm = T))
    
    
    NA_alpha_u_avg
    NA_alpha_l_avg
    
  }
  
  
  ## Identity link - treatment group
  {
    
    
    link_type = "identity"
    sample_size = 200
    group = "treat"
    
    if(group == "treat"){
      grp = "trt"
    }else{
      grp = "crl"
    }
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
    
    M <- 500
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_corr1 <- rep(NA, M)
    Cov_perc_alpha_u_trt_corr1 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_corr1 <-  rep(NA, M)
    Cov_alpha_CI_corr1 <- rep(NA, M)
    
    CI_perc_upper_alpha_u_trt_orig1 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_orig1 <- rep(NA, M)
    Cov_perc_alpha_u_trt_orig1 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_orig1 <-  rep(NA, M)
    Cov_alpha_CI_orig1 <- rep(NA, M)
    
    # Wald CI
    CI_perc_upper_alpha_u_trt_corr2 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_corr2 <- rep(NA, M)
    Cov_perc_alpha_u_trt_corr2 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_corr2 <-  rep(NA, M)
    Cov_alpha_CI_corr2 <- rep(NA, M)
    
    CI_perc_upper_alpha_u_trt_orig2 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_orig2 <- rep(NA, M)
    Cov_perc_alpha_u_trt_orig2 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_orig2 <-  rep(NA, M)
    Cov_alpha_CI_orig2 <- rep(NA, M)
    
    NA_alpha_u_avg <- rep(NA, M)
    NA_alpha_l_avg <- rep(NA, M)
    
    B_outer <- 400
    B_inner <- 50
    
    # Regarding these simulated data in log, the original optimization process cannot converge
    # Thus, only the outer bootstrap has the results, the inner bootstrap all failed.
    

    for(i in c(1:M)){
      
      print(i)
      
      Orig <- try(readRDS(file = paste0('orig_point_est/', 
                                        link_type, "_", sample_size,
                                        "/",  
                                        group, 
                                        "/select_", grp, "_", sample_size, "_", link_type, "_", i)), silent = T)
      Double <-  try(readRDS(file = paste0('double_boot_400_50/',
                                           link_type, "_", sample_size,
                                           "/", 
                                           group,  
                                           "/", 
                                           "adaptive_", grp, "_", sample_size, "_", link_type, "_", i, "_double")), silent = T)
      
      
      if(length(Double) > 1){
        
        alpha_u_boot_outer_orig <- rep(NA, B_outer)
        alpha_l_boot_outer_orig <- rep(NA, B_outer)
        alpha_u_boot_outer_corr <- rep(NA, B_outer)
        alpha_l_boot_outer_corr <- rep(NA, B_outer)
        NA_count_u_i <-  rep(NA, B_outer)
        NA_count_l_i <-  rep(NA, B_outer)
        
        for(j in 1:B_outer){
          
          temp <-  Double[[j]]
          alpha_l_boot_outer_orig[j] <- temp[[1]][1]
          alpha_u_boot_outer_orig[j] <- temp[[1]][2]
          alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[2]], na.rm = T)
          alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[3]], na.rm = T)
          
          NA_count_u_i[j] <- temp[[4]][1]
          NA_count_l_i[j] <- temp[[4]][2]
          
        }
        
        Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
        Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
        
        # Percentile CI
        CI_perc_upper_alpha_u_trt_corr1[i] <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
        CI_perc_low_alpha_l_trt_corr1[i] <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
        
        CI_perc_upper_alpha_u_trt_orig1[i] <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
        CI_perc_low_alpha_l_trt_orig1[i] <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
        
        Cov_perc_alpha_u_trt_corr1[i] <- ifelse(alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_corr1[i], 1, 0)
        Cov_perc_alpha_l_trt_corr1[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_corr1[i], 1, 0)
        Cov_perc_alpha_u_trt_orig1[i] <- ifelse(alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_orig1[i], 1, 0)
        Cov_perc_alpha_l_trt_orig1[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_orig1[i], 1, 0)
        
        Cov_alpha_CI_corr1[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_corr1[i] & alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_corr1[i], 1, 0)
        Cov_alpha_CI_orig1[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_orig1[i] & alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_orig1[i], 1, 0)
        
        # Wald CI
        CI_perc_low_alpha_l_trt_corr2[i] <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
        CI_perc_upper_alpha_u_trt_corr2[i] <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
        
        CI_perc_low_alpha_l_trt_orig2[i] <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
        CI_perc_upper_alpha_u_trt_orig2[i] <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
        
        Cov_perc_alpha_u_trt_corr2[i] <- ifelse(alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_corr2[i], 1, 0)
        Cov_perc_alpha_l_trt_corr2[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_corr2[i], 1, 0)
        Cov_perc_alpha_u_trt_orig2[i] <- ifelse(alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_orig2[i], 1, 0)
        Cov_perc_alpha_l_trt_orig2[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_orig2[i], 1, 0)
        
        Cov_alpha_CI_corr2[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_corr2[i] & alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_corr2[i], 1, 0)
        Cov_alpha_CI_orig2[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_orig2[i] & alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_orig2[i], 1, 0)
        
        NA_alpha_u_avg[i] <- mean(NA_count_u_i, na.rm = T)
        NA_alpha_l_avg[i] <- mean(NA_count_l_i, na.rm = T)
        
      }
      
      # }
      
    }
    
    # check NA values
    which(is.na(Cov_alpha_CI_corr1) == T)
    which(is.na(Cov_alpha_CI_corr1) == T) |> length()
    
    
    # Percentile CI
    c(mean(Cov_perc_alpha_l_trt_orig1, na.rm = T), 
      mean(Cov_perc_alpha_u_trt_orig1, na.rm = T), 
      mean(Cov_alpha_CI_orig1, na.rm = T))
    c( mean(Cov_perc_alpha_l_trt_corr1, na.rm = T),
       mean(Cov_perc_alpha_u_trt_corr1, na.rm = T),
       mean(Cov_alpha_CI_corr1, na.rm = T))
    
    
    # Wald CI
    c(mean(Cov_perc_alpha_l_trt_orig2, na.rm = T),
      mean(Cov_perc_alpha_u_trt_orig2, na.rm = T),
      mean(Cov_alpha_CI_orig2, na.rm = T))
    c(mean(Cov_perc_alpha_l_trt_corr2, na.rm = T),
      mean(Cov_perc_alpha_u_trt_corr2, na.rm = T),
      mean(Cov_alpha_CI_corr2, na.rm = T))
    
    
    NA_alpha_u_avg
    NA_alpha_l_avg
    
  }
  
  
  ## Identity link - control group
  {
    
    
    link_type = "identity"
    sample_size = 200
    group = "control"
    
    if(group == "treat"){
      grp = "trt"
    }else{
      grp = "crl"
    }
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
    
    M <- 500
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_corr1 <- rep(NA, M)
    Cov_perc_alpha_u_trt_corr1 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_corr1 <-  rep(NA, M)
    Cov_alpha_CI_corr1 <- rep(NA, M)
    
    CI_perc_upper_alpha_u_trt_orig1 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_orig1 <- rep(NA, M)
    Cov_perc_alpha_u_trt_orig1 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_orig1 <-  rep(NA, M)
    Cov_alpha_CI_orig1 <- rep(NA, M)
    
    # Wald CI
    CI_perc_upper_alpha_u_trt_corr2 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_corr2 <- rep(NA, M)
    Cov_perc_alpha_u_trt_corr2 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_corr2 <-  rep(NA, M)
    Cov_alpha_CI_corr2 <- rep(NA, M)
    
    CI_perc_upper_alpha_u_trt_orig2 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_orig2 <- rep(NA, M)
    Cov_perc_alpha_u_trt_orig2 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_orig2 <-  rep(NA, M)
    Cov_alpha_CI_orig2 <- rep(NA, M)
    
    NA_alpha_u_avg <- rep(NA, M)
    NA_alpha_l_avg <- rep(NA, M)
    
    B_outer <- 400
    B_inner <- 50
    
    # Regarding these simulated data in log, the original optimization process cannot converge
    # Thus, only the outer bootstrap has the results, the inner bootstrap all failed.
    
    
    for(i in c(1:M)){
      
      print(i)
      
      Orig <- try(readRDS(file = paste0('orig_point_est/', 
                                        link_type, "_", sample_size,
                                        "/",  
                                        group, 
                                        "/select_", grp, "_", sample_size, "_", link_type, "_", i)), silent = T)
      Double <-  try(readRDS(file = paste0('double_boot_400_50/',
                                           link_type, "_", sample_size,
                                           "/", 
                                           group,  
                                           "/", 
                                           "adaptive_", grp, "_", sample_size, "_", link_type, "_", i, "_double")), silent = T)
      
      
      if(length(Double) > 1){
        
        alpha_u_boot_outer_orig <- rep(NA, B_outer)
        alpha_l_boot_outer_orig <- rep(NA, B_outer)
        alpha_u_boot_outer_corr <- rep(NA, B_outer)
        alpha_l_boot_outer_corr <- rep(NA, B_outer)
        NA_count_u_i <-  rep(NA, B_outer)
        NA_count_l_i <-  rep(NA, B_outer)
        
        for(j in 1:B_outer){
          
          temp <-  Double[[j]]
          alpha_l_boot_outer_orig[j] <- temp[[1]][1]
          alpha_u_boot_outer_orig[j] <- temp[[1]][2]
          alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[2]], na.rm = T)
          alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[3]], na.rm = T)
          
          NA_count_u_i[j] <- temp[[4]][1]
          NA_count_l_i[j] <- temp[[4]][2]
          
        }
        
        Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
        Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
        
        # Percentile CI
        CI_perc_upper_alpha_u_trt_corr1[i] <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
        CI_perc_low_alpha_l_trt_corr1[i] <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
        
        CI_perc_upper_alpha_u_trt_orig1[i] <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
        CI_perc_low_alpha_l_trt_orig1[i] <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
        
        Cov_perc_alpha_u_trt_corr1[i] <- ifelse(alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_corr1[i], 1, 0)
        Cov_perc_alpha_l_trt_corr1[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_corr1[i], 1, 0)
        Cov_perc_alpha_u_trt_orig1[i] <- ifelse(alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_orig1[i], 1, 0)
        Cov_perc_alpha_l_trt_orig1[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_orig1[i], 1, 0)
        
        Cov_alpha_CI_corr1[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_corr1[i] & alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_corr1[i], 1, 0)
        Cov_alpha_CI_orig1[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_orig1[i] & alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_orig1[i], 1, 0)
        
        # Wald CI
        CI_perc_low_alpha_l_trt_corr2[i] <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
        CI_perc_upper_alpha_u_trt_corr2[i] <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
        
        CI_perc_low_alpha_l_trt_orig2[i] <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
        CI_perc_upper_alpha_u_trt_orig2[i] <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
        
        Cov_perc_alpha_u_trt_corr2[i] <- ifelse(alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_corr2[i], 1, 0)
        Cov_perc_alpha_l_trt_corr2[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_corr2[i], 1, 0)
        Cov_perc_alpha_u_trt_orig2[i] <- ifelse(alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_orig2[i], 1, 0)
        Cov_perc_alpha_l_trt_orig2[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_orig2[i], 1, 0)
        
        Cov_alpha_CI_corr2[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_corr2[i] & alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_corr2[i], 1, 0)
        Cov_alpha_CI_orig2[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_orig2[i] & alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_orig2[i], 1, 0)
        
        NA_alpha_u_avg[i] <- mean(NA_count_u_i, na.rm = T)
        NA_alpha_l_avg[i] <- mean(NA_count_l_i, na.rm = T)
        
      }
      
      # }
      
    }
    
    # check NA values
    which(is.na(Cov_alpha_CI_corr1) == T)
    which(is.na(Cov_alpha_CI_corr1) == T) |> length()
    
    
    # Percentile CI
    c(mean(Cov_perc_alpha_l_trt_orig1, na.rm = T), 
      mean(Cov_perc_alpha_u_trt_orig1, na.rm = T), 
      mean(Cov_alpha_CI_orig1, na.rm = T))
    c( mean(Cov_perc_alpha_l_trt_corr1, na.rm = T),
       mean(Cov_perc_alpha_u_trt_corr1, na.rm = T),
       mean(Cov_alpha_CI_corr1, na.rm = T))
    
    
    # Wald CI
    c(mean(Cov_perc_alpha_l_trt_orig2, na.rm = T),
      mean(Cov_perc_alpha_u_trt_orig2, na.rm = T),
      mean(Cov_alpha_CI_orig2, na.rm = T))
    c(mean(Cov_perc_alpha_l_trt_corr2, na.rm = T),
      mean(Cov_perc_alpha_u_trt_corr2, na.rm = T),
      mean(Cov_alpha_CI_corr2, na.rm = T))
    
    
    NA_alpha_u_avg
    NA_alpha_l_avg
    
  }
  
  
  # below is the true value for logit link
  alpha_crl_l_true <- -0.2632813
  alpha_crl_u_true <-  0.9632812
  alpha_trt_l_true <- -0.5685547
  alpha_trt_u_true <-  0.7760938
  

  ## Logit link - treatment group, weight1 done, 
  {
    
    
    link_type = "logit"
    sample_size = 200
    group = "treat"
    IF_type = "weight2"
    
    if(group == "treat"){
      grp = "trt"
    }else{
      grp = "crl"
    }
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
    
    M <- 500
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_corr1 <- rep(NA, M)
    Cov_perc_alpha_u_trt_corr1 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_corr1 <-  rep(NA, M)
    Cov_alpha_CI_corr1 <- rep(NA, M)
    
    CI_perc_upper_alpha_u_trt_orig1 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_orig1 <- rep(NA, M)
    Cov_perc_alpha_u_trt_orig1 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_orig1 <-  rep(NA, M)
    Cov_alpha_CI_orig1 <- rep(NA, M)
    
    # Wald CI
    CI_perc_upper_alpha_u_trt_corr2 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_corr2 <- rep(NA, M)
    Cov_perc_alpha_u_trt_corr2 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_corr2 <-  rep(NA, M)
    Cov_alpha_CI_corr2 <- rep(NA, M)
    
    CI_perc_upper_alpha_u_trt_orig2 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_orig2 <- rep(NA, M)
    Cov_perc_alpha_u_trt_orig2 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_orig2 <-  rep(NA, M)
    Cov_alpha_CI_orig2 <- rep(NA, M)
    
    NA_alpha_u_avg <- rep(NA, M)
    NA_alpha_l_avg <- rep(NA, M)
    
    B_outer <- 400
    B_inner <- 50
    
    # Regarding these simulated data in log, the original optimization process cannot converge
    # Thus, only the outer bootstrap has the results, the inner bootstrap all failed.
    
    for(i in c(1:M)){
      
      print(i)
      
      Orig <- try(readRDS(file = paste0('orig_point_est/', 
                                        link_type, "_", sample_size,
                                        "/",  
                                        group, "_", IF_type, 
                                        "/select_", grp, "_", sample_size, "_", link_type, "_", IF_type, "_", i)), silent = T)
      Double <-  try(readRDS(file = paste0('double_boot_400_50/',
                                           link_type, "_", sample_size,
                                           "/", 
                                           group, "_", IF_type, 
                                           "/", 
                                           "double_", grp, "_", sample_size, "_", link_type, "_", IF_type,"_", i, "_corr")), silent = T)
      
      
      if(length(Double) > 1){
        
        alpha_u_boot_outer_orig <- rep(NA, B_outer)
        alpha_l_boot_outer_orig <- rep(NA, B_outer)
        alpha_u_boot_outer_corr <- rep(NA, B_outer)
        alpha_l_boot_outer_corr <- rep(NA, B_outer)
        NA_count_u_i <-  rep(NA, B_outer)
        NA_count_l_i <-  rep(NA, B_outer)
        
        for(j in 1:B_outer){
          
          temp <-  Double[[j]]
          alpha_l_boot_outer_orig[j] <- temp[[1]][1]
          alpha_u_boot_outer_orig[j] <- temp[[1]][2]
          alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[2]], na.rm = T)
          alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[3]], na.rm = T)
          
          NA_count_u_i[j] <- temp[[4]][1]
          NA_count_l_i[j] <- temp[[4]][2]
          
        }
        
        Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
        Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
        
        # Percentile CI
        CI_perc_upper_alpha_u_trt_corr1[i] <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
        CI_perc_low_alpha_l_trt_corr1[i] <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
        
        CI_perc_upper_alpha_u_trt_orig1[i] <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
        CI_perc_low_alpha_l_trt_orig1[i] <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
        
        Cov_perc_alpha_u_trt_corr1[i] <- ifelse(alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_corr1[i], 1, 0)
        Cov_perc_alpha_l_trt_corr1[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_corr1[i], 1, 0)
        Cov_perc_alpha_u_trt_orig1[i] <- ifelse(alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_orig1[i], 1, 0)
        Cov_perc_alpha_l_trt_orig1[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_orig1[i], 1, 0)
        
        Cov_alpha_CI_corr1[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_corr1[i] & alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_corr1[i], 1, 0)
        Cov_alpha_CI_orig1[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_orig1[i] & alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_orig1[i], 1, 0)
        
        # Wald CI
        CI_perc_low_alpha_l_trt_corr2[i] <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
        CI_perc_upper_alpha_u_trt_corr2[i] <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
        
        CI_perc_low_alpha_l_trt_orig2[i] <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
        CI_perc_upper_alpha_u_trt_orig2[i] <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
        
        Cov_perc_alpha_u_trt_corr2[i] <- ifelse(alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_corr2[i], 1, 0)
        Cov_perc_alpha_l_trt_corr2[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_corr2[i], 1, 0)
        Cov_perc_alpha_u_trt_orig2[i] <- ifelse(alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_orig2[i], 1, 0)
        Cov_perc_alpha_l_trt_orig2[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_orig2[i], 1, 0)
        
        Cov_alpha_CI_corr2[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_corr2[i] & alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_corr2[i], 1, 0)
        Cov_alpha_CI_orig2[i] <- ifelse(alpha_trt_l_true >= CI_perc_low_alpha_l_trt_orig2[i] & alpha_trt_u_true <= CI_perc_upper_alpha_u_trt_orig2[i], 1, 0)
        
        NA_alpha_u_avg[i] <- mean(NA_count_u_i, na.rm = T)
        NA_alpha_l_avg[i] <- mean(NA_count_l_i, na.rm = T)
        
      }
      
      
    }
    
    # check NA values
    which(is.na(Cov_alpha_CI_corr1) == T)
    which(is.na(Cov_alpha_CI_corr1) == T) |> length()
    
    
    # Percentile CI
    c(mean(Cov_perc_alpha_l_trt_orig1, na.rm = T), 
      mean(Cov_perc_alpha_u_trt_orig1, na.rm = T), 
      mean(Cov_alpha_CI_orig1, na.rm = T))
    c( mean(Cov_perc_alpha_l_trt_corr1, na.rm = T),
       mean(Cov_perc_alpha_u_trt_corr1, na.rm = T),
       mean(Cov_alpha_CI_corr1, na.rm = T))
    
    
    # Wald CI
    c(mean(Cov_perc_alpha_l_trt_orig2, na.rm = T),
      mean(Cov_perc_alpha_u_trt_orig2, na.rm = T),
      mean(Cov_alpha_CI_orig2, na.rm = T))
    c(mean(Cov_perc_alpha_l_trt_corr2, na.rm = T),
      mean(Cov_perc_alpha_u_trt_corr2, na.rm = T),
      mean(Cov_alpha_CI_corr2, na.rm = T))
    
    
    NA_alpha_u_avg
    NA_alpha_l_avg
    
  }
  
  
  ## Logit link - control group, weight1 done, 
  {
    
    
    link_type = "logit"
    sample_size = 200
    group = "control"
    IF_type = "weight2"
    
    if(group == "treat"){
      grp = "trt"
    }else{
      grp = "crl"
    }
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new")
    
    M <- 500
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_corr1 <- rep(NA, M)
    Cov_perc_alpha_u_trt_corr1 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_corr1 <-  rep(NA, M)
    Cov_alpha_CI_corr1 <- rep(NA, M)
    
    CI_perc_upper_alpha_u_trt_orig1 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_orig1 <- rep(NA, M)
    Cov_perc_alpha_u_trt_orig1 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_orig1 <-  rep(NA, M)
    Cov_alpha_CI_orig1 <- rep(NA, M)
    
    # Wald CI
    CI_perc_upper_alpha_u_trt_corr2 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_corr2 <- rep(NA, M)
    Cov_perc_alpha_u_trt_corr2 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_corr2 <-  rep(NA, M)
    Cov_alpha_CI_corr2 <- rep(NA, M)
    
    CI_perc_upper_alpha_u_trt_orig2 <- rep(NA, M)
    CI_perc_low_alpha_l_trt_orig2 <- rep(NA, M)
    Cov_perc_alpha_u_trt_orig2 <-  rep(NA, M)
    Cov_perc_alpha_l_trt_orig2 <-  rep(NA, M)
    Cov_alpha_CI_orig2 <- rep(NA, M)
    
    NA_alpha_u_avg <- rep(NA, M)
    NA_alpha_l_avg <- rep(NA, M)
    
    B_outer <- 400
    B_inner <- 50
    
    # Regarding these simulated data in log, the original optimization process cannot converge
    # Thus, only the outer bootstrap has the results, the inner bootstrap all failed.
    
    for(i in c(1:M)){
      
      print(i)
      
      Orig <- try(readRDS(file = paste0('orig_point_est/', 
                                        link_type, "_", sample_size,
                                        "/",  
                                        group, "_", IF_type, 
                                        "/select_", grp, "_", sample_size, "_", link_type, "_", IF_type, "_", i)), silent = T)
      Double <-  try(readRDS(file = paste0('double_boot_400_50/',
                                           link_type, "_", sample_size,
                                           "/", 
                                           group, "_", IF_type, 
                                           "/", 
                                           "double_", grp, "_", sample_size, "_", link_type, "_", IF_type,"_", i, "_corr")), silent = T)
      
      
      if(length(Double) > 1){
        
        alpha_u_boot_outer_orig <- rep(NA, B_outer)
        alpha_l_boot_outer_orig <- rep(NA, B_outer)
        alpha_u_boot_outer_corr <- rep(NA, B_outer)
        alpha_l_boot_outer_corr <- rep(NA, B_outer)
        NA_count_u_i <-  rep(NA, B_outer)
        NA_count_l_i <-  rep(NA, B_outer)
        
        for(j in 1:B_outer){
          
          temp <-  Double[[j]]
          alpha_l_boot_outer_orig[j] <- temp[[1]][1]
          alpha_u_boot_outer_orig[j] <- temp[[1]][2]
          alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[2]], na.rm = T)
          alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[3]], na.rm = T)
          
          NA_count_u_i[j] <- temp[[4]][1]
          NA_count_l_i[j] <- temp[[4]][2]
          
        }
        
        Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
        Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
        
        # Percentile CI
        CI_perc_upper_alpha_u_trt_corr1[i] <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
        CI_perc_low_alpha_l_trt_corr1[i] <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
        
        CI_perc_upper_alpha_u_trt_orig1[i] <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
        CI_perc_low_alpha_l_trt_orig1[i] <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
        
        Cov_perc_alpha_u_trt_corr1[i] <- ifelse(alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_corr1[i], 1, 0)
        Cov_perc_alpha_l_trt_corr1[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_corr1[i], 1, 0)
        Cov_perc_alpha_u_trt_orig1[i] <- ifelse(alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_orig1[i], 1, 0)
        Cov_perc_alpha_l_trt_orig1[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_orig1[i], 1, 0)
        
        Cov_alpha_CI_corr1[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_corr1[i] & alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_corr1[i], 1, 0)
        Cov_alpha_CI_orig1[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_orig1[i] & alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_orig1[i], 1, 0)
        
        # Wald CI
        CI_perc_low_alpha_l_trt_corr2[i] <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
        CI_perc_upper_alpha_u_trt_corr2[i] <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
        
        CI_perc_low_alpha_l_trt_orig2[i] <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
        CI_perc_upper_alpha_u_trt_orig2[i] <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
        
        Cov_perc_alpha_u_trt_corr2[i] <- ifelse(alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_corr2[i], 1, 0)
        Cov_perc_alpha_l_trt_corr2[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_corr2[i], 1, 0)
        Cov_perc_alpha_u_trt_orig2[i] <- ifelse(alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_orig2[i], 1, 0)
        Cov_perc_alpha_l_trt_orig2[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_orig2[i], 1, 0)
        
        Cov_alpha_CI_corr2[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_corr2[i] & alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_corr2[i], 1, 0)
        Cov_alpha_CI_orig2[i] <- ifelse(alpha_crl_l_true >= CI_perc_low_alpha_l_trt_orig2[i] & alpha_crl_u_true <= CI_perc_upper_alpha_u_trt_orig2[i], 1, 0)
        
        NA_alpha_u_avg[i] <- mean(NA_count_u_i, na.rm = T)
        NA_alpha_l_avg[i] <- mean(NA_count_l_i, na.rm = T)
        
      }
      
      
    }
    
    # check NA values
    which(is.na(Cov_alpha_CI_corr1) == T)
    which(is.na(Cov_alpha_CI_corr1) == T) |> length()
    
    
    # Percentile CI
    c(mean(Cov_perc_alpha_l_trt_orig1, na.rm = T), 
      mean(Cov_perc_alpha_u_trt_orig1, na.rm = T), 
      mean(Cov_alpha_CI_orig1, na.rm = T))
    c( mean(Cov_perc_alpha_l_trt_corr1, na.rm = T),
       mean(Cov_perc_alpha_u_trt_corr1, na.rm = T),
       mean(Cov_alpha_CI_corr1, na.rm = T))
    
    # Wald CI
    c(mean(Cov_perc_alpha_l_trt_orig2, na.rm = T),
      mean(Cov_perc_alpha_u_trt_orig2, na.rm = T),
      mean(Cov_alpha_CI_orig2, na.rm = T))
    c(mean(Cov_perc_alpha_l_trt_corr2, na.rm = T),
      mean(Cov_perc_alpha_u_trt_corr2, na.rm = T),
      mean(Cov_alpha_CI_corr2, na.rm = T))
    
    
    NA_alpha_u_avg
    NA_alpha_l_avg
    
  }
  
  
}


