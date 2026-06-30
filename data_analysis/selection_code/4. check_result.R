rm(list=ls())
library(dplyr)
library(orthogonalsplinebasis)
library(MASS)
library(ggplot2)
library(monotone) # unimonotone 



## onestep - point estimation
{
  
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm/select_alpha/onestep_adaptive")
  Boot_num <- 200
  
  # search_method = "grid"
  search_method = "adaptive"
  
  # log, weight1 
  {
      
    # trt
    orig_trt_log_weight1 <- readRDS("orig_trt_log_weight1.RData")
    
    Alpha_l_trt_point_orig <- orig_trt_log_weight1[[1]][1]
    Alpha_u_trt_point_orig <- orig_trt_log_weight1[[1]][2]
    
    corr_trt_log_weight1 <- readRDS(paste0('bias_correct_', Boot_num, '/corr_trt_log_weight1_', search_method, '.RData'))
    
    Boot_est_l <- rep(NA, Boot_num)
    Boot_est_u <- rep(NA, Boot_num)
    Type <- rep(NA, Boot_num)
    
    for(j in 1:Boot_num){
      
      temp <- corr_trt_log_weight1[[j]]
      Boot_est_l[j] <- temp[[1]][1]
      Boot_est_u[j] <- temp[[1]][2]
      Type[j] <- temp[[2]]
      
    }
    
    Alpha_l_trt_point_corr <- 2 * Alpha_l_trt_point_orig - mean(Boot_est_l, na.rm = T)
    Alpha_u_trt_point_corr <- 2 * Alpha_u_trt_point_orig - mean(Boot_est_u, na.rm = T)
    
    trt_log_weight1_orig <- c(Alpha_l_trt_point_orig, Alpha_u_trt_point_orig)
    trt_log_weight1_corr <- c(Alpha_l_trt_point_corr, Alpha_u_trt_point_corr)
    
    which(is.na(Boot_est_l) == T) |> length()
    which(is.na(Boot_est_u) == T) |> length()
    
   
    # crl
    orig_crl_log_weight1 <- readRDS("orig_crl_log_weight1.RData")
    
    Alpha_l_crl_point_orig <- orig_crl_log_weight1[[1]][1]
    Alpha_u_crl_point_orig <- orig_crl_log_weight1[[1]][2]
    
    corr_crl_log_weight1 <- readRDS(paste0('bias_correct_', Boot_num, '/corr_crl_log_weight1_', search_method, '.RData'))
    
    Boot_est_l <- rep(NA, Boot_num)
    Boot_est_u <- rep(NA, Boot_num)
    Type <- rep(NA, Boot_num)
    
    for(j in 1:Boot_num){
      
      temp <- corr_crl_log_weight1[[j]]
      Boot_est_l[j] <- temp[[1]][1]
      Boot_est_u[j] <- temp[[1]][2]
      Type[j] <- temp[[2]]
      
    }
    
    Alpha_l_crl_point_corr <- 2 * Alpha_l_crl_point_orig - mean(Boot_est_l, na.rm = T)
    Alpha_u_crl_point_corr <- 2 * Alpha_u_crl_point_orig - mean(Boot_est_u, na.rm = T)
    
    crl_log_weight1_orig <-  c(Alpha_l_crl_point_orig, Alpha_u_crl_point_orig)
    crl_log_weight1_corr <- c(Alpha_l_crl_point_corr, Alpha_u_crl_point_corr)
    
    which(is.na(Boot_est_l) == T) |> length()
    which(is.na(Boot_est_u) == T) |> length() 
    
  }
  
  # log, weight2
  {
    
    # trt
    orig_trt_log_weight2 <- readRDS("orig_trt_log_weight2.RData")
    
    Alpha_l_trt_point_orig <- orig_trt_log_weight2[[1]][1]
    Alpha_u_trt_point_orig <- orig_trt_log_weight2[[1]][2]
    
    corr_trt_log_weight2 <- readRDS(paste0('bias_correct_', Boot_num, '/corr_trt_log_weight2_', search_method, '.RData'))

    Boot_est_l <- rep(NA, Boot_num)
    Boot_est_u <- rep(NA, Boot_num)
    Type <- rep(NA, Boot_num)
    
    for(j in 1:Boot_num){
      
      temp <- corr_trt_log_weight2[[j]]
      Boot_est_l[j] <- temp[[1]][1]
      Boot_est_u[j] <- temp[[1]][2]
      Type[j] <- temp[[2]]
      
    }
    
    Alpha_l_trt_point_corr <- 2 * Alpha_l_trt_point_orig - mean(Boot_est_l, na.rm = T)
    Alpha_u_trt_point_corr <- 2 * Alpha_u_trt_point_orig - mean(Boot_est_u, na.rm = T)
    
    trt_log_weight2_orig <-  c(Alpha_l_trt_point_orig, Alpha_u_trt_point_orig)
    trt_log_weight2_corr <-  c(Alpha_l_trt_point_corr, Alpha_u_trt_point_corr)
    
    which(is.na(Boot_est_l) == T) |> length()
    which(is.na(Boot_est_u) == T) |> length()
    
    
    # crl
    orig_crl_log_weight2 <- readRDS("orig_crl_log_weight2.RData")
    
    Alpha_l_crl_point_orig <- orig_crl_log_weight2[[1]][1]
    Alpha_u_crl_point_orig <- orig_crl_log_weight2[[1]][2]
    
    corr_crl_log_weight2 <- readRDS(paste0('bias_correct_', Boot_num, '/corr_crl_log_weight2_', search_method, '.RData'))

    Boot_est_l <- rep(NA, Boot_num)
    Boot_est_u <- rep(NA, Boot_num)
    Type <- rep(NA, Boot_num)
    
    for(j in 1:Boot_num){
      
      temp <- corr_crl_log_weight2[[j]]
      Boot_est_l[j] <- temp[[1]][1]
      Boot_est_u[j] <- temp[[1]][2]
      Type[j] <- temp[[2]]
      
    }
    
    Alpha_l_crl_point_corr <- 2 * Alpha_l_crl_point_orig - mean(Boot_est_l, na.rm = T)
    Alpha_u_crl_point_corr <- 2 * Alpha_u_crl_point_orig - mean(Boot_est_u, na.rm = T)
    
    crl_log_weight2_orig <-   c(Alpha_l_crl_point_orig, Alpha_u_crl_point_orig)
    crl_log_weight2_corr <-   c(Alpha_l_crl_point_corr, Alpha_u_crl_point_corr)
    
    which(is.na(Boot_est_l) == T) |> length()
    which(is.na(Boot_est_u) == T) |> length() 
    
  }
  

  trt_log_weight1_orig
  trt_log_weight2_orig
  crl_log_weight1_orig
  crl_log_weight2_orig
  
  trt_log_weight1_corr
  trt_log_weight2_corr
  crl_log_weight1_corr
  crl_log_weight2_corr
  
  
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm/select_alpha/onestep_adaptive")
  Boot_num <- 200
  
  # search_method = "grid"
  search_method = "adaptive"
  
  
  # logit, weight1 
  {
    
    # trt
    orig_trt_logit_weight1 <- readRDS("orig_trt_logit_weight1.RData")
    
    Alpha_l_trt_point_orig <- orig_trt_logit_weight1[[1]][1]
    Alpha_u_trt_point_orig <- orig_trt_logit_weight1[[1]][2]
    
    corr_trt_logit_weight1 <- readRDS(paste0('bias_correct_', Boot_num, '/corr_trt_logit_weight1_', search_method, '.RData'))
    
    Boot_est_l <- rep(NA, Boot_num)
    Boot_est_u <- rep(NA, Boot_num)
    Type <- rep(NA, Boot_num)
    
    for(j in 1:Boot_num){
      
      temp <- corr_trt_logit_weight1[[j]]
      Boot_est_l[j] <- temp[[1]][1]
      Boot_est_u[j] <- temp[[1]][2]
      Type[j] <- temp[[2]]
      
    }
    
    Alpha_l_trt_point_corr <- 2 * Alpha_l_trt_point_orig - mean(Boot_est_l, na.rm = T)
    Alpha_u_trt_point_corr <- 2 * Alpha_u_trt_point_orig - mean(Boot_est_u, na.rm = T)
    
    trt_logit_weight1_orig <- c(Alpha_l_trt_point_orig, Alpha_u_trt_point_orig)
    trt_logit_weight1_corr <- c(Alpha_l_trt_point_corr, Alpha_u_trt_point_corr)
    
    which(is.na(Boot_est_l) == T) |> length()
    which(is.na(Boot_est_u) == T) |> length()
    
    
    # crl
    orig_crl_logit_weight1 <- readRDS("orig_crl_logit_weight1.RData")
    
    Alpha_l_crl_point_orig <- orig_crl_logit_weight1[[1]][1]
    Alpha_u_crl_point_orig <- orig_crl_logit_weight1[[1]][2]
    
    corr_crl_logit_weight1 <- readRDS(paste0('bias_correct_', Boot_num, '/corr_crl_logit_weight1_', search_method, '.RData'))
    
    Boot_est_l <- rep(NA, Boot_num)
    Boot_est_u <- rep(NA, Boot_num)
    Type <- rep(NA, Boot_num)
    
    for(j in 1:Boot_num){
      
      temp <- corr_crl_logit_weight1[[j]]
      Boot_est_l[j] <- temp[[1]][1]
      Boot_est_u[j] <- temp[[1]][2]
      Type[j] <- temp[[2]]
      
    }
    
    Alpha_l_crl_point_corr <- 2 * Alpha_l_crl_point_orig - mean(Boot_est_l, na.rm = T)
    Alpha_u_crl_point_corr <- 2 * Alpha_u_crl_point_orig - mean(Boot_est_u, na.rm = T)
    
    crl_logit_weight1_orig <-  c(Alpha_l_crl_point_orig, Alpha_u_crl_point_orig)
    crl_logit_weight1_corr <- c(Alpha_l_crl_point_corr, Alpha_u_crl_point_corr)
    
    which(is.na(Boot_est_l) == T) |> length()
    which(is.na(Boot_est_u) == T) |> length() 
    
  }
  
  
  # logit, weight2
  {
    
    # trt
    orig_trt_logit_weight2 <- readRDS("orig_trt_logit_weight2.RData")
    
    Alpha_l_trt_point_orig <- orig_trt_logit_weight2[[1]][1]
    Alpha_u_trt_point_orig <- orig_trt_logit_weight2[[1]][2]
    
    corr_trt_logit_weight2 <- readRDS(paste0('bias_correct_', Boot_num, '/corr_trt_logit_weight2_', search_method, '.RData'))
    
    Boot_est_l <- rep(NA, Boot_num)
    Boot_est_u <- rep(NA, Boot_num)
    Type <- rep(NA, Boot_num)
    
    for(j in 1:Boot_num){
      
      temp <- corr_trt_logit_weight2[[j]]
      Boot_est_l[j] <- temp[[1]][1]
      Boot_est_u[j] <- temp[[1]][2]
      Type[j] <- temp[[2]]
      
    }
    
    Alpha_l_trt_point_corr <- 2 * Alpha_l_trt_point_orig - mean(Boot_est_l, na.rm = T)
    Alpha_u_trt_point_corr <- 2 * Alpha_u_trt_point_orig - mean(Boot_est_u, na.rm = T)
    
    trt_logit_weight2_orig <- c(Alpha_l_trt_point_orig, Alpha_u_trt_point_orig)
    trt_logit_weight2_corr <- c(Alpha_l_trt_point_corr, Alpha_u_trt_point_corr)
    
    which(is.na(Boot_est_l) == T) |> length()
    which(is.na(Boot_est_u) == T) |> length()
    
    
    # crl
    orig_crl_logit_weight2 <- readRDS("orig_crl_logit_weight2.RData")
    
    Alpha_l_crl_point_orig <- orig_crl_logit_weight2[[1]][1]
    Alpha_u_crl_point_orig <- orig_crl_logit_weight2[[1]][2]
    
    corr_crl_logit_weight2 <- readRDS(paste0('bias_correct_', Boot_num, '/corr_crl_logit_weight2_', search_method, '.RData'))
    
    Boot_est_l <- rep(NA, Boot_num)
    Boot_est_u <- rep(NA, Boot_num)
    Type <- rep(NA, Boot_num)
    
    for(j in 1:Boot_num){
      
      temp <- corr_crl_logit_weight2[[j]]
      Boot_est_l[j] <- temp[[1]][1]
      Boot_est_u[j] <- temp[[1]][2]
      Type[j] <- temp[[2]]
      
    }
    
    Alpha_l_crl_point_corr <- 2 * Alpha_l_crl_point_orig - mean(Boot_est_l, na.rm = T)
    Alpha_u_crl_point_corr <- 2 * Alpha_u_crl_point_orig - mean(Boot_est_u, na.rm = T)
    
    crl_logit_weight2_orig <-  c(Alpha_l_crl_point_orig, Alpha_u_crl_point_orig)
    crl_logit_weight2_corr <- c(Alpha_l_crl_point_corr, Alpha_u_crl_point_corr)
    
    which(is.na(Boot_est_l) == T) |> length()
    which(is.na(Boot_est_u) == T) |> length() 
    
  }
  

  trt_logit_weight1_orig
  trt_logit_weight2_orig
  crl_logit_weight1_orig
  crl_logit_weight2_orig
  
  trt_logit_weight1_corr
  trt_logit_weight2_corr
  crl_logit_weight1_corr
  crl_logit_weight2_corr
  
  
  
  # identity - not run
  {
    
    # trt
    orig_trt_identity <- readRDS("orig_trt_identity.RData")
    
    Alpha_l_trt_point_orig <- orig_trt_identity[[1]][1]
    Alpha_u_trt_point_orig <- orig_trt_identity[[1]][2]
    
    corr_trt_identity <- readRDS(paste0('bias_correct_', Boot_num, '/corr_trt_identity.RData'))
    
    Boot_est_l <- rep(NA, Boot_num)
    Boot_est_u <- rep(NA, Boot_num)
    Type <- rep(NA, Boot_num)
    
    for(j in 1:Boot_num){
      
      temp <- corr_trt_identity[[j]]
      Boot_est_l[j] <- temp[[1]][1]
      Boot_est_u[j] <- temp[[1]][2]
      Type[j] <- temp[[2]]
      
    }
    
    Alpha_l_trt_point_corr <- 2 * Alpha_l_trt_point_orig - mean(Boot_est_l, na.rm = T)
    Alpha_u_trt_point_corr <- 2 * Alpha_u_trt_point_orig - mean(Boot_est_u, na.rm = T)
    
    trt_identity_orig <-  c(Alpha_l_trt_point_orig, Alpha_u_trt_point_orig)
    trt_identity_corr <-  c(Alpha_l_trt_point_corr, Alpha_u_trt_point_corr)
    
    which(is.na(Boot_est_l) == T) |> length()
    which(is.na(Boot_est_u) == T) |> length()
    
    
    # crl
    orig_crl_identity <- readRDS("orig_crl_identity.RData")
    
    Alpha_l_crl_point_orig <- orig_crl_identity[[1]][1]
    Alpha_u_crl_point_orig <- orig_crl_identity[[1]][2]
    
    corr_crl_identity <- readRDS(paste0('bias_correct_', Boot_num, '/corr_crl_identity.RData'))
    
    Boot_est_l <- rep(NA, Boot_num)
    Boot_est_u <- rep(NA, Boot_num)
    Type <- rep(NA, Boot_num)
    
    for(j in 1:Boot_num){
      
      temp <- corr_crl_identity[[j]]
      Boot_est_l[j] <- temp[[1]][1]
      Boot_est_u[j] <- temp[[1]][2]
      Type[j] <- temp[[2]]
      
    }
    
    Alpha_l_crl_point_corr <- 2 * Alpha_l_crl_point_orig - mean(Boot_est_l, na.rm = T)
    Alpha_u_crl_point_corr <- 2 * Alpha_u_crl_point_orig - mean(Boot_est_u, na.rm = T)
    
    crl_identity_orig <-   c(Alpha_l_crl_point_orig, Alpha_u_crl_point_orig)
    crl_identity_corr <-   c(Alpha_l_crl_point_corr, Alpha_u_crl_point_corr)
    
    which(is.na(Boot_est_l) == T) |> length()
    which(is.na(Boot_est_u) == T) |> length() 
    
  }
  
  trt_identity_orig
  trt_identity_corr
  crl_identity_orig
  crl_identity_corr
  
}


## onestep (binary search or adaptive search) - confidence set 
{
  
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm/select_alpha/onestep_adaptive")
  B_outer <- 400
  B_inner <- 50
  
  
  # log, treat
  {
    
    # trt
    # Orig <- readRDS("orig_trt_log_weight1.RData")
    # temp <-  readRDS("double_400_50/double_trt_log_weight1_adaptive.RData")
    # # temp <-  readRDS("double_400_50/double_trt_log_weight1_grid.RData")
  
    Orig <- readRDS("orig_trt_log_weight2.RData")
    temp <-  readRDS("double_400_50/double_trt_log_weight2_adaptive.RData")
    # temp <-  readRDS("double_400_50/double_trt_log_weight2_grid.RData")
    
    
    alpha_u_boot_outer_orig <- rep(NA, B_outer)
    alpha_l_boot_outer_orig <- rep(NA, B_outer)
    alpha_u_boot_outer_corr <- rep(NA, B_outer)
    alpha_l_boot_outer_corr <- rep(NA, B_outer)
    NA_count_u_i <-  rep(NA, B_outer)
    NA_count_l_i <-  rep(NA, B_outer)
    Type <-  rep(NA, B_outer)
    
    for(j in 1:B_outer){
      
      alpha_l_boot_outer_orig[j] <- temp[[j]][[1]][1]
      alpha_u_boot_outer_orig[j] <- temp[[j]][[1]][2]
      alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[j]][[2]], na.rm = T)
      alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[j]][[3]], na.rm = T)
      
      NA_count_u_i[j] <- temp[[j]][[4]][1]
      NA_count_l_i[j] <- temp[[j]][[4]][2]
      
      Type[j] <- temp[[j]][[5]]
      
    }
      
    Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
    Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
      
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_corr1 <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
    
    CI_perc_upper_alpha_u_trt_orig1 <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_orig1 <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T)
    
    # Wald CI
    CI_perc_low_alpha_l_trt_corr2 <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_corr2 <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
    
    CI_perc_low_alpha_l_trt_orig2 <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_orig2 <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
    
    
    # percentile CI
    c(CI_perc_low_alpha_l_trt_orig1, CI_perc_upper_alpha_u_trt_orig1)
    c(CI_perc_low_alpha_l_trt_corr1, CI_perc_upper_alpha_u_trt_corr1)
    
    # Wald CI
    c(CI_perc_low_alpha_l_trt_orig2, CI_perc_upper_alpha_u_trt_orig2)
    c(CI_perc_low_alpha_l_trt_corr2, CI_perc_upper_alpha_u_trt_corr2)
    
    
    max(NA_count_u_i)
    max(NA_count_l_i)
    
    
    # crl
    Orig <- readRDS("orig_crl_log_weight1.RData")
    temp <-  readRDS("double_400_50/double_crl_log_weight1_adaptive.RData")
    # temp <-  readRDS("double_400_50/double_crl_log_weight1_grid.RData")
    
    alpha_u_boot_outer_orig <- rep(NA, B_outer)
    alpha_l_boot_outer_orig <- rep(NA, B_outer)
    alpha_u_boot_outer_corr <- rep(NA, B_outer)
    alpha_l_boot_outer_corr <- rep(NA, B_outer)
    NA_count_u_i <-  rep(NA, B_outer)
    NA_count_l_i <-  rep(NA, B_outer)
    Type <-  rep(NA, B_outer)
    
    for(j in 1:B_outer){
      
      alpha_l_boot_outer_orig[j] <- temp[[j]][[1]][1]
      alpha_u_boot_outer_orig[j] <- temp[[j]][[1]][2]
      alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[j]][[2]], na.rm = T)
      alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[j]][[3]], na.rm = T)
      
      NA_count_u_i[j] <- temp[[j]][[4]][1]
      NA_count_l_i[j] <- temp[[j]][[4]][2]
      
      Type[j] <- temp[[j]][[5]]
    }
    
    Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
    Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_corr1 <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
    
    CI_perc_upper_alpha_u_trt_orig1 <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_orig1 <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
    
    # Wald CI
    CI_perc_low_alpha_l_trt_corr2 <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_corr2 <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
    
    CI_perc_low_alpha_l_trt_orig2 <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_orig2 <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
    
    c(CI_perc_low_alpha_l_trt_orig1, CI_perc_upper_alpha_u_trt_orig1)
    c(CI_perc_low_alpha_l_trt_corr1, CI_perc_upper_alpha_u_trt_corr1)
    
    max(NA_count_u_i)
    max(NA_count_l_i)
    
  }
    

  # log, control
  {
    
    # crl
    # Orig <- readRDS("orig_crl_log_weight1.RData")
    # temp <-  readRDS("double_400_50/double_crl_log_weight1_adaptive.RData")
    # # temp <-  readRDS("double_400_50/double_crl_log_weight1_grid.RData")
    
    Orig <- readRDS("orig_crl_log_weight2.RData")
    temp <-  readRDS("double_400_50/double_crl_log_weight2_adaptive.RData")
    # temp <-  readRDS("double_400_50/double_crl_log_weight2_grid.RData")
    
    alpha_u_boot_outer_orig <- rep(NA, B_outer)
    alpha_l_boot_outer_orig <- rep(NA, B_outer)
    alpha_u_boot_outer_corr <- rep(NA, B_outer)
    alpha_l_boot_outer_corr <- rep(NA, B_outer)
    NA_count_u_i <-  rep(NA, B_outer)
    NA_count_l_i <-  rep(NA, B_outer)
    Type <-  rep(NA, B_outer)
    
    for(j in 1:B_outer){
      
      alpha_l_boot_outer_orig[j] <- temp[[j]][[1]][1]
      alpha_u_boot_outer_orig[j] <- temp[[j]][[1]][2]
      alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[j]][[2]], na.rm = T)
      alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[j]][[3]], na.rm = T)
      
      NA_count_u_i[j] <- temp[[j]][[4]][1]
      NA_count_l_i[j] <- temp[[j]][[4]][2]
     
      Type[j] <- temp[[j]][[5]] 
    }
    
    Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
    Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_corr1 <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
    
    CI_perc_upper_alpha_u_trt_orig1 <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_orig1 <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
    
    # Wald CI
    CI_perc_low_alpha_l_trt_corr2 <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_corr2 <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
    
    CI_perc_low_alpha_l_trt_orig2 <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_orig2 <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
    
    # percentile CI
    c(CI_perc_low_alpha_l_trt_orig1, CI_perc_upper_alpha_u_trt_orig1)
    c(CI_perc_low_alpha_l_trt_corr1, CI_perc_upper_alpha_u_trt_corr1)
    
    # Wald CI
    c(CI_perc_low_alpha_l_trt_orig2, CI_perc_upper_alpha_u_trt_orig2)
    c(CI_perc_low_alpha_l_trt_corr2, CI_perc_upper_alpha_u_trt_corr2)
    
    max(NA_count_u_i)
    max(NA_count_l_i)
    
    
    
    # crl
    Orig <- readRDS("orig_crl_log_weight2.RData")
    temp <-  readRDS("double_400_50/double_crl_log_weight2_adaptive.RData")
    # temp <-  readRDS("double_400_50/double_crl_log_weight2_grid.RData")
    
    alpha_u_boot_outer_orig <- rep(NA, B_outer)
    alpha_l_boot_outer_orig <- rep(NA, B_outer)
    alpha_u_boot_outer_corr <- rep(NA, B_outer)
    alpha_l_boot_outer_corr <- rep(NA, B_outer)
    NA_count_u_i <-  rep(NA, B_outer)
    NA_count_l_i <-  rep(NA, B_outer)
    Type <-  rep(NA, B_outer)
    
    for(j in 1:B_outer){
      
      alpha_l_boot_outer_orig[j] <- temp[[j]][[1]][1]
      alpha_u_boot_outer_orig[j] <- temp[[j]][[1]][2]
      alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[j]][[2]], na.rm = T)
      alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[j]][[3]], na.rm = T)
      
      NA_count_u_i[j] <- temp[[j]][[4]][1]
      NA_count_l_i[j] <- temp[[j]][[4]][2]
      
      Type[j] <- temp[[j]][[5]]
    }
    
    Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
    Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_corr1 <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
    
    CI_perc_upper_alpha_u_trt_orig1 <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_orig1 <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
    
    # Wald CI
    CI_perc_low_alpha_l_trt_corr2 <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_corr2 <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
    
    CI_perc_low_alpha_l_trt_orig2 <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_orig2 <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
    
    c(CI_perc_low_alpha_l_trt_orig1, CI_perc_upper_alpha_u_trt_orig1)
    c(CI_perc_low_alpha_l_trt_corr1, CI_perc_upper_alpha_u_trt_corr1)
    
    max(NA_count_u_i)
    max(NA_count_l_i)
    
    
  }
  
  
  # identity - not run
  {
    
    # trt
    Orig <- readRDS("orig_trt_identity.RData")
    temp <-  readRDS("double_400_50/double_trt_identity.RData")
    
    alpha_u_boot_outer_orig <- rep(NA, B_outer)
    alpha_l_boot_outer_orig <- rep(NA, B_outer)
    alpha_u_boot_outer_corr <- rep(NA, B_outer)
    alpha_l_boot_outer_corr <- rep(NA, B_outer)
    NA_count_u_i <-  rep(NA, B_outer)
    NA_count_l_i <-  rep(NA, B_outer)
    Type <-  rep(NA, B_outer)
    
    for(j in 1:B_outer){
      
      alpha_l_boot_outer_orig[j] <- temp[[j]][[1]][1]
      alpha_u_boot_outer_orig[j] <- temp[[j]][[1]][2]
      alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[j]][[2]], na.rm = T)
      alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[j]][[3]], na.rm = T)
      
      NA_count_u_i[j] <- temp[[j]][[4]][1]
      NA_count_l_i[j] <- temp[[j]][[4]][2]
      
      Type[j] <- temp[[j]][[5]]
    }
    
    Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
    Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_corr1 <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
    
    CI_perc_upper_alpha_u_trt_orig1 <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_orig1 <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
    
    # Wald CI
    CI_perc_low_alpha_l_trt_corr2 <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_corr2 <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
    
    CI_perc_low_alpha_l_trt_orig2 <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_orig2 <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
    
    c(CI_perc_low_alpha_l_trt_orig1, CI_perc_upper_alpha_u_trt_orig1)
    c(CI_perc_low_alpha_l_trt_corr1, CI_perc_upper_alpha_u_trt_corr1)
    
    max(NA_count_u_i)
    max(NA_count_l_i)
    
    
    
    
    # crl
    Orig <- readRDS("orig_crl_identity.RData")
    temp <-  readRDS("double_400_50/double_crl_identity.RData")
    
    alpha_u_boot_outer_orig <- rep(NA, B_outer)
    alpha_l_boot_outer_orig <- rep(NA, B_outer)
    alpha_u_boot_outer_corr <- rep(NA, B_outer)
    alpha_l_boot_outer_corr <- rep(NA, B_outer)
    NA_count_u_i <-  rep(NA, B_outer)
    NA_count_l_i <-  rep(NA, B_outer)
    Type <-  rep(NA, B_outer)
    
    for(j in 1:B_outer){
      
      alpha_l_boot_outer_orig[j] <- temp[[j]][[1]][1]
      alpha_u_boot_outer_orig[j] <- temp[[j]][[1]][2]
      alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[j]][[2]], na.rm = T)
      alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[j]][[3]], na.rm = T)
      
      NA_count_u_i[j] <- temp[[j]][[4]][1]
      NA_count_l_i[j] <- temp[[j]][[4]][2]
      
      Type[j] <- temp[[j]][[5]]
      
    }
    
    Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
    Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_corr1 <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
    
    CI_perc_upper_alpha_u_trt_orig1 <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_orig1 <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
    
    # Wald CI
    CI_perc_low_alpha_l_trt_corr2 <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_corr2 <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
    
    CI_perc_low_alpha_l_trt_orig2 <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_orig2 <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
    
    c(CI_perc_low_alpha_l_trt_orig1, CI_perc_upper_alpha_u_trt_orig1)
    c(CI_perc_low_alpha_l_trt_corr1, CI_perc_upper_alpha_u_trt_corr1)
    
    max(NA_count_u_i)
    max(NA_count_l_i)
    
    
  }
  
  
  # logit, trt
  {
    
    # trt
    # Orig <- readRDS("orig_trt_logit_weight1.RData")
    # temp <-  readRDS("double_400_50/double_trt_logit_weight1_adaptive.RData")
    # # temp <-  readRDS("double_400_50/double_trt_logit_weight1_grid.RData")
    
    Orig <- readRDS("orig_trt_logit_weight2.RData")
    temp <-  readRDS("double_400_50/double_trt_logit_weight2_adaptive.RData")
    # temp <-  readRDS("double_400_50/double_trt_logit_weight2_grid.RData")
    
    alpha_u_boot_outer_orig <- rep(NA, B_outer)
    alpha_l_boot_outer_orig <- rep(NA, B_outer)
    alpha_u_boot_outer_corr <- rep(NA, B_outer)
    alpha_l_boot_outer_corr <- rep(NA, B_outer)
    NA_count_u_i <-  rep(NA, B_outer)
    NA_count_l_i <-  rep(NA, B_outer)
    Type <-  rep(NA, B_outer)
    
    for(j in 1:B_outer){
      
      alpha_l_boot_outer_orig[j] <- temp[[j]][[1]][1]
      alpha_u_boot_outer_orig[j] <- temp[[j]][[1]][2]
      alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[j]][[2]], na.rm = T)
      alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[j]][[3]], na.rm = T)
      
      NA_count_u_i[j] <- temp[[j]][[4]][1]
      NA_count_l_i[j] <- temp[[j]][[4]][2]
      
      Type[j] <- temp[[j]][[5]]
      
    }
    
    Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
    Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_corr1 <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
    
    CI_perc_upper_alpha_u_trt_orig1 <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_orig1 <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
    
    # Wald CI
    CI_perc_low_alpha_l_trt_corr2 <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_corr2 <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
    
    CI_perc_low_alpha_l_trt_orig2 <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_orig2 <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
    
    # percentile CI
    c(CI_perc_low_alpha_l_trt_orig1, CI_perc_upper_alpha_u_trt_orig1)
    c(CI_perc_low_alpha_l_trt_corr1, CI_perc_upper_alpha_u_trt_corr1)
    
    # Wald CI
    c(CI_perc_low_alpha_l_trt_orig2, CI_perc_upper_alpha_u_trt_orig2)
    c(CI_perc_low_alpha_l_trt_corr2, CI_perc_upper_alpha_u_trt_corr2)
    
    
    max(NA_count_u_i)
    max(NA_count_l_i)
    
    
    # crl
    Orig <- readRDS("orig_crl_logit_weight1.RData")
    temp <-  readRDS("double_400_50/double_crl_logit_weight1_adaptive.RData")
    # temp <-  readRDS("double_400_50/double_crl_logit_weight1_grid.RData")
    
    alpha_u_boot_outer_orig <- rep(NA, B_outer)
    alpha_l_boot_outer_orig <- rep(NA, B_outer)
    alpha_u_boot_outer_corr <- rep(NA, B_outer)
    alpha_l_boot_outer_corr <- rep(NA, B_outer)
    NA_count_u_i <-  rep(NA, B_outer)
    NA_count_l_i <-  rep(NA, B_outer)
    Type <-  rep(NA, B_outer)
    
    for(j in 1:B_outer){
      
      alpha_l_boot_outer_orig[j] <- temp[[j]][[1]][1]
      alpha_u_boot_outer_orig[j] <- temp[[j]][[1]][2]
      alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[j]][[2]], na.rm = T)
      alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[j]][[3]], na.rm = T)
      
      NA_count_u_i[j] <- temp[[j]][[4]][1]
      NA_count_l_i[j] <- temp[[j]][[4]][2]
      
      Type[j] <- temp[[j]][[5]]
    }
    
    Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
    Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_corr1 <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
    
    CI_perc_upper_alpha_u_trt_orig1 <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_orig1 <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
    
    # Wald CI
    CI_perc_low_alpha_l_trt_corr2 <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_corr2 <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
    
    CI_perc_low_alpha_l_trt_orig2 <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_orig2 <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
    
    # percentile CI
    c(CI_perc_low_alpha_l_trt_orig1, CI_perc_upper_alpha_u_trt_orig1)
    c(CI_perc_low_alpha_l_trt_corr1, CI_perc_upper_alpha_u_trt_corr1)
    
    # Wald CI
    c(CI_perc_low_alpha_l_trt_orig2, CI_perc_upper_alpha_u_trt_orig2)
    c(CI_perc_low_alpha_l_trt_corr2, CI_perc_upper_alpha_u_trt_corr2)
    
    max(NA_count_u_i)
    max(NA_count_l_i)
    
  }
  
  # logit, crl
  {
    
    # trt
    # Orig <- readRDS("orig_crl_logit_weight1.RData")
    # temp <-  readRDS("double_400_50/double_crl_logit_weight1_adaptive.RData")
    # # temp <-  readRDS("double_400_50/double_crl_logit_weight1_grid.RData")

    Orig <- readRDS("orig_crl_logit_weight2.RData")
    temp <-  readRDS("double_400_50/double_crl_logit_weight2_adaptive.RData")
    # temp <-  readRDS("double_400_50/double_crl_logit_weight2_grid.RData")
    
    alpha_u_boot_outer_orig <- rep(NA, B_outer)
    alpha_l_boot_outer_orig <- rep(NA, B_outer)
    alpha_u_boot_outer_corr <- rep(NA, B_outer)
    alpha_l_boot_outer_corr <- rep(NA, B_outer)
    NA_count_u_i <-  rep(NA, B_outer)
    NA_count_l_i <-  rep(NA, B_outer)
    Type <-  rep(NA, B_outer)
    
    for(j in 1:B_outer){
      
      alpha_l_boot_outer_orig[j] <- temp[[j]][[1]][1]
      alpha_u_boot_outer_orig[j] <- temp[[j]][[1]][2]
      alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[j]][[2]], na.rm = T)
      alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[j]][[3]], na.rm = T)
      
      NA_count_u_i[j] <- temp[[j]][[4]][1]
      NA_count_l_i[j] <- temp[[j]][[4]][2]
      
      Type[j] <- temp[[j]][[5]] 
    }
    
    Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
    Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_corr1 <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
    
    CI_perc_upper_alpha_u_trt_orig1 <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_orig1 <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
    
    # Wald CI
    CI_perc_low_alpha_l_trt_corr2 <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_corr2 <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
    
    CI_perc_low_alpha_l_trt_orig2 <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_orig2 <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
    
    # percentile CI
    c(CI_perc_low_alpha_l_trt_orig1, CI_perc_upper_alpha_u_trt_orig1)
    c(CI_perc_low_alpha_l_trt_corr1, CI_perc_upper_alpha_u_trt_corr1)
    
    # Wald CI
    c(CI_perc_low_alpha_l_trt_orig2, CI_perc_upper_alpha_u_trt_orig2)
    c(CI_perc_low_alpha_l_trt_corr2, CI_perc_upper_alpha_u_trt_corr2)
    
    max(NA_count_u_i)
    max(NA_count_l_i)
    
    
    
    # crl
    Orig <- readRDS("orig_crl_logit_weight2.RData")
    temp <-  readRDS("double_400_50/double_crl_logit_weight2_adaptive.RData")
    # temp <-  readRDS("double_400_50/double_crl_logit_weight2_grid.RData")
    
    alpha_u_boot_outer_orig <- rep(NA, B_outer)
    alpha_l_boot_outer_orig <- rep(NA, B_outer)
    alpha_u_boot_outer_corr <- rep(NA, B_outer)
    alpha_l_boot_outer_corr <- rep(NA, B_outer)
    NA_count_u_i <-  rep(NA, B_outer)
    NA_count_l_i <-  rep(NA, B_outer)
    Type <-  rep(NA, B_outer)
    
    for(j in 1:B_outer){
      
      alpha_l_boot_outer_orig[j] <- temp[[j]][[1]][1]
      alpha_u_boot_outer_orig[j] <- temp[[j]][[1]][2]
      alpha_l_boot_outer_corr[j] <- 2 * alpha_l_boot_outer_orig[j] - mean(temp[[j]][[2]], na.rm = T)
      alpha_u_boot_outer_corr[j] <- 2 * alpha_u_boot_outer_orig[j] - mean(temp[[j]][[3]], na.rm = T)
      
      NA_count_u_i[j] <- temp[[j]][[4]][1]
      NA_count_l_i[j] <- temp[[j]][[4]][2]
      
      Type[j] <- temp[[j]][[5]]
    }
    
    Orig_corr_l <- 2 * Orig[[1]][1] - mean(alpha_l_boot_outer_orig, na.rm = T)
    Orig_corr_u <- 2 * Orig[[1]][2] - mean(alpha_u_boot_outer_orig, na.rm = T)
    
    # Percentile CI
    CI_perc_upper_alpha_u_trt_corr1 <- quantile(alpha_u_boot_outer_corr, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_corr1 <- quantile(alpha_l_boot_outer_corr, 0.025, na.rm = T)
    
    CI_perc_upper_alpha_u_trt_orig1 <- quantile(alpha_u_boot_outer_orig, 0.975, na.rm = T)
    CI_perc_low_alpha_l_trt_orig1 <- quantile(alpha_l_boot_outer_orig, 0.025, na.rm = T) 
    
    # Wald CI
    CI_perc_low_alpha_l_trt_corr2 <- Orig_corr_l - qnorm(0.975) * sd(alpha_l_boot_outer_corr, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_corr2 <- Orig_corr_u + qnorm(0.975) * sd(alpha_u_boot_outer_corr, na.rm = T) # extra
    
    CI_perc_low_alpha_l_trt_orig2 <- Orig[[1]][1]  - qnorm(0.975) * sd(alpha_l_boot_outer_orig, na.rm = T) # extra
    CI_perc_upper_alpha_u_trt_orig2 <- Orig[[1]][2]  + qnorm(0.975) * sd(alpha_u_boot_outer_orig, na.rm = T) # extra
    
    c(CI_perc_low_alpha_l_trt_orig1, CI_perc_upper_alpha_u_trt_orig1)
    c(CI_perc_low_alpha_l_trt_corr1, CI_perc_upper_alpha_u_trt_corr1)
    
    max(NA_count_u_i)
    max(NA_count_l_i)
    
    
  }
  
}
  
  

