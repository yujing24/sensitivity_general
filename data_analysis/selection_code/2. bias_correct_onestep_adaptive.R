rm(list=ls())

library(dplyr)
library(survival)
library(MAVE)
library(dfoptim)
library(parallel)
library(orthogonalsplinebasis) 
library(BB)
library(statmod)  # For Gauss-Legendre quadrature



## One-step Multiplier Bias-correction on beta
{

  # first use the binary search method, then use the grid search method
  f_onestep_bias_adaptive <- function(seed,
                                      orig_input = orig_input,
                                      # Alpha_seq = seq(-1.5, 1.5, by = 0.01),
                                      Alpha_seq = seq(-3, 3, by = 0.01),
                                      weight_type = "rademacher", # rademacher, normal, exp
                                      link_type = "log", # "log", "logit", "identity"
                                      search_method = "adaptive", # "binary", "grid" or "adaptive"
                                      range = c(0.2, 0.5) # c(0.2, 0.5) or c(1.2, 3)
  ){
    
    library(dplyr)
    library(orthogonalsplinebasis) 
    
    Beta_all = orig_input[[2]]
    IF_mat_all = orig_input[[3]]
    IF_derive_inverse_all = orig_input[[4]]
    Convergence_all = orig_input[[5]]
    
    N <- dim(IF_mat_all[[1]])[2]
    
    
    ### Define functions
    {
      s_identity <- function(z) z
      
      s_logit <- function(x) 1 / (1 + exp(-x))
      
      s_log <- function(z) exp(z)
      
      if(link_type == "identity"){
        s_link = s_identity
      }else if(link_type == "log"){
        s_link = s_log
      }else{
        s_link = s_logit
      }
      
      # rademacher, normal, exponential
      if(weight_type == "rademacher"){
        ## Rademacher weights
        f_weight <- function(N){
          return(sample(c(-1,1), N, replace = T, prob = c(0.5,0.5))) 
        }
      }else if(weight_type == "normal"){
        ## standard normal weights
        f_weight <- function(N){
          return(rnorm(N, mean = 0, sd = 1))
        }
      }else{
        ## Exp(1) weights
        f_weight <- function(N){
          return(rexp(N, rate = 1) - 1)
        }
      }
      
      spline_seq_orig <- seq(60, 460, by=1)
      basis <- SplineBasis(knots = c(60,60,60,60,260,460,460,460,460))
      B_t_matrix_orig <- evaluate(basis, spline_seq_orig)
      
      f_onestep_beta_boot <- function(alpha_index){
        
        Beta_orig <- Beta_all[[alpha_index]]
        IF_mat <- IF_mat_all[[alpha_index]]
        IF_mat_deriv_inverse <- IF_derive_inverse_all[[alpha_index]]
        
        error <- rep(0, 5)
        for(j in 1:N){
          error <- error + IF_mat_deriv_inverse %*% ((W_all[j]) * IF_mat[, j]) / N
        }
        Beta_boot <-  Beta_orig - error
        E_Yt_boot <- s_link(B_t_matrix_orig %*% Beta_boot)
        Max_E_Yt_boot <- max(E_Yt_boot)
        Min_E_Yt_boot <- min(E_Yt_boot)
        
        return(c(Max_E_Yt_boot, Min_E_Yt_boot))
      }
      
    }
    
    set.seed(seed)
    W_all <- f_weight(N) # generate weights for all alpha
    # W_all <- (W_all - mean(W_all)) / sd(W_all)
    length_alpha <- length(Alpha_seq)
    
    if(search_method == "adaptive"){
      
      if(max(Convergence_all) > 0){
        
        # use the grid search method 
        
        res_outer <- lapply(1:length(Alpha_seq), f_onestep_beta_boot)
        Max_E_Yt_boot <-  sapply(res_outer, function(x) x[1])
        Min_E_Yt_boot <- sapply(res_outer, function(x) x[2])
        
        ########## Find alpha_u and alpha_l ##########
        res_u_boot <- Alpha_seq[max(which(Max_E_Yt_boot < range[2]))]
        res_l_boot <- Alpha_seq[min(which(Min_E_Yt_boot > range[1]))]
        
        return(list(c(res_l_boot, res_u_boot), 
                    "grid (cvg)"))
        
      }else{
        
        # use binary search method first 
        
        find_alpha_u <- function(Alpha_seq, threshold = range[2]){
          left <- 1
          right <- length_alpha
          
          val_left <- f_onestep_beta_boot(left)[1]
          val_right <- f_onestep_beta_boot(right)[1]
          
          if (val_left >= threshold) {
            return(NA)
          }
          
          if (val_right < threshold) {
            return(Alpha_seq[right])
          }
          
          while (right - left > 1) {
            mid <- floor((left + right) / 2)
            val_mid <- f_onestep_beta_boot(mid)[1]
            
            if (val_mid < threshold) {
              left <- mid
            } else {
              right <- mid
            }
          }
          
          Alpha_seq[left]
        }
        
        find_alpha_l <- function(Alpha_seq, threshold = range[1]){
          left <- 1
          right <- length_alpha
          
          val_left <- f_onestep_beta_boot(left)[2]
          val_right <- f_onestep_beta_boot(right)[2]
          
          if (val_right <= threshold) {
            return(NA)
          }
          
          if (val_left > threshold) {
            return(Alpha_seq[left])
          }
          
          while (right - left > 1) {
            mid <- floor((left + right) / 2)
            val_mid <- f_onestep_beta_boot(mid)[2]
            
            if (val_mid > threshold) {
              right <- mid
            } else {
              left <- mid
            }
          }
          
          Alpha_seq[right]
        }
        
        res_u_boot <- find_alpha_u(Alpha_seq = Alpha_seq, threshold = range[2])
        
        res_l_boot <- find_alpha_l(Alpha_seq = Alpha_seq, threshold = range[1])
        
        if((is.na(res_u_boot) == F) && (is.na(res_l_boot) == F)){
          
          return(list(c(res_l_boot, res_u_boot), 
                      "binary"))
          
        }else{
          
          # rerun this case using the grid search method
          
          res_outer <- lapply(1:length(Alpha_seq), f_onestep_beta_boot)
          Max_E_Yt_boot <-  sapply(res_outer, function(x) x[1])
          Min_E_Yt_boot <- sapply(res_outer, function(x) x[2])
          
          ########## Find alpha_u and alpha_l ##########
          res_u_boot <- Alpha_seq[max(which(Max_E_Yt_boot < range[2]))]
          res_l_boot <- Alpha_seq[min(which(Min_E_Yt_boot > range[1]))]
          
          return(list(c(res_l_boot, res_u_boot), 
                      "grid (NA)"))
          
        }
        
      }
      
      
    }
    
    if(search_method == "binary"){
      
      # use binary search method first 
      
      find_alpha_u <- function(Alpha_seq, threshold = range[2]){
        left <- 1
        right <- length_alpha
        
        val_left <- f_onestep_beta_boot(left)[1]
        val_right <- f_onestep_beta_boot(right)[1]
        
        if (val_left >= threshold) {
          return(NA)
        }
        
        if (val_right < threshold) {
          return(Alpha_seq[right])
        }
        
        while (right - left > 1) {
          mid <- floor((left + right) / 2)
          val_mid <- f_onestep_beta_boot(mid)[1]
          
          if (val_mid < threshold) {
            left <- mid
          } else {
            right <- mid
          }
        }
        
        Alpha_seq[left]
      }
      
      find_alpha_l <- function(Alpha_seq, threshold = range[1]){
        left <- 1
        right <- length_alpha
        
        val_left <- f_onestep_beta_boot(left)[2]
        val_right <- f_onestep_beta_boot(right)[2]
        
        if (val_right <= threshold) {
          return(NA)
        }
        
        if (val_left > threshold) {
          return(Alpha_seq[left])
        }
        
        while (right - left > 1) {
          mid <- floor((left + right) / 2)
          val_mid <- f_onestep_beta_boot(mid)[2]
          
          if (val_mid > threshold) {
            right <- mid
          } else {
            left <- mid
          }
        }
        
        Alpha_seq[right]
      }
      
      res_u_boot <- find_alpha_u(Alpha_seq = Alpha_seq, threshold = range[2])
      
      res_l_boot <- find_alpha_l(Alpha_seq = Alpha_seq, threshold = range[1])
      
      return(list(c(res_l_boot, res_u_boot), 
                  "binary"))
      
    }
    
    if(search_method == "grid"){
      
      # use the grid search method 
      
      res_outer <- lapply(1:length(Alpha_seq), f_onestep_beta_boot)
      Max_E_Yt_boot <-  sapply(res_outer, function(x) x[1])
      Min_E_Yt_boot <- sapply(res_outer, function(x) x[2])
      
      ########## Find alpha_u and alpha_l ##########
      res_u_boot <- Alpha_seq[max(which(Max_E_Yt_boot < range[2]))]
      res_l_boot <- Alpha_seq[min(which(Min_E_Yt_boot > range[1]))]
      
      return(list(c(res_l_boot, res_u_boot), 
                  "grid"))
      
    }
    
  }
  
}


Boot_num <- 200

setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm/select_alpha/onestep_adaptive")



# log link
{
  
  search_method = "grid"
  # search_method = "adaptive"
  
  # trt, weight1
  orig_trt_log_weight1 <- readRDS("orig_trt_log_weight1.RData")
  
  set.seed(123)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_trt_log_weight1,
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "log", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(1.2, 3),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_trt_log_weight1_", search_method, ".RData"))
  
  
  # trt, weight2
  orig_trt_log_weight2 <- readRDS("orig_trt_log_weight2.RData")
  
  set.seed(234)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_trt_log_weight2,
                                    # Alpha_seq = seq(-1.5, 1.5, by = 0.01),
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "log", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(1.2, 3),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_trt_log_weight2_", search_method, ".RData"))
  
  
  # crl, weight1
  orig_crl_log_weight1 <- readRDS("orig_crl_log_weight1.RData")
  
  set.seed(345)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_crl_log_weight1,
                                    # Alpha_seq = seq(-1.5, 1.5, by = 0.01),
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "log", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(1.2, 3),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_crl_log_weight1_", search_method, ".RData"))
  
  
  # crl, weight2
  orig_crl_log_weight2 <- readRDS("orig_crl_log_weight2.RData")
  
  set.seed(456)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_crl_log_weight2,
                                    # Alpha_seq = seq(-1.5, 1.5, by = 0.01),
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "log", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(1.2, 3),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_crl_log_weight2_", search_method, ".RData"))
  
  
}


# logit link 
{
  
  # search_method = "grid"
  search_method = "adaptive"
  
  # trt, weight1
  orig_trt_log_weight1 <- readRDS("orig_trt_logit_weight1.RData")
  
  set.seed(123)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_trt_log_weight1,
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "logit", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(0.2, 0.5),
                                    mc.cores = 64)
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_trt_logit_weight1_", search_method, ".RData"))
  
  
  # trt, weight2
  orig_trt_log_weight2 <- readRDS("orig_trt_logit_weight2.RData")
  
  set.seed(234)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_trt_log_weight2,
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "logit", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(0.2, 0.5),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_trt_logit_weight2_", search_method, ".RData"))
  
  
  # crl, weight1
  orig_crl_log_weight1 <- readRDS("orig_crl_logit_weight1.RData")
  
  set.seed(345)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_crl_log_weight1,
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "logit", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(0.2, 0.5),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_crl_logit_weight1_", search_method, ".RData"))
  
  
  # crl, weight2
  orig_crl_log_weight2 <- readRDS("orig_crl_logit_weight2.RData")
  
  set.seed(456)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_crl_log_weight2,
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "logit", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(0.2, 0.5),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_crl_logit_weight2_", search_method, ".RData"))
  
}


# identity link - not run
{
  
  # trt
  orig_trt_identity <- readRDS("orig_trt_identity.RData")
  
  set.seed(567)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_trt_identity,
                                    Alpha_seq = seq(-1.5, 1.5, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "identity", # "log", "logit", "identity"
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_trt_identity.RData"))
  
  
  # crl
  orig_crl_identity <- readRDS("orig_crl_identity.RData")
  
  set.seed(678)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_crl_identity,
                                    Alpha_seq = seq(-1.5, 1.5, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "identity", # "log", "logit", "identity"
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_crl_identity.RData"))
  
  
}





