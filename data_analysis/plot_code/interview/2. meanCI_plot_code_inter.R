
########## generate the bootstrap t confidence interval ##########
# load the bootstrap dataset
rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis")


######## Load the original simulated data ########
{
  
  ARC_data <- readRDS("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_data.rds")
  
  orig_simu_trt <- filter(ARC_data, Trt == "home_visits")
  orig_simu_crl <- filter(ARC_data, Trt == "control")
  
}


Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
spline_seq <- seq(60, 460, by = 1)
range_value <- c(1.2, 3)


# calculate bootstrap t confidence interval or JN confidence interval
Wald_CI <- function(alpha, refer_mean, var_paraboot, month){
  
  Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
  mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
  
  CI_left <- mean_orig - qnorm(0.975) * sqrt(var_paraboot)
  CI_right <- mean_orig + qnorm(0.975) * sqrt(var_paraboot)
  
  return(c(mean_orig,  CI_left, CI_right))
}

paraboot_var <- function(data, alpha, month){
  
  Month <- c(3, 6, 9, 12)
  Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
  index_month <- which(Month == month) + 1
  index_alpha <- which(Alpha_seq == alpha)
  
  Est_paraboot <- rep(NA, length(data))
  for(j in 1:length(data)){
    temp <- data[[j]][[1]]
    Est_paraboot[j] <- log(temp[index_alpha, index_month])
  }
  var_paraboot <- var(Est_paraboot)
  
  return(var_paraboot)
}


## log, weight1 
{
  
  # Load the paraboot estimation results 
  data_paraboot_trt <- readRDS("result/CI_plot_inter/paraboot_trt_log_weight1")
  data_paraboot_crl <- readRDS("result/CI_plot_inter/paraboot_crl_log_weight1")
  
  # load the estimation results from original ARC_data (no bootstrap)
  Orig_res_trt <- readRDS("result/CI_plot_inter/result_trt_log_weight1.RData")
  Orig_res_crl <- readRDS("result/CI_plot_inter/result_crl_log_weight1.RData")
  
  
  ##### Control #####
  {
    
    # get the estimation of E[Y(t)] of the original data for specific month
    SIM_control_03 <- Orig_res_crl[[1]][, 2]
    SIM_control_06 <- Orig_res_crl[[1]][, 3]
    SIM_control_09 <- Orig_res_crl[[1]][, 4]
    SIM_control_12 <- Orig_res_crl[[1]][, 5]
    
    
    ##### month = 3 #####
    month <- 3
    alpha <- -0.6
    refer_mean <- SIM_control_03 
    var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month)
    t_SIM_control_03  <- exp(Wald_CI(alpha, refer_mean = log(SIM_control_03), var_paraboot, month))
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month)
      temp <- exp(Wald_CI(alpha, refer_mean = log(SIM_control_03), var_paraboot, month))
      t_SIM_control_03  <- rbind(t_SIM_control_03, temp)
    }
    
    plot_SIM_control_df3 <- data.frame(Alpha = c("-0.6","-0.3","0","0.3","0.6"),
                                       Time = c(80, 85, 90, 95, 100),
                                       Estimate = t_SIM_control_03[, 1],
                                       Lower_ci = t_SIM_control_03[, 2],
                                       Upper_ci = t_SIM_control_03[, 3])
    
    
    ##### month = 6 #####
    month <- 6
    alpha <- -0.6
    refer_mean <- SIM_control_06 
    var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month)
    t_SIM_control_06  <- exp(Wald_CI(alpha, refer_mean = log(SIM_control_06), var_paraboot, month))
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month)
      temp <- exp(Wald_CI(alpha, refer_mean = log(SIM_control_06), var_paraboot, month))
      t_SIM_control_06  <- rbind(t_SIM_control_06, temp)
    }
    
    plot_SIM_control_df6 <- data.frame(Alpha = c("-0.6","-0.3","0","0.3","0.6"),
                                       Time = c(170, 175, 180, 185, 190),
                                       Estimate = t_SIM_control_06[, 1],
                                       Lower_ci = t_SIM_control_06[, 2],
                                       Upper_ci = t_SIM_control_06[, 3])
    
    
    ##### month = 9 #####
    month <- 9
    alpha <- -0.6
    refer_mean <- SIM_control_09 
    var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month)
    t_SIM_control_09  <- exp(Wald_CI(alpha, refer_mean = log(SIM_control_09), var_paraboot, month))
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month)
      temp <- exp(Wald_CI(alpha, refer_mean = log(SIM_control_09), var_paraboot, month))
      t_SIM_control_09  <- rbind(t_SIM_control_09, temp)
    }
    
    plot_SIM_control_df9 <- data.frame(Alpha = c("-0.6","-0.3","0","0.3","0.6"),
                                       Time = c(260, 265, 270, 275, 280),
                                       Estimate = t_SIM_control_09[, 1],
                                       Lower_ci = t_SIM_control_09[, 2],
                                       Upper_ci = t_SIM_control_09[, 3])
    
    ##### month = 12 #####
    month <- 12
    alpha <- -0.6
    refer_mean <- SIM_control_12
    var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month)
    t_SIM_control_12  <-  exp(Wald_CI(alpha, refer_mean = log(SIM_control_12), var_paraboot, month))
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month)
      temp <- exp(Wald_CI(alpha, refer_mean = log(SIM_control_12), var_paraboot, month))
      t_SIM_control_12  <- rbind(t_SIM_control_12, temp)
    }
    
    plot_SIM_control_df12 <- data.frame(Alpha = c("-0.6","-0.3","0","0.3","0.6"),
                                        Time = c(350, 355, 360, 365, 370),
                                        Estimate = t_SIM_control_12[, 1],
                                        Lower_ci = t_SIM_control_12[, 2],
                                        Upper_ci = t_SIM_control_12[, 3])
    
    
    plot_df_SIM_control <- rbind(plot_SIM_control_df3,
                                 plot_SIM_control_df6,
                                 plot_SIM_control_df9,
                                 plot_SIM_control_df12)
    plot_df_SIM_control
    
    
    CI_plot_crl <- ggplot(plot_df_SIM_control, aes(x = Time, y = Estimate, colour = Alpha))+
      geom_point(size = 3)+
      # geom_point(data = data_orig_crl, size = 0.5)+
      geom_errorbar(aes(ymin = Lower_ci,ymax = Upper_ci),size = 1)+
      theme_bw()+
      expand_limits(x=c(30,480),y=c(0,5))+
      scale_x_continuous(breaks=seq(90,450,by=90))+
      scale_colour_manual(values = c( # "#253494",
                                      "#41b6c4", 
                                     "#a1dab4", "#ffffb2",
                                     "#fd8d3c", "#bd0026"),
                          # "#f03b20",# "#fecc5c",# "#2c7fb8", # "grey50"),
                          breaks=c("0.6", "0.3", "0", "-0.3", "-0.6")) + 
      labs(x="Time (Days)", y="ACQ Score", color=expression(alpha)) +
      labs(color=expression(alpha))+
      theme(axis.title=element_text(size=25),
            legend.title=element_text(size=20),
            axis.text=element_text(size=18),
            legend.text=element_text(size=18))+
      annotate("text", x = 400, y = 4.5,label="PT Arm",size=6)+
      geom_hline(yintercept=range_value[1],lty="dotted",size=2)+
      geom_hline(yintercept=range_value[2],lty="dotted",size=2)
    
    CI_plot_crl
    
    ggsave(plot = CI_plot_crl, 
           width = 8, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots_inter/CI_plot/CI_crl_log_weight1.pdf")
    
  }
  
  
  ##### Treat #####
  {
    
    
    # get the estimation of E[Y(t)] of the original data for specific month
    SIM_treat_03 <- Orig_res_trt[[1]][, 2]
    SIM_treat_06 <- Orig_res_trt[[1]][, 3]
    SIM_treat_09 <- Orig_res_trt[[1]][, 4]
    SIM_treat_12 <- Orig_res_trt[[1]][, 5]
    
    
    ##### month = 3 #####
    month <- 3
    alpha <- -0.6
    refer_mean <- SIM_treat_03 
    var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month)
    t_SIM_treat_03  <- exp(Wald_CI(alpha, refer_mean = log(SIM_treat_03), var_paraboot, month))
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month)
      temp <- exp(Wald_CI(alpha, refer_mean = log(SIM_treat_03), var_paraboot, month))
      t_SIM_treat_03  <- rbind(t_SIM_treat_03, temp)
    }
    
    plot_SIM_treat_df3 <- data.frame(Alpha = c("-0.6","-0.3","0","0.3","0.6"),
                                       Time = c(80, 85, 90, 95, 100),
                                       Estimate = t_SIM_treat_03[, 1],
                                       Lower_ci = t_SIM_treat_03[, 2],
                                       Upper_ci = t_SIM_treat_03[, 3])
    
    
    ##### month = 6 #####
    month <- 6
    alpha <- -0.6
    refer_mean <- SIM_treat_06 
    var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month)
    t_SIM_treat_06  <- exp(Wald_CI(alpha, refer_mean = log(SIM_treat_06), var_paraboot, month))
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month)
      temp <- exp(Wald_CI(alpha, refer_mean = log(SIM_treat_06), var_paraboot, month))
      t_SIM_treat_06  <- rbind(t_SIM_treat_06, temp)
    }
    
    plot_SIM_treat_df6 <- data.frame(Alpha = c("-0.6","-0.3","0","0.3","0.6"),
                                       Time = c(170, 175, 180, 185, 190),
                                       Estimate = t_SIM_treat_06[, 1],
                                       Lower_ci = t_SIM_treat_06[, 2],
                                       Upper_ci = t_SIM_treat_06[, 3])
    
    
    ##### month = 9 #####
    month <- 9
    alpha <- -0.6
    refer_mean <- SIM_treat_09 
    var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month)
    t_SIM_treat_09  <- exp(Wald_CI(alpha, refer_mean = log(SIM_treat_09), var_paraboot, month))
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month)
      temp <- exp(Wald_CI(alpha, refer_mean = log(SIM_treat_09), var_paraboot, month))
      t_SIM_treat_09  <- rbind(t_SIM_treat_09, temp)
    }
    
    plot_SIM_treat_df9 <- data.frame(Alpha = c("-0.6","-0.3","0","0.3","0.6"),
                                       Time = c(260, 265, 270, 275, 280),
                                       Estimate = t_SIM_treat_09[, 1],
                                       Lower_ci = t_SIM_treat_09[, 2],
                                       Upper_ci = t_SIM_treat_09[, 3])
    
    ##### month = 12 #####
    month <- 12
    alpha <- -0.6
    refer_mean <- SIM_treat_12
    var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month)
    t_SIM_treat_12  <-  exp(Wald_CI(alpha, refer_mean = log(SIM_treat_12), var_paraboot, month))
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month)
      temp <- exp(Wald_CI(alpha, refer_mean = log(SIM_treat_12), var_paraboot, month))
      t_SIM_treat_12  <- rbind(t_SIM_treat_12, temp)
    }
    
    plot_SIM_treat_df12 <- data.frame(Alpha = c("-0.6","-0.3","0","0.3","0.6"),
                                        Time = c(350, 355, 360, 365, 370),
                                        Estimate = t_SIM_treat_12[, 1],
                                        Lower_ci = t_SIM_treat_12[, 2],
                                        Upper_ci = t_SIM_treat_12[, 3])
    
    
    plot_df_SIM_treat <- rbind(plot_SIM_treat_df3,
                                 plot_SIM_treat_df6,
                                 plot_SIM_treat_df9,
                                 plot_SIM_treat_df12)
    plot_df_SIM_treat
    
    
    CI_plot_trt <- ggplot(plot_df_SIM_treat, aes(x = Time, y = Estimate, colour = Alpha))+
      geom_point(size = 3)+
      # geom_point(data = data_orig_trt, size = 0.5)+
      geom_errorbar(aes(ymin = Lower_ci,ymax = Upper_ci),size = 1)+
      theme_bw()+
      expand_limits(x=c(30,480),y=c(0,5))+
      scale_x_continuous(breaks=seq(90,450,by=90))+
      scale_colour_manual(values = c( # "#253494",
                                      "#41b6c4", 
                                     "#a1dab4", "#ffffb2",
                                     "#fd8d3c", "#bd0026"),
                              # "#f03b20",# "#fecc5c",# "#2c7fb8", # "grey50"),
                          breaks=c("0.6", "0.3", "0", "-0.3", "-0.6")) + 
      # "original data"))+
      labs(x="Time (Days)", y="ACQ Score", color=expression(alpha)) +
      labs(color=expression(alpha))+
      theme(axis.title=element_text(size=25),
            legend.title=element_text(size=20),
            axis.text=element_text(size=18),
            legend.text=element_text(size=18))+
      theme(legend.position = "none")+
      annotate("text", x = 400, y = 4.5,label="PT+HV Arm",size=6)+
      geom_hline(yintercept=range_value[1],lty="dotted",size=2)+
      geom_hline(yintercept=range_value[2],lty="dotted",size=2)
    
    CI_plot_trt
    
    ggsave(plot = CI_plot_trt, 
           width = 7, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots_inter/CI_plot/CI_trt_log_weight1.pdf")
    
  }

}

