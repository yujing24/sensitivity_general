
########## generate the bootstrap t confidence interval ##########
# load the bootstrap dataset
rm(list=ls())
library(ggplot2)
library(dplyr)

source("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm/1. function_inference.R")
setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")

Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2)
spline_seq <- seq(60, 460, by = 1)
range_value <- c(0.2, 0.5)


## logit, weight1 
{
  
  # Load the paraboot estimation results 
  data_paraboot_trt <- readRDS("result/CI_plot/paraboot_trt_logit_weight1")
  data_paraboot_crl <- readRDS("result/CI_plot/paraboot_crl_logit_weight1")
  
  # load the estimation results from original ARC_data (no bootstrap)
  Orig_res_trt <- readRDS("result/CI_plot/result_trt_logit_weight1.RData")
  Orig_res_crl <- readRDS("result/CI_plot/result_crl_logit_weight1.RData")
  
  
  ##### Control #####
  {

    # get the estimation of E[Y(t)] of the original data for specific month
    SIM_control_03 <- Orig_res_crl[[1]][, 2]
    SIM_control_06 <- Orig_res_crl[[1]][, 3]
    SIM_control_09 <- Orig_res_crl[[1]][, 4]
    SIM_control_12 <- Orig_res_crl[[1]][, 5]
    
    ##### month = 3 #####
    month <- 3
    alpha <- -0.8
    refer_mean <- SIM_control_03 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_control_03  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_control_03  <- rbind(t_SIM_control_03, temp)
    }
    
    plot_SIM_control_df3 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                       Time = c(75, 80, 85, 90, 95, 100),
                                       Estimate = t_SIM_control_03[, 1],
                                       Lower_ci = t_SIM_control_03[, 2],
                                       Upper_ci = t_SIM_control_03[, 3])
    
    
    ##### month = 6 #####
    month <- 6
    alpha <- -0.8
    refer_mean <- SIM_control_06 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_control_06  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_control_06  <- rbind(t_SIM_control_06, temp)
    }
    
    plot_SIM_control_df6 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                       Time = c(165, 170, 175, 180, 185, 190),
                                       Estimate = t_SIM_control_06[, 1],
                                       Lower_ci = t_SIM_control_06[, 2],
                                       Upper_ci = t_SIM_control_06[, 3])
    
    
    ##### month = 9 #####
    month <- 9
    alpha <- -0.8
    refer_mean <- SIM_control_09 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_control_09  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_control_09  <- rbind(t_SIM_control_09, temp)
    }
    
    plot_SIM_control_df9 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                       Time = c(255, 260, 265, 270, 275, 280),
                                       Estimate = t_SIM_control_09[, 1],
                                       Lower_ci = t_SIM_control_09[, 2],
                                       Upper_ci = t_SIM_control_09[, 3])
    
    ##### month = 12 #####
    month <- 12
    alpha <- -0.8
    refer_mean <- SIM_control_12 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_control_12  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_control_12  <- rbind(t_SIM_control_12, temp)
    }
    
    plot_SIM_control_df12 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                        Time = c(345, 350, 355, 360, 365, 370),
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
      expand_limits(x=c(30,480),y=c(0,1))+
      scale_x_continuous(breaks=seq(90,450,by=90))+
      scale_colour_manual(values = c("#253494", "#41b6c4", 
                                     "#a1dab4", "#ffffb2",
                                     "#fd8d3c", "#bd0026"),
                          # "#f03b20",# "#fecc5c",# "#2c7fb8", # "grey50"),
                          breaks=c("1.2", "0.8", "0.4", "0", "-0.4", "-0.8")) + 
      labs(x="Time (Days)", y="Rate of ACQ Score < 1.5", color=expression(alpha)) +
      labs(color=expression(alpha))+
      theme(axis.title=element_text(size=25),
            legend.title=element_text(size=20),
            axis.text=element_text(size=18),
            legend.text=element_text(size=18))+
      annotate("text", x = 400, y = 0.85,label="PT Arm",size=6)+
      geom_hline(yintercept=range_value[1],lty="dotted",size=2)+
      geom_hline(yintercept=range_value[2],lty="dotted",size=2)
    
    CI_plot_crl
    
    ggsave(plot = CI_plot_crl, 
           width = 8, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots/CI_plot/CI_crl_logit_weight1.pdf")
    
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
    alpha <- -0.8
    refer_mean <- SIM_treat_03 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_treat_03  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_treat_03  <- rbind(t_SIM_treat_03, temp)
    }
    
    plot_SIM_treat_df3 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                       Time = c(75, 80, 85, 90, 95, 100),
                                       Estimate = t_SIM_treat_03[, 1],
                                       Lower_ci = t_SIM_treat_03[, 2],
                                       Upper_ci = t_SIM_treat_03[, 3])
    
    
    ##### month = 6 #####
    month <- 6
    alpha <- -0.8
    refer_mean <- SIM_treat_06 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_treat_06  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_treat_06  <- rbind(t_SIM_treat_06, temp)
    }
    
    plot_SIM_treat_df6 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                       Time = c(165, 170, 175, 180, 185, 190),
                                       Estimate = t_SIM_treat_06[, 1],
                                       Lower_ci = t_SIM_treat_06[, 2],
                                       Upper_ci = t_SIM_treat_06[, 3])
    
    
    ##### month = 9 #####
    month <- 9
    alpha <- -0.8
    refer_mean <- SIM_treat_09 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_treat_09  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_treat_09  <- rbind(t_SIM_treat_09, temp)
    }
    
    plot_SIM_treat_df9 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                       Time = c(255, 260, 265, 270, 275, 280),
                                       Estimate = t_SIM_treat_09[, 1],
                                       Lower_ci = t_SIM_treat_09[, 2],
                                       Upper_ci = t_SIM_treat_09[, 3])
    
    ##### month = 12 #####
    month <- 12
    alpha <- -0.8
    refer_mean <- SIM_treat_12 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_treat_12  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_treat_12  <- rbind(t_SIM_treat_12, temp)
    }
    
    plot_SIM_treat_df12 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                        Time = c(345, 350, 355, 360, 365, 370),
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
      expand_limits(x=c(30,480),y=c(0,1))+
      scale_x_continuous(breaks=seq(90,450,by=90))+
      scale_colour_manual(values = c("#253494", "#41b6c4", 
                                     "#a1dab4", "#ffffb2",
                                     "#fd8d3c", "#bd0026"),
                              # "#f03b20",# "#fecc5c",# "#2c7fb8", # "grey50"),
                          breaks=c("1.2", "0.8", "0.4", "0", "-0.4", "-0.8")) + 
      # "original data"))+
      labs(x="Time (Days)", y="Rate of ACQ Score < 1.5", color=expression(alpha)) +
      labs(color=expression(alpha))+
      theme(axis.title=element_text(size=25),
            legend.title=element_text(size=20),
            axis.text=element_text(size=18),
            legend.text=element_text(size=18))+
      theme(legend.position = "none")+
      annotate("text", x = 400, y = 0.85,label="PT+HV Arm",size=6)+
      geom_hline(yintercept=range_value[1],lty="dotted",size=2)+
      geom_hline(yintercept=range_value[2],lty="dotted",size=2)
    
    CI_plot_trt
    
    ggsave(plot = CI_plot_trt, 
           width = 7, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots/CI_plot/CI_trt_logit_weight1.pdf")
    
  }

  
}


## logit, weight2
{
  
  # Load the paraboot estimation results 
  data_paraboot_trt <- readRDS("result/CI_plot/paraboot_trt_logit_weight2")
  data_paraboot_crl <- readRDS("result/CI_plot/paraboot_crl_logit_weight2")
  
  # load the estimation results from original ARC_data (no bootstrap)
  Orig_res_trt <- readRDS("result/CI_plot/result_trt_logit_weight2.RData")
  Orig_res_crl <- readRDS("result/CI_plot/result_crl_logit_weight2.RData")
  
  
  ##### Control #####
  {
    
    # get the estimation of E[Y(t)] of the original data for specific month
    SIM_control_03 <- Orig_res_crl[[1]][, 2]
    SIM_control_06 <- Orig_res_crl[[1]][, 3]
    SIM_control_09 <- Orig_res_crl[[1]][, 4]
    SIM_control_12 <- Orig_res_crl[[1]][, 5]
    
    
    ##### month = 3 #####
    month <- 3
    alpha <- -0.8
    refer_mean <- SIM_control_03 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_control_03  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_control_03  <- rbind(t_SIM_control_03, temp)
    }
    
    plot_SIM_control_df3 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                       Time = c(75, 80, 85, 90, 95, 100),
                                       Estimate = t_SIM_control_03[, 1],
                                       Lower_ci = t_SIM_control_03[, 2],
                                       Upper_ci = t_SIM_control_03[, 3])
    
    
    ##### month = 6 #####
    month <- 6
    alpha <- -0.8
    refer_mean <- SIM_control_06 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_control_06  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_control_06  <- rbind(t_SIM_control_06, temp)
    }
    
    plot_SIM_control_df6 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                       Time = c(165, 170, 175, 180, 185, 190),
                                       Estimate = t_SIM_control_06[, 1],
                                       Lower_ci = t_SIM_control_06[, 2],
                                       Upper_ci = t_SIM_control_06[, 3])
    
    
    ##### month = 9 #####
    month <- 9
    alpha <- -0.8
    refer_mean <- SIM_control_09 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_control_09  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_control_09  <- rbind(t_SIM_control_09, temp)
    }
    
    plot_SIM_control_df9 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                       Time = c(255, 260, 265, 270, 275, 280),
                                       Estimate = t_SIM_control_09[, 1],
                                       Lower_ci = t_SIM_control_09[, 2],
                                       Upper_ci = t_SIM_control_09[, 3])
    
    ##### month = 12 #####
    month <- 12
    alpha <- -0.8
    refer_mean <- SIM_control_12 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_control_12  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_crl, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_control_12  <- rbind(t_SIM_control_12, temp)
    }
    
    plot_SIM_control_df12 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                        Time = c(345, 350, 355, 360, 365, 370),
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
      expand_limits(x=c(30,480),y=c(0,1))+
      scale_x_continuous(breaks=seq(90,450,by=90))+
      scale_colour_manual(values = c("#253494", "#41b6c4", 
                                     "#a1dab4", "#ffffb2",
                                     "#fd8d3c", "#bd0026"),
                          # "#f03b20",# "#fecc5c",# "#2c7fb8", # "grey50"),
                          breaks=c("1.2", "0.8", "0.4", "0", "-0.4", "-0.8")) + 
      labs(x="Time (Days)", y="Rate of ACQ Score < 1.5", color=expression(alpha)) +
      labs(color=expression(alpha))+
      theme(axis.title=element_text(size=25),
            legend.title=element_text(size=20),
            axis.text=element_text(size=18),
            legend.text=element_text(size=18))+
      annotate("text", x = 400, y = 0.85,label="PT Arm",size=6)+
      geom_hline(yintercept=range_value[1],lty="dotted",size=2)+
      geom_hline(yintercept=range_value[2],lty="dotted",size=2)
    
    CI_plot_crl
    
    ggsave(plot = CI_plot_crl, 
           width = 8, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots/CI_plot/CI_crl_logit_weight2.pdf")
    
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
    alpha <- -0.8
    refer_mean <- SIM_treat_03 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_treat_03  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_treat_03  <- rbind(t_SIM_treat_03, temp)
    }
    
    plot_SIM_treat_df3 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                     Time = c(75, 80, 85, 90, 95, 100),
                                     Estimate = t_SIM_treat_03[, 1],
                                     Lower_ci = t_SIM_treat_03[, 2],
                                     Upper_ci = t_SIM_treat_03[, 3])
    
    
    ##### month = 6 #####
    month <- 6
    alpha <- -0.8
    refer_mean <- SIM_treat_06 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_treat_06  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_treat_06  <- rbind(t_SIM_treat_06, temp)
    }
    
    plot_SIM_treat_df6 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                     Time = c(165, 170, 175, 180, 185, 190),
                                     Estimate = t_SIM_treat_06[, 1],
                                     Lower_ci = t_SIM_treat_06[, 2],
                                     Upper_ci = t_SIM_treat_06[, 3])
    
    
    ##### month = 9 #####
    month <- 9
    alpha <- -0.8
    refer_mean <- SIM_treat_09 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_treat_09  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_treat_09  <- rbind(t_SIM_treat_09, temp)
    }
    
    plot_SIM_treat_df9 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                     Time = c(255, 260, 265, 270, 275, 280),
                                     Estimate = t_SIM_treat_09[, 1],
                                     Lower_ci = t_SIM_treat_09[, 2],
                                     Upper_ci = t_SIM_treat_09[, 3])
    
    ##### month = 12 #####
    month <- 12
    alpha <- -0.8
    refer_mean <- SIM_treat_12 |> qlogis()
    
    var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
    mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
    t_SIM_treat_12  <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
    
    
    for(alpha in Alpha_seq[-1]){
      var_paraboot <- paraboot_var(data = data_paraboot_trt, alpha, month, Alpha_seq = Alpha_seq, link = "logit")
      mean_orig <- refer_mean[which(Alpha_seq == alpha)] 
      temp <- Wald_CI(mean_orig, var_orig = var_paraboot) |> plogis()
      t_SIM_treat_12  <- rbind(t_SIM_treat_12, temp)
    }
    
    plot_SIM_treat_df12 <- data.frame(Alpha = c("-0.8","-0.4","0","0.4","0.8","1.2"),
                                      Time = c(345, 350, 355, 360, 365, 370),
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
      expand_limits(x=c(30,480),y=c(0,1))+
      scale_x_continuous(breaks=seq(90,450,by=90))+
      scale_colour_manual(values = c("#253494", "#41b6c4", 
                                     "#a1dab4", "#ffffb2",
                                     "#fd8d3c", "#bd0026"),
                          # "#f03b20",# "#fecc5c",# "#2c7fb8", # "grey50"),
                          breaks=c("1.2", "0.8", "0.4", "0", "-0.4", "-0.8")) + 
      # "original data"))+
      labs(x="Time (Days)", y="Rate of ACQ Score < 1.5", color=expression(alpha)) +
      labs(color=expression(alpha))+
      theme(axis.title=element_text(size=25),
            legend.title=element_text(size=20),
            axis.text=element_text(size=18),
            legend.text=element_text(size=18))+
      theme(legend.position = "none")+
      annotate("text", x = 400, y = 0.85,label="PT+HV Arm",size=6)+
      geom_hline(yintercept=range_value[1],lty="dotted",size=2)+
      geom_hline(yintercept=range_value[2],lty="dotted",size=2)
    
    CI_plot_trt
    
    ggsave(plot = CI_plot_trt, 
           width = 7, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots/CI_plot/CI_trt_logit_weight2.pdf")
    
  }
  
  
}


