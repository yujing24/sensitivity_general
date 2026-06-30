rm(list=ls())
library(ggplot2)
library(metR)
library(isoband)
library(dplyr)

setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")



### weight1
{
  ##############   Month 6  ############## 
  {
    
    df_m6 <- readRDS("result/contour_plot/Contour_plot_month6_SIM_logit.rds")
    
    SIM_p1 <- ggplot(df_m6, aes(x = Alpha_crl, y = Alpha_trt, z = Estimate_vec))+
      geom_contour_fill(binwidth=0.06, show.legend = F) +
      scale_fill_divergent(mid="white",low="forestgreen",high="orange") +
      geom_contour2(aes(z=Estimate_vec,label=..level..),
                    colour="black", size=0.5, binwidth=0.06, skip=0.06,
                    margin=grid::unit(c(5,5,5,5),"pt"),
                    label.placer = label_placer_fraction(frac=0.38,
                                                         rot_adjuster=angle_fixed()),
                    label_size=7) + # in package we will remove this geom_contour2() function
      # geom_contour2(data = df2_m6_restrict, aes(x=Alpha_crl, y=Alpha_trt, z=CI_info, label=..level..)) +
      theme_bw() +
      scale_colour_divergent(mid="yellow",low="forestgreen",high="orange")+
      labs(x="Alpha for PT Arm",y="Alpha for PT+HV Arm")+
      theme(axis.title = element_text(size=25),
            axis.text = element_text(size=18))+
      geom_vline(xintercept=0,lty="dotted",size=1)+
      geom_hline(yintercept=0,lty="dotted",size=1)+
      geom_point(x=0, y=0, size=3)+
      scale_x_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))+
      scale_y_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))
    
    SIM_p1
    
    ggsave(plot = SIM_p1, 
           width = 7, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots/contour_plot/contour_month6_mean_logit.pdf")
    
    
    df2_m6 <- mutate(df_m6, CI_info = Lower_ci*(Lower_ci>0) + Upper_ci*(Upper_ci<0))
    
    min(df2_m6$CI_info)
    max(df2_m6$CI_info)
    
    SIM_p2 <- ggplot(df2_m6, aes(x = Alpha_crl, y=Alpha_trt, z=CI_info))+
      geom_contour_fill(breaks = c(seq(-0.33, 0.33, by = 0.06))) +
      scale_fill_divergent(mid="white",low="forestgreen",high="orange",
                           midpoint = 0) +
      theme_bw()+
      theme(legend.title = element_blank())+
      labs(x="Alpha for PT Arm",y="Alpha for PT+HV Arm")+
      theme(axis.title = element_text(size=25),
            axis.text = element_text(size=18),
            legend.text=element_text(size=14),
            legend.key.size = unit(1,"cm"))+
      geom_vline(xintercept=0,lty="dotted",size=1)+
      geom_hline(yintercept=0,lty="dotted",size=1)+
      geom_point(x=0,y=0,size=3)+
      scale_x_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))+
      scale_y_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))
    
    SIM_p2
    
    ggsave(plot = SIM_p2, 
           width = 8, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots/contour_plot/contour_month6_diff_logit.pdf")
    
    
  }
  
  
  ##############   Month 12  ############## 
  {
    
    
    #################################
    df_m12 <- readRDS("result/contour_plot/Contour_plot_month12_SIM_logit.rds")
    
    
    SIM_p3 <- ggplot(df_m12,aes(x=Alpha_crl, y=Alpha_trt, z=Estimate_vec))+
      geom_contour_fill(binwidth=0.06, show.legend = F) +
      scale_fill_divergent(mid="white",low="forestgreen",high="orange") +
      geom_contour2(aes(z=Estimate_vec,label=..level..),
                    colour="black", size=0.5, binwidth=0.06, skip=0.06,
                    margin=grid::unit(c(5,5,5,5),"pt"),
                    label.placer = label_placer_fraction(frac=0.38,
                                                         rot_adjuster=angle_fixed()),
                    label_size=7) + 
      theme_bw()+
      scale_colour_divergent(mid="yellow",low="forestgreen",high="orange")+
      labs(x="Alpha for PT Arm",y="Alpha for PT+HV Arm")+
      theme(axis.title = element_text(size=25),
            axis.text = element_text(size=18))+
      geom_vline(xintercept=0,lty="dotted",size=1)+
      geom_hline(yintercept=0,lty="dotted",size=1)+
      geom_point(x=0, y=0, size=3)+
      scale_x_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))+
      scale_y_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))
    
    SIM_p3
    
    ggsave(plot = SIM_p3, 
           width = 7, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots/contour_plot/contour_month12_mean_logit.pdf")
    
    
    df2_m12 <- mutate(df_m12, CI_info = Lower_ci*(Lower_ci>0) + Upper_ci*(Upper_ci<0))
    
    min(df2_m12$CI_info)
    max(df2_m12$CI_info)
    
    SIM_p4 <- ggplot(df2_m12,aes(x=Alpha_crl, y=Alpha_trt, z=CI_info))+
      geom_contour_fill(breaks = c(seq(-0.33, 0.33, by = 0.06))) +
      scale_fill_divergent(mid="white",low="forestgreen",high="orange",
                           midpoint = 0) +
      theme_bw()+
      theme(legend.title = element_blank())+
      labs(x="Alpha for PT Arm",y="Alpha for PT+HV Arm")+
      theme(axis.title = element_text(size=25),
            axis.text = element_text(size=18),
            legend.text=element_text(size=14),
            legend.key.size = unit(1,"cm"))+
      geom_vline(xintercept=0,lty="dotted",size=1)+
      geom_hline(yintercept=0,lty="dotted",size=1)+
      geom_point(x=0,y=0,size=3)+
      scale_x_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))+
      scale_y_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))
    
    SIM_p4
    
    ggsave(plot = SIM_p4, 
           width = 8, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots/contour_plot/contour_month12_diff_logit.pdf")
    
    
  }
}



### weight2
{
  ##############   Month 6  ############## 
  {
    
    df_m6 <- readRDS("result/contour_plot/Contour_plot_month6_SIM_logit_weight2.rds")
    
    SIM_p1 <- ggplot(df_m6, aes(x = Alpha_crl, y = Alpha_trt, z = Estimate_vec))+
      geom_contour_fill(binwidth=0.06, show.legend = F) +
      scale_fill_divergent(mid="white",low="forestgreen",high="orange") +
      geom_contour2(aes(z=Estimate_vec,label=..level..),
                    colour="black", size=0.5, binwidth=0.06, skip=0.06,
                    margin=grid::unit(c(5,5,5,5),"pt"),
                    label.placer = label_placer_fraction(frac=0.38,
                                                         rot_adjuster=angle_fixed()),
                    label_size=7) + # in package we will remove this geom_contour2() function
      # geom_contour2(data = df2_m6_restrict, aes(x=Alpha_crl, y=Alpha_trt, z=CI_info, label=..level..)) +
      theme_bw() +
      scale_colour_divergent(mid="yellow",low="forestgreen",high="orange")+
      labs(x="Alpha for PT Arm",y="Alpha for PT+HV Arm")+
      theme(axis.title = element_text(size=25),
            axis.text = element_text(size=18))+
      geom_vline(xintercept=0,lty="dotted",size=1)+
      geom_hline(yintercept=0,lty="dotted",size=1)+
      geom_point(x=0, y=0, size=3)+
      scale_x_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))+
      scale_y_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))
    
    SIM_p1
    
    ggsave(plot = SIM_p1, 
           width = 7, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots/contour_plot/contour_month6_mean_logit_weight2.pdf")
    
    
    df2_m6 <- mutate(df_m6, CI_info = Lower_ci*(Lower_ci>0) + Upper_ci*(Upper_ci<0))
    
    min(df2_m6$CI_info)
    max(df2_m6$CI_info)
    
    SIM_p2 <- ggplot(df2_m6, aes(x = Alpha_crl, y=Alpha_trt, z=CI_info))+
      geom_contour_fill(breaks = c(seq(-0.33, 0.33, by = 0.06))) +
      scale_fill_divergent(mid="white",low="forestgreen",high="orange",
                           midpoint = 0) +
      theme_bw()+
      theme(legend.title = element_blank())+
      labs(x="Alpha for PT Arm",y="Alpha for PT+HV Arm")+
      theme(axis.title = element_text(size=25),
            axis.text = element_text(size=18),
            legend.text=element_text(size=14),
            legend.key.size = unit(1,"cm"))+
      geom_vline(xintercept=0,lty="dotted",size=1)+
      geom_hline(yintercept=0,lty="dotted",size=1)+
      geom_point(x=0,y=0,size=3)+
      scale_x_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))+
      scale_y_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))
    
    SIM_p2
    
    ggsave(plot = SIM_p2, 
           width = 8, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots/contour_plot/contour_month6_diff_logit_weight2.pdf")
    
    
  }
  
  
  ##############   Month 12  ############## 
  {
    
    
    #################################
    df_m12 <- readRDS("result/contour_plot/Contour_plot_month12_SIM_logit_weight2.rds")
    
    
    SIM_p3 <- ggplot(df_m12,aes(x=Alpha_crl, y=Alpha_trt, z=Estimate_vec))+
      geom_contour_fill(binwidth=0.06, show.legend = F) +
      scale_fill_divergent(mid="white",low="forestgreen",high="orange") +
      geom_contour2(aes(z=Estimate_vec,label=..level..),
                    colour="black", size=0.5, binwidth=0.06, skip=0.06,
                    margin=grid::unit(c(5,5,5,5),"pt"),
                    label.placer = label_placer_fraction(frac=0.38,
                                                         rot_adjuster=angle_fixed()),
                    label_size=7) + 
      theme_bw()+
      scale_colour_divergent(mid="yellow",low="forestgreen",high="orange")+
      labs(x="Alpha for PT Arm",y="Alpha for PT+HV Arm")+
      theme(axis.title = element_text(size=25),
            axis.text = element_text(size=18))+
      geom_vline(xintercept=0,lty="dotted",size=1)+
      geom_hline(yintercept=0,lty="dotted",size=1)+
      geom_point(x=0, y=0, size=3)+
      scale_x_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))+
      scale_y_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))
    
    SIM_p3
    
    ggsave(plot = SIM_p3, 
           width = 7, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots/contour_plot/contour_month12_mean_logit_weight2.pdf")
    
    
    df2_m12 <- mutate(df_m12, CI_info = Lower_ci*(Lower_ci>0) + Upper_ci*(Upper_ci<0))
    
    min(df2_m12$CI_info)
    max(df2_m12$CI_info)
    
    SIM_p4 <- ggplot(df2_m12,aes(x=Alpha_crl, y=Alpha_trt, z=CI_info))+
      geom_contour_fill(breaks = c(seq(-0.33, 0.33, by = 0.06))) +
      scale_fill_divergent(mid="white",low="forestgreen",high="orange",
                           midpoint = 0) +
      theme_bw()+
      theme(legend.title = element_blank())+
      labs(x="Alpha for PT Arm",y="Alpha for PT+HV Arm")+
      theme(axis.title = element_text(size=25),
            axis.text = element_text(size=18),
            legend.text=element_text(size=14),
            legend.key.size = unit(1,"cm"))+
      geom_vline(xintercept=0,lty="dotted",size=1)+
      geom_hline(yintercept=0,lty="dotted",size=1)+
      geom_point(x=0,y=0,size=3)+
      scale_x_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))+
      scale_y_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8, 1.2))
    
    SIM_p4
    
    ggsave(plot = SIM_p4, 
           width = 8, 
           height = 6, 
           dpi = 500, 
           filename = "result/plots/contour_plot/contour_month12_diff_logit_weight2.pdf")
    
    
  }
}
