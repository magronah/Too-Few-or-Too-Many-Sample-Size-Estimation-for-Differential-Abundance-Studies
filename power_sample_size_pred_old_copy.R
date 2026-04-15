library(tidyverse)
library(dplyr)
library(power.nb)
library(patchwork)
library(readr)
library(rlist)
library(latex2exp)
library(scam)
library(ggnewscale)
source("functions.R")
source("param_vals.R")
library(scales)
library(grid)
#############################################################################
gam_fit_list           =   readRDS("results/gam_fit_list.rds")
logmeanFit_list        =   readRDS("results/logmeanFit_list.rds")
dispersionFit_list     =   readRDS("results/dispersionFit_list.rds")
logfoldchangeFit_list  =   readRDS("results/logfoldchangeFit_list.rds")
##############################################################################
data_names       =   c("Obesity","Blueberry Soil","Glass-Plastic")
nsim = 5; notu = 1000
nsample_vec =  seq(10,150,20)
grid_len =  20  

power_dd_list  =  pred_dd_list = list()
logfoldchange_sim_list  =  logmean_sim_list  = list()

for(i in 1:length(data_names)){
  ##############################################################################
  dispersion_param      =  dispersionFit_list[[i]]$param
  logmean_param         =  logmeanFit_list[[i]]$param
  logfoldchange_param   =  logfoldchangeFit_list[[i]]
  ##############################################################################
  countdata_sims_list  =  list()
  for(j in 1:length(nsample_vec)){
    countdata_sims_list[[j]]  =  countdata_sim_fun(logmean_param,
                                                   logfoldchange_param,
                                                   dispersion_param,
                                                   nsamp_per_group = nsample_vec[j],
                                                   ncont  = NULL,
                                                   ntreat = NULL,
                                                   notu   = notu,
                                                   nsim   = nsim,
                                                   disp_scale = disp_scale,
                                                   max_lfc    = max_lfc,
                                                   maxlfc_iter = maxlfc_iter,
                                                   seed = 131)
  }
  names(countdata_sims_list) = paste0("sample_",nsample_vec)
  ##############################################################################
  logfoldchange_list  =   read_data(countdata_sims_list,"logfoldchange_list")
  logmean_list        =   read_data(countdata_sims_list,"logmean_list")
  logfoldchange_sim_list[[i]]  =  logfoldchange_list
  logmean_sim_list[[i]]        =  logmean_list
  ##############################################################################
  desq_est_list   =  list()
  for(k in 1:length(countdata_sims_list)){
    countdata_list       =    countdata_sims_list[[k]]$countdata_list
    metadata_list        =    countdata_sims_list[[k]]$metadata_list
    desq_est_list[[k]]   =    deseq_fun_est(metadata_list  =  metadata_list,
                                            countdata_list =  countdata_list,
                                            alpha_level    =  alpha_level,
                                            group_colname  =  "Groups",
                                            sample_colname =  "Samples",
                                            num_cores      =   num_cores,
                                            ref_name       =  "control")
    
  }
  names(desq_est_list) = paste0("sample_",nsample_vec)
  ##############################################################################
  deseq_list = lapply(desq_est_list, function(x){
    read_data(x, "deseq_estimate")
  })
  
  pval_est_list <- lapply(deseq_list, function(sample_list) {
    lapply(sample_list, function(sim_df) {
      sim_df$padj
    })
  })
  
  p_val         =   unname(unlist(pval_est_list))
  pval_reject   =   (!is.na(p_val) & p_val < alpha_level)
  ##############################################################################
  sample_size         =   rep(nsample_vec,
                              times = sapply(lapply(logfoldchange_list, unlist),
                                             length))
  ##############################################################################
  pred_data = data.frame(logmean  =  unlist(logmean_list),
                         abs_lfc  =  abs(unlist(logfoldchange_list)),
                         sample_size   =  sample_size,
                         logsample_size = log(sample_size))
  
  pred_data$pvalue_reject    =   factor(as.numeric(pval_reject))
  rownames(pred_data)        =   NULL
  ##############################################################################
  model                  =    gam_fit_list[[i]]$gam_mod
  power_estimate         =    power_grid(pred_data, model, grid_len)
  pred_data$dataset      =    data_names[i]
  power_estimate$dataset =    data_names[i]
  power_dd_list[[i]]     =    power_estimate
  pred_dd_list[[i]]      =    pred_data
}

names(power_dd_list)    =   data_names
names(pred_dd_list)     =   data_names
names(logfoldchange_sim_list)  =  data_names
names(logmean_sim_list)  =  data_names
################################################################
saveRDS(power_dd_list, file = "results/power_dd_list.rds")
saveRDS(pred_dd_list, file = "results/pred_dd_list.rds")
saveRDS(logfoldchange_sim_list, file = "results/logfoldchange_sim_list.rds")
saveRDS(logmean_sim_list, file = "results/logmean_sim_list.rds")
################################################################
nsample_vec
sub_samples = c(30, 70, 150)
pred_df <- do.call(rbind, lapply(pred_dd_list, function(df) {
  df[df$sample_size %in% sub_samples, ]
}))

power_df <- do.call(rbind, lapply(power_dd_list, function(df) {
  df[df$sample_size %in% sub_samples, ]
}))
#############################################################################
cont_breaks <- seq(0.2, 0.8, 0.2)
cont_breaks <- seq(0.2, 0.8, 0.5)
ddplt  = plot_fun2(pred_df,
                   power_df,
                   cont_breaks, facet_labs)
ddplt
#############################################################################
target_powers <- seq(0.1, 0.99, 0.05)
lfc_vals <- c(1.3, 1.5, 2)
ss_all <- do.call(
  rbind,
  lapply(seq_along(gam_fit_list), function(i) {
    sample_size_est(target_powers, lfc_vals, gam_fit_list[[i]]$gam_mod) |>
      transform(dataset = data_names[[i]])
  })
)


ss_var_lmc <- function(target_powers,lmc_vals, mod){
  dd_pow <- do.call(rbind, lapply(lmc_vals, function(lmc) {
    
    ss_vals <- sapply(target_powers, function(tp) {
      ss_solver(
        target_power = tp,
        logmean = lmc,
        abs_lfc = 1.5,
        model = mod,
        xmin = log2(5),
        xmax = log2(500)
      )
    })
    
    data.frame(
      sample_size = ss_vals,
      power = target_powers,
      logmean = lmc
    )
  }))
  
}
#############################################################################
actual_n <- data.frame(
  dataset  =  data_names,
  n_actual =  c(25, 40, 18)
)
#############################################################################
p1 <- plot_power_curve(
  data = ss_all,
  group_var = abs_lfc,
  actual_n  =  actual_n,
  legend_label = TeX("$|\\log_2(fold \\ change)$|")
)

p2 <- plot_power_curve(
  data = ss_all_lmc,
  actual_n  =  actual_n,
  group_var = logmean,
  legend_label = TeX("$\\log_2(mean \\ count)$")
)

p1/p2

#####################################################################
lmc_vals <- c(-1, 0, 5)



ss_lmc1 = ss_var_lmc(target_powers,lmc_vals, gam_fit_list[[1]]$gam_mod)
ss_lmc1$dataset =  data_names[[1]]


ss_lmc2 = ss_var_lmc(target_powers,lmc_vals, gam_fit_list[[2]]$gam_mod)
ss_lmc2$dataset =  data_names[[2]]

ss_lmc3 = ss_var_lmc(target_powers,lmc_vals, gam_fit_list[[3]]$gam_mod)
ss_lmc3$dataset =  data_names[[3]]


ss_all_lmc = rbind(ss_lmc1, ss_lmc2, ss_lmc3)

nn= 15
#####################################################################


p1 + scale_y_continuous(trans = "logit")

p2 = ggplot(dd_pow, aes((sample_size), power, color = factor(abs_lfc))) +
  geom_point(size=2) +
  geom_line(lwd =1) +
  theme_bw() +
  labs(
    x = "Sample size per group",
    y = "Target power",
    color = "log2(mean count)"
  )+
  theme(
    #strip.text = element_text(size = nn, color = "black"),
    legend.title = element_text(size = nn, color = "black"),
    legend.text = element_text(size = nn, color = "black"),
    axis.text.x = element_text(size = nn, color = "black"),
    axis.text.y = element_text(size = nn, color = "black"),
    axis.title.x = element_text(size = nn, color = "black"),
    axis.title.y = element_text(size = nn, color = "black")
  )

p2
p1|p2
target_powers <- seq(0.1,0.99,0.05)

ss_vals <- sapply(target_powers, function(tp) {
  ss_solver(
    target_power = tp,
    logmean = log2(7),
    abs_lfc = log2(2.5),
    model = mod,
    xmin = log2(5),
    xmax = log2(500)
  )
})

dd_pow = data.frame(sample_size = ss_vals, 
                    power = target_powers)

ggplot(dd_pow, aes(sample_size, power)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "Sample size per group",
       y = "Target power") 



#################################################################
ggplot(dd_pow, aes(sample_size, power)) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  theme_bw()


cont_breaks2 = 0.6

power_estimate$sample_size <- factor(power_estimate$sample_size,
                                     levels = c(25, 75, 125, 150))

dd$pvalue_reject <- factor(dd$pvalue_reject,
                           levels = c(0, 1),
                           labels = c("No", "Yes"))

#ggrastr::rasterise(
#  geom_point(aes(colour = pvalue_reject), alpha = 0.35),
#  dpi = 300
#) +
# scale_colour_manual(
#  values = c("No" = "black", "Yes" = "red"),
#  name = "Significant Taxa"
#) +


one_break <- 0.2
gg_2dimc_one <- ggplot(dd, aes(logmean, abs_lfc)) +
  ggnewscale::new_scale_colour() +
  geom_contour(
    data = power_estimate,
    aes(
      x = logmean,
      y = abs_lfc,
      z = power,
      colour = factor(sample_size)
    ),
    breaks = one_break,
    linewidth = 1
  ) +
  metR::geom_text_contour(
    data = power_estimate,
    aes(
      x = logmean,
      y = abs_lfc,
      z = power,
      label = after_stat(level),
      colour = factor(sample_size)),
    breaks = one_break,
    stroke = 0.15,
    check_overlap = FALSE,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "25"  = "#F8766D",
      "75"  = "#7CAE00",
      "125" = "#00BFC4",
      "150" = "#C77CFF"
    ),
    name = "Sample size",
    labels = c(
      "25"  = "25 samples per group",
      "75"  = "75 samples per group",
      "125" = "125 samples per group",
      "150" = "150 samples per group"
    )
  ) +
  xlab(TeX("$\\log_2$(mean counts)")) +
  ylab(TeX("|$\\log_2$(fold change)|")) + 
  theme_bw(base_size = 16)

gg_2dimc_one



