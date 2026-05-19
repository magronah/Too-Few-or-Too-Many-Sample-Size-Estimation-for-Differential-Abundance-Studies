library(tidyverse)
library(dplyr)
library(power.nb)
library(patchwork)
library(readr)
library(rlist)
library(latex2exp)
library(scam)
library(purrr)
library(scales)
library(grid)
library(ggnewscale)
#############################################################################
source("functions.R")
source("param_vals.R")
#############################################################################
gam_fit_list           =   readRDS("results/gam_fit_list.rds")
logmeanFit_list        =   readRDS("results/logmeanFit_list.rds")
dispersionFit_list     =   readRDS("results/dispersionFit_list.rds")
logfoldchangeFit_list  =   readRDS("results/logfoldchangeFit_list.rds")
##############################################################################
data_names       =   c("obesity","blueberry soil","glass-plastic")
nsim = 5; notu = 1000
nsample_vec =  seq(10,150,20)
grid_len =  20  
#############################################################################
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
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}
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
cont_breaks <- c(0.2, 0.6,0.8) #seq(0.2, 0.8, 0.2)
ddplt  = plot_fun2(pred_df,
                   power_df,
                   cont_breaks, facet_labs)

ddplt

#############################################################################
target_powers <- seq(0.1, 0.99, 0.05)
lfc_vals <- c(1.3, 1.5, 2)
lmc_vals <- c(-1, 0, 5)

ss_vary_lfc <- do.call(
  rbind,
  lapply(seq_along(gam_fit_list), function(i) {
    sample_size_grid(
      target_powers = target_powers,
      vary_vals = lfc_vals,
      mod = gam_fit_list[[i]]$gam_mod,
      vary = "abs_lfc",
      fixed_logmean = 5
    ) |>
      transform(dataset = data_names[[i]])
  })
)

ss_vary_lmc <- do.call(
  rbind,
  lapply(seq_along(gam_fit_list), function(i) {
    sample_size_grid(
      target_powers = target_powers,
      vary_vals = lmc_vals,
      mod = gam_fit_list[[i]]$gam_mod,
      vary = "logmean",
      fixed_abs_lfc = 1.5
    ) |>
      transform(dataset = data_names[[i]])
  })
)

#############################################################################
file_names       =   c("ob_ross","Blueberry",
                       "glass_plastic_oberbeckmann")

ref_name_vec     =   c("H","Soil","plastic")
#############################################################################
actual_n <- map_dfr(file_names, function(data_path) {
  
  metadata <- read.table(
    file.path(data_path, paste0(data_path, "_metadata.tsv")),
    header = TRUE, sep = "\t",
    check.names = FALSE, comment.char = ""
  )
  
  counts <- table(metadata$comparison)
  
   data.frame(file_name = data_path, min_n = min(counts),
                  max_n = max(counts))
})

actual_n$dataset = data_names
#############################################################################
p1 <- plot_power_curve(
  data = ss_vary_lfc,
  group_var = abs_lfc,
  nn = 14,
  actual_n  =  actual_n,
  legend_label = TeX("$|\\log_2(fold \\ change)$|"))

p2 <- plot_power_curve(
  data = ss_vary_lmc,
  actual_n  =  actual_n,
  nn = 14,
  group_var = logmean,
  legend_label = TeX("$\\log_2(mean \\ count)$")
)

ggsave(
  filename = "figures/ss1.pdf",
  plot = p1,
  width = 17,
  height = 5,
  units = "in",
  device = cairo_pdf
)


ggsave(
  filename = "figures/ss2.pdf",
  plot = p2,
  width = 17,
  height = 5,
  units = "in",
  device = cairo_pdf
)
p1/p2
#####################################################################




