library(tidyverse)
library(dplyr)
library(here)
#remove.packages("power.nb")
#library(devtools)
#install_github("magronah/power.nb")
library(power.nb)
library(patchwork)
#####################################################
dir_vec    =  c("Blueberry","glass_plastic_oberbeckmann",
                "ob_goodrich","ob_ross", "ob_turnbaugh", "ob_zhu")

data_list  =  metadata_list  =  list()

for(i in 1:length(dir_vec)){
  
  data_list[[i]]       =  read.table(paste0(dir_vec[i], "/",
                                            dir_vec[i], "_ASVs_table.tsv"),
                                     header = TRUE, sep = "\t",
                                     check.names = FALSE, comment.char = "")
  
  
  metadata_list[[i]]   =  read.table(paste0(dir_vec[i], "/",
                                            dir_vec[i], "_metadata.tsv"),
                                     header = TRUE, sep = "\t",
                                     check.names = FALSE, comment.char = "")
  
}

names(data_list)  =  names(metadata_list)  =  dir_vec
lapply(data_list, dim)

dim(data_list[[2]])
sum(metadata_list[[2]]$comparison == "plastic")

#####################################################
### Filter out low abundant taxa
filter_data_list  =  list()
for(i in 1:length(metadata_list)){
  filter_data_list[[i]]  = filter_low_count(countdata = data_list[[i]],
                                            metadata  = metadata_list[[i]],
                                            abund_thresh = 10,
                                            sample_thresh = 3,
                                            sample_colname = "sampleid",
                                            group_colname  = "comparison")
}

lapply(data_list, dim)

names(filter_data_list)  =  dir_vec
#####################################################
#Fit log mean count
logmean_list    =  logmeanFit_list  =  list()

for(i in 1:length(filter_data_list)){
  logmean_list[[i]]    =  log(rowMeans(filter_data_list[[i]] ))
  logmeanFit_list[[i]] =  logmean_fit(logmean_list[[i]] , sig = 0.05,
                                      max.comp = 4, max.boot = 100)
}

names(logmean_list) = names(logmeanFit_list) =  dir_vec
saveRDS(logmeanFit_list, file = "results/logmeanParam_list.rds")

#Estimate foldchange from Deseq
foldchange_est_list =  list()
for(i in 1:length(logmean_list)){
  foldchange_est_list[[i]] <- deseqfun(countdata = filter_data_list[[i]],
                                       metadata  = metadata_list[[i]],
                                       alpha_level = 0.1,
                                       ref_name  = NULL,
                                       group_colname = "comparison",
                                       sample_colname = "sampleid")
}

names(foldchange_est_list) =  dir_vec

logfoldchange_list  =  logfoldchangeFit_list = list()
dispersionFit_list  =  dispersion_list = list()

i=1
dispersion_fit(dispersion_list[[i]], logmean_list[[i]]) 
for(i in 1:length(logmean_list)){
  
  logfoldchange_list[[i]] =  foldchange_est_list[[i]]$deseq_estimate$log2FoldChange
  dispersion_list[[i]]    =  foldchange_est_list[[i]]$dispersion
  dispersionFit_list[[i]] =  dispersion_fit(dispersion_list[[i]], logmean_list[[i]])
  
  
  ## Fit foldchange from Deseq
  logmean = logmean_list[[i]]
  logfoldchange = logfoldchange_list[[i]]
  logfoldchangeFit_list[[i]] <- logfoldchange_fit(logmean,
                                                  logfoldchange,
                                                  ncore = 3,
                                                  max_sd_ord = 2,
                                                  max_np = 5,
                                                  minval = -5,
                                                  maxval = 5,
                                                  itermax = 100,
                                                  NP = 800,
                                                  seed = 100)
 }

foldchange_list  =  lapply(foldchange_est_list, function(x){
                                             x$logfoldchange})

names(foldchange_list)    =  names(foldchange_est_list) = dir_vec
names(dispersionFit_list) =  names(logfoldchangeFit_list) =  dir_vec

saveRDS(foldchange_est_list, file = "results/foldchange_all_est_list.rds")
saveRDS(foldchange_list, file = "results/foldchange_est_list.rds")
saveRDS(logfoldchangeFit_list, file = "results/logfoldchangeParam_list.rds")
saveRDS(dispersionFit_list, file = "results/dispersionParam_list.rds")
###############################################################
logmeanParam_list     =  readRDS("results/logmeanParam_list.rds")
logfoldchangeParam_list =  readRDS("results/logfoldchangeParam_list.rds")
dispersionParam_list    =  readRDS("results/dispersionParam_list.rds")
foldchange_est_list   =  readRDS("results/foldchange_est_list.rds")
###########################################################
## Simualate data
nsim  =  100
notu_list = lapply(filter_data_list, nrow)
countdata_sims_list = list()

for(i in 1:length(metadata_list)){
  
  grps   =  as.numeric(table(metadata_list[[i]]$comparison))
  ncont  =  grps[[1]] 
  ntreat =  grps[[2]]
  
  notu  = notu_list[[i]]
  logmean_param        =  logmeanParam_list[[i]]$param
  dispersion_param     =  dispersionParam_list[[i]]$param
  logfoldchange_param  =  logfoldchangeParam_list[[i]]


  countdata_sims_list[[i]] = countdata_sim_fun(logmean_param,
                                               logfoldchange_param,
                                               dispersion_param,
                                               nsamp_per_group = NULL,
                                               ncont = ncont,
                                               ntreat = ntreat,
                                               notu,
                                               nsim = nsim,
                                               disp_scale = 0.3,
                                               max_lfc = 15,
                                               maxlfc_iter = 5000,
                                               seed = 121)
  
}


names(countdata_sims_list) = dir_vec

head(countdata_sims_list[[1]]$countdata_list$sim_1)

View(countdata_sims_list$Blueberry$logfoldchange_list$sim_1)
View(countdata_sims_list$Blueberry$logmean_list$sim_1)
View(countdata_sims_list$Blueberry$control_countdata_list$sim_1)
View(countdata_sims_list$Blueberry$treat_countdata_list$sim_1)
###############################################################
## Inspect the plot
plt_lfc_list  =   plt_mean_list  = list()
for(i in 1:length(metadata_list)){
  
  countdata_sims = countdata_sims_list[[i]]
  
  notu = notu_list[[i]]
  dd  = data.frame(logmean  = c(logmean_list[[i]],
                                countdata_sims$logmean_list$sim_1),
                   logfoldchange  = c(foldchange_est_list[[i]],
                                      countdata_sims$logfoldchange_list$sim_1),
                   type  = rep(c("observation","simulation"), each = notu) )
  
  plt_mean_list[[i]] = ggplot(dd, aes(x = logmean, group = type, 
                                      colour = type)) +
                              geom_density()  +
                              theme_bw()
  
  plt_lfc_list[[i]] = ggplot(dd, aes(x = logfoldchange, group = type, 
                                     colour = type)) +
                             geom_density()  +
                               theme_bw()
}


p = list()
for(i in 1:length(plt_lfc_list)){
  p[[i]] =(plt_mean_list[[i]]|plt_lfc_list[[i]]) + plot_layout(guides = "collect")
}
(p[[1]]/p[[2]]/p[[3]]) + plot_layout(guides = "collect")
(p[[4]]/p[[5]]/p[[6]]) + plot_layout(guides = "collect")
###############################################################
true_lmean_list = deseq_est_list = list()
true_lfoldchange_list = desq_est = list()

for(j in 1:length(filter_data_list)){
  countdata_sims  =   countdata_sims_list[[j]]
  countdata_list  =   countdata_sims$countdata_list
  metadata_list   =   countdata_sims$metadata_list
  desq_est[[j]]   =   deseq_fun_est(metadata_list =  metadata_list,
                               countdata_list =  countdata_list,
                               alpha_level    =  0.1,
                               group_colname = "Groups",
                               sample_colname = "Samples",
                               num_cores      =  4,
                               ref_name       = NULL)
  ###############################################################
  true_lmean_list[[j]]    =    countdata_sims$logmean_list
  deseq_est_list[[j]]     =    lapply(desq_est[[j]], function(x){x$deseq_estimate})
  true_lfoldchange_list[[j]] = countdata_sims$logfoldchange_list
}

###############################################################
gamFit_list = list()
for(i in 1:length(filter_data_list)){
  gamFit_list[[i]] <- gam_fit(deseq_est_list[[i]],
                              true_lfoldchange_list[[i]], 
                              true_lmean_list[[i]],
                              grid_len = 50,
                              alpha_level=0.1)
}

names(gamFit_list) = dir_vec
saveRDS(gamFit_list, file = "results/gamFit_list.rds")
###############################################################
cont_breaks     =  seq(0,1,0.1)
for(i in 1:length(filter_data_list)){
  
  combined_data   =  gamFit[[i]]$combined_data
  power_estimate  =  gamFit[[i]]$power_estimate
  
  contour_plot_list[[i]] <- contour_plot_fun(combined_data,
                                             power_estimate,
                                             cont_breaks)
}

((contour_plot_list[[1]]|contour_plot_list[[2]]|contour_plot_list[[3]])/
  (contour_plot_list[[4]]|contour_plot_list[[5]]|contour_plot_list[[6]])) + 
    plot_layout(guides = "collect")
###########################################################
#### Sample size calculation
nsim = 5
nsamp_vec = c(20, 40, 60, 80)
notu = 1000
countdata_sims_list  =  list()
i = 1
for(j in 1:length(nsamp_vec)){
  
  logmean_param        =  logmeanParam_list[[i]]$param
  dispersion_param     =  dispersionParam_list[[i]]$param
  logfoldchange_param  =  logfoldchangeParam_list[[i]]
  
  #for(j in 1:length(nsamp_vec)){
    
  countdata_sims_list[[j]]  =  countdata_sim_fun(logmean_param,
                                                 logfoldchange_param,
                                                 dispersion_param,
                                                 nsamp_per_group = nsamp_vec[j],
                                                 notu,
                                                 nsim = nsim,
                                                 disp_scale = 0.3,
                                                 max_lfc = 15,
                                                 maxlfc_iter = 5000,
                                                 seed = 121)
#}

}
names(countdata_sims_list) = paste0("sample_",nsamp_vec)
###############################################################
desq_est_list  =  list()
for(i in 1:length(countdata_sims_list)){
  
  countdata_list       =   countdata_sims_list[[i]]$countdata_list
  metadata_list        =   countdata_sims_list[[i]]$metadata_list
  desq_est_list[[i]]   =   deseq_fun_est(metadata_list =  metadata_list,
                               countdata_list =  countdata_list,
                               alpha_level    =  0.1,
                               group_colname  = "Groups",
                               sample_colname = "Samples",
                               num_cores      =  4,
                               ref_name       = "control")
  
}

names(desq_est_list) = paste0("sample_",nsamp_vec)
###############################################################
sample_size_vec   <- rep(nsamp_vec, each = notu*nsim)


true_logmean   <- unlist(read_data(countdata_sims_list,"logmean_list"), 
                           recursive = TRUE, use.names = FALSE)
true_logfoldchange  <- unlist(read_data(countdata_sims_list,"logfoldchange_list"), 
                        recursive = TRUE, use.names = FALSE)

est_list      <- do.call("c", desq_est_list)
p_val <- do.call("c", lapply(est_list, function(est) est$deseq_estimate$padj))
pval_reject   =   (!is.na(p_val) & p_val < 0.1)

comb   =   tibble(lmean_abund  =   true_logmean,
                  abs_lfc      =   abs(true_logfoldchange),
                  pval_reject  =   as.numeric(pval_reject),
                  sample_size  =   sample_size_vec)

library(scam)
###############################################################
#df =  length(sample_vec) -1 # degrees of freedom
#comb$lss = log2(comb$sample_size)

fit_3d <- scam(pval_reject ~ s(lmean_abund, abs_lfc,bs="tedmi") +
                 s(sample_size,lmean_abund,bs="tedmi") +
                 s(sample_size,abs_lfc,bs="tedmi"),
                data = comb, family = binomial)

list(combined_data=comb, gam_mod = fit_3d)
concurvity(fit_3d, full = TRUE)

grid_len = 50
cont_breaks = 0.7#seq(0.5,1, 0.1)

            contour_plot_fun(combined_data=comb,
                             power_estimate,
                             cont_breaks)
  
pp      =   with(comb,
                 expand.grid(lmean_abund = seq(min(lmean_abund),
                                               max(lmean_abund),
                                               length  = grid_len),
                             abs_lfc   =  seq(min(abs_lfc),
                                              max(abs_lfc),
                                              length  =  grid_len),
                             sample_size = seq(min(sample_size),
                                               max(sample_size),
                                               length  =  grid_len)))
pp
pp$power <- predict(fit_3d, newdata = pp,type = "response")
power_estimate = pp
combined_data=comb 
combined_data$pvalue_reject <- factor(combined_data$pval_reject)

power_estimate1 <- power_estimate %>%
  filter(sample_size %in% nsamp_vec)

cont_breaks = c(0.01,0.2,0.5,0.7)#seq(0.1,1, 0.1)

head(power_estimate)

library(dplyr)

# 1) Remove any NA/Inf that would drop all points in a panel
combined_data2 <- combined_data %>%
  filter(is.finite(lmean_abund), is.finite(abs_lfc))

# 2) First, make sure it plots WITHOUT rasterisation
gg_2dimc <- ggplot(combined_data2, aes(lmean_abund, abs_lfc)) +
  geom_point(alpha = 0.5) +
  scale_colour_manual(values = c("black", "red")) +
  geom_contour(
    data   = power_estimate1,
    aes(z = power, color = factor(sample_size)),
    linewidth = 0.8,
    breaks = cont_breaks
  ) +
  scale_colour_manual(values = c("blue", "red")) +
  metR::geom_label_contour(
    data = power_estimate1,
    aes(z = power, label = sprintf("%.3f", after_stat(level))),
    breaks = cont_breaks,
    check_overlap = TRUE
  ) +
  theme_bw()

gg_2dimc




gg_2dimc <- (ggplot(combined_data)
             + aes(lmean_abund, abs_lfc)
             + ggrastr::rasterise(geom_point( alpha = 0.5))
             #+ xlab(TeX("$\\log_2$(mean counts)"))
             #+ ylab(TeX("|$\\log_2$(fold change)|"))
             #+ scale_colour_manual(values = c("black", "red"))
             + geom_contour(data = power_estimate,
                            aes(z=power,color = as.factor(sample_size)),lwd=1,
                            breaks = cont_breaks)
             + metR::geom_label_contour(data = power_estimate,
                                        aes(z= power,label = sprintf("%.3f", after_stat(level))),
                                        breaks = cont_breaks)
             + theme_bw() 
             
)

gg_2dimc

