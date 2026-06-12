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
#####################################################
#install.packages("devtools")
#devtools::install_github("magronah/power.nb")
#####################################################
filter_data_list       =   logmeanFit_list     =  list()
gam_fit_list           =   dispersionFit_list  =  list()
logfoldchangeFit_list  =   list()
########################################################################
dirs <- list.dirs(path = getwd(), 
                  full.names = FALSE, 
                  recursive = FALSE)[-c(1,2)]


data_names  <- setdiff(dirs,"results")
data_names  <- setdiff(data_names,"results_other_data")
data_names  <- setdiff(data_names,"figures")


ref_name_vec     =   c("Control","Soil","H", "H", "Dalhousie", 
                       "H", "H","H", "control",  
                       "plastic", "North America","COLD",
                       "H" ,"H" ,"H" , "nonIBD", "WTP","Y" ,"H" , 
                       "H" , "H" ,"H" ,"H" ,"H" , "plastic","plastic",
                       "plastic","H" ,"H" ,"plastic")
#############################################################################
for(i in 1:length(data_names)){

  ##Read count data and metadata
  data_path = data_names[i]
  
  data <- read.table(
    file.path(data_path, paste0(data_path,"_ASVs_table.tsv")),
    header = TRUE, sep = "\t",
    check.names = FALSE, comment.char = "")
  
  metadata <- read.table(
    file.path(data_path, paste0(data_path,"_metadata.tsv")),
    header = TRUE, sep = "\t",
    check.names = FALSE, comment.char = "")
  
  metadata <- metadata %>%
    setNames(c("SampleID", "Groups"))
  
  sample_sizes = table(metadata$Groups)
  #############################################################################
  # ###### Pre-filtering low abundant taxa
  filter_data_list[[i]]  =  filter_low_count(countdata      =  data,
                                             metadata       =  metadata,
                                             abund_thresh   =  abund_thresh,
                                             sample_thresh  =  sample_thresh,
                                             sample_colname =  "SampleID",
                                             group_colname  =  "Groups")
  # #############################################################################
  # ########## Fold change and Dispersion Estimation ############
  foldchange_est <- deseqfun(countdata      =   filter_data_list[[i]],
                             metadata       =   metadata,
                             alpha_level    =   alpha_level,
                             ref_name       =   ref_name_vec[i],
                             group_colname  =   "Groups",
                             sample_colname =   "SampleID")

  logfoldchange =  foldchange_est$deseq_estimate$log2FoldChange
  # 
  # # #############################################################################
  logmean    =  log2(rowMeans(filter_data_list[[i]]))
  quiet <- function(expr) {
    out <- suppressWarnings(suppressMessages(
      capture.output(res <- eval.parent(substitute(expr)))
    ))
    res}

  logmeanFit_list[[i]]  =   quiet(logmean_fit(logmean,
                                              sig      =  sig,
                                              max.comp =  max.comp,
                                              max.boot =  max.boot))
  # # #############################################################################
  logfoldchangeFit_list[[i]] <- quiet(logfoldchange_fit(logmean,
                                                        logfoldchange,
                                                        ncore      =  ncore,
                                                        max_sd_ord =  max_sd_ord,
                                                        max_np     =  max_np,
                                                        minval     =  minval,
                                                        maxval     =  maxval,
                                                        itermax    =  itermax,
                                                        NP         =  NP,
                                                        seed       =  seed))
  # # #############################################################################
  # # ######## Modelling the dispersion estimates   #######
  dispersion    =  foldchange_est$dispersion
  dispersionFit_list[[i]] =  dispersion_fit(dispersion, logmean)
  # # #############################################################################
  dispersion_param    =  dispersionFit_list[[i]]$param
  logmean_param       =  logmeanFit_list[[i]]$param
  logfoldchange_param =  logfoldchangeFit_list[[i]]
  # # #############################################################################
  countdata_sims_list  =  list()
  for(k in 1:length(nsample_vec)){
    countdata_sims_list[[k]]  =  countdata_sim_fun(logmean_param,
                                                   logfoldchange_param,
                                                   dispersion_param,
                                                   nsamp_per_group = nsample_vec[k],
                                                   ncont       =  NULL,
                                                   ntreat      =  NULL,
                                                   notu        =  notu,
                                                   nsim        =  nsim,
                                                   disp_scale  =  disp_scale,
                                                   max_lfc     =  max_lfc,
                                                   maxlfc_iter =  maxlfc_iter,
                                                   seed        =  seed)
  }

  names(countdata_sims_list) = paste0("sample_",nsample_vec)
  # # #############################################################################
  # # ## Estimate p-values associated to fold changes for each taxa 
  # # ## for simulated data per sample size
  desq_est_list  =  list()
  for(j in 1:length(countdata_sims_list)){
    countdata_list       =   countdata_sims_list[[j]]$countdata_list
    metadata_list        =   countdata_sims_list[[j]]$metadata_list
    desq_est_list[[j]]   =   deseq_fun_est(metadata_list  =   metadata_list,
                                           countdata_list =   countdata_list,
                                           alpha_level    =   alpha_level,
                                           group_colname  =   "Groups",
                                           sample_colname =   "Samples",
                                           num_cores      =   num_cores,
                                           ref_name       =   "control")

  }
  names(desq_est_list) = paste0("sample_",nsample_vec)
  # # #############################################################################
  # # #### Fit Generalized Additive Model (GAM) for power estimation 
  deseq_list = lapply(desq_est_list, function(x){
    read_data(x, "deseq_estimate") })

  pval_est_list <- lapply(deseq_list, function(sample_list) {
    lapply(sample_list, function(sim_df) {
      sim_df$padj}) })


  logfoldchange_list  =   read_data(countdata_sims_list,"logfoldchange_list")
  logmean_list        =   read_data(countdata_sims_list,"logmean_list")

  gam_fit_list[[i]]   =   power_fun_ss(pval_est_list,
                                       logmean_list,
                                       nsample_vec  = nsample_vec,
                                       logfoldchange_list,
                                       alpha_level  =  alpha_level)

  saveRDS(gam_fit_list[[i]],
          file = paste0("results_other_data/gam_fit_",
                        data_path,".rds"))

  saveRDS(filter_data_list[[i]],
          file = paste0("results_other_data/filter_data_",
                        data_path,".rds"))


  saveRDS(logmeanFit_list[[i]],
          file = paste0("results_other_data/logmeanFit_",
                        data_path,".rds"))

  saveRDS(dispersionFit_list[[i]],
          file = paste0("results_other_data/dispersionFit_",
                        data_path,".rds"))

  saveRDS(logfoldchangeFit_list[[i]],
          file = paste0("results_other_data/logfoldchangeFit_",
                        data_path,".rds"))

  saveRDS(logfoldchange,
          file = paste0("results_other_data/logfoldchangeEst_",
                        data_path,".rds"))

  saveRDS(logmean,
          file = paste0("results_other_data/logmeanEst_",
                        data_path,".rds"))
  
  saveRDS(sample_sizes, 
          file = paste0("results_other_data/sample_size_",
                        data_path,".rds"))
  
}

saveRDS(ref_name_vec, 
        file = paste0("results_other_data/reference_names.rds"))
#############################################################################
names(gam_fit_list)           =   data_names
names(filter_data_list)       =   data_names
names(logmeanFit_list)        =   data_names
names(dispersionFit_list)     =   data_names
names(logfoldchangeFit_list)  =   data_names
#############################################################################
saveRDS(gam_fit_list, file = "results_other_data/gam_fit_list.rds")
saveRDS(filter_data_list, file = "results_other_data/filter_data_list.rds")
saveRDS(logmeanFit_list, file = "results_other_data/logmeanFit_list.rds")
saveRDS(dispersionFit_list, file = "results_other_data/dispersionFit_list.rds")
saveRDS(logfoldchangeFit_list, file = "results_other_data/logfoldchangeFit_list.rds")
#############################################################################
