library(tidyverse)
library(purrr)
library(ggplot2)
library(dplyr)
library(scam)
library(forcats)
library(power.nb)
library(patchwork)
library(readr)
library(rlist)
library(latex2exp)
library(scales)
library(grid)
library(ggnewscale)

source("functions.R")
path <- file.path(getwd(), "results_other_data")
####################################################################
logmean_df  <-  read_to_df("^logmeanEst_", "logmeanEst_")
logmean_df  <-  reorder_by_median(logmean_df)

logfc_df    <-  read_to_df("^logfoldchangeEst_", "logfoldchangeEst_")
logfc_df    <-  reorder_by_median(logfc_df)

logmeanFit_list       <- read_to_list("^logmeanFit_", "logmeanFit_")
logfoldchangeFit_list <- read_to_list("^logfoldchangeFit_", "logfoldchangeFit_")
dispersionFit_list    <- read_to_list("^dispersionFit_", "dispersionFit_")
gam_fit_list          <- read_to_list("^gam_fit_", "gam_fit_")
sample_size_list      <- read_to_list("^sample_size_", "sample_size_")
reference_names      <-  readRDS(paste0(path,"/reference_names.rds"))

data_names    =    names(dispersionFit_list)

####################################################################
##' The distributions of median of mean count 
##' and fold change across  30 datasets
##' 

lmc_dist <- ggplot(logmean_df, aes(Dataset, value)) +
  geom_boxplot(fill = "grey85", colour = "black", alpha = 0.1,outlier.size = 1.5) +
  labs(x = "Dataset", y = expression(log[2]("mean abundance"))) +
  my_theme


abslfc_dist <- ggplot(logfc_df, aes(Dataset, abs(value))) +
  geom_boxplot(fill = "grey85", colour = "black", alpha = 0.1, outlier.size = 1.5) +
  labs(x = "Dataset", y = expression("|" * log[2]("fold change") * "|")) +
  scale_y_continuous() +
  my_theme


ggsave(
  filename = "figures/lmc_dist.eps",
  plot = lmc_dist,
  device = cairo_ps,  
  width = 10,
  height = 8,
  units = "in"
)

ggsave(
  filename = "figures/abslfc_dist.eps",
  plot = abslfc_dist,
  device = cairo_ps,  
  width = 10,
  height = 8,
  units = "in"
)


##########################################################################
logmean_med <- logmean_df %>%
  group_by(Dataset) %>%
  summarise(median = median(value, na.rm = TRUE)) %>%
  mutate(type = "logmean")

logfc_med <- logfc_df %>%
  group_by(Dataset) %>%
  summarise(median = median(value, na.rm = TRUE)) %>%
  mutate(type = "logfoldchange")
##########################################################################
med_df <- bind_rows(logmean_med, logfc_med)

pp_med <- ggplot(med_df, aes(x = type, y = median)) +
  geom_boxplot(
    width = 0.45,
    fill = "grey85",
    colour = "black",
    linewidth = 0.8,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 2,
    alpha = 0.7,
    colour = "black"
  ) +
  scale_x_discrete(
    labels = c(
      logmean = expression(log[2]("mean abundance")),
      logfoldchange = expression("|" * log[2]("fold change") * "|")
    )
  ) +
  labs(
    x = NULL,
    y = "Median value across datasets"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 15, colour = "black"),
    axis.title.x = element_text(size = 15, colour = "black"),
    
    
    axis.line = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    
    plot.margin = margin(10, 10, 10, 10)
  ) + scale_y_log10(
    breaks = c(0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2),
    labels = scales::label_number(accuracy = 0.01)
  )


ggsave(
  filename = "figures/med_both.eps",
  plot = pp_med,
  device = cairo_ps,  # important for proper rendering
  width = 9,
  height = 7,
  units = "in"
)
###########################################################################
abslfc <- summary(abs(logfc_med$median))
lmc    <- summary(logmean_med$median)

summary_abslfc_med <- data.frame(
  Min    =  as.numeric(abslfc["Min."]),
  Median =  as.numeric(abslfc["Median"]),
  Max = as.numeric(abslfc["Max."]))

summary_lmc_med <- data.frame(
  Min    =  as.numeric(lmc["Min."]),
  Median =  as.numeric(lmc["Median"]),
  Max = as.numeric(lmc["Max."]))
###########################################################################
saveRDS(summary_abslfc_med,
        file = paste0("results_other_data/summary_abslfc_median.rds"))
saveRDS(summary_lmc_med,
        file = paste0("results_other_data/summary_lmc_median.rds"))
###########################################################################
### plot sample sizes per group
control_n <- mapply(function(x, ref) {
  # match control group
  if (ref %in% names(x)) {
    return(x[ref])
  } else {
    return(NA)  # in case mismatch
  }
}, sample_size_list, reference_names)

treatment_n <- mapply(function(x, ref) {
  # take the other group
  other <- setdiff(names(x), ref)
  if (length(other) > 0) {
    return(x[other][1])
  } else {
    return(NA)
  }
}, sample_size_list, reference_names)

sample_size_df <- data.frame(
  Dataset = names(sample_size_list),
  min_n = sapply(sample_size_list, min),
  max_n = sapply(sample_size_list, max), 
  control = control_n,
  treatment = treatment_n,
  row.names = NULL
)

sample_size_long <- sample_size_df %>%
  pivot_longer(
    cols = c(min_n, max_n, control, treatment),
    names_to = "type",
    values_to = "value")

sample_size_long <- sample_size_long %>%
  group_by(Dataset) %>%
  mutate(total_n = sum(value)) %>%
  ungroup()

sample_size_long$Dataset <- reorder(sample_size_long$Dataset, 
                                    sample_size_long$total_n)

sample_size_cont_treat = sample_size_long[sample_size_long$type %in% 
                                     c("control","treatment"),]


pp_ss <- ggplot(sample_size_cont_treat, aes(x = Dataset, y = value, fill = type)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.75,
    color = "black",
    linewidth = 0.3,
    alpha = 0.9
  ) +
  scale_fill_manual(values = c(
    "control" = "#bdbdbd",
    "treatment" = "#252525"
  )) +
  scale_y_continuous(trans = "log10") +
  labs(
    x = "Dataset",
    y = "Sample size per group",
    fill = "Group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_line(color = "grey85")
  )


ggsave(
  filename = "figures/cont_treat_ss.eps",
  plot = pp_ss,
  device = cairo_ps,   
  width = 10,
  height = 8,
  units = "in"
)

######################################################################
sample_size_max_min = sample_size_long[sample_size_long$type %in% 
                                            c("min_n","max_n"),]

pp_sample <- ggplot(sample_size_max_min, aes(x = type, y = (value))) +
  geom_boxplot(
    width = 0.45,
    fill = "grey85",
    colour = "black",
    linewidth = 0.8,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 2,
    alpha = 0.7,
    colour = "black"
  ) +
  scale_x_discrete(
    labels = c(
      min_n = "Minimum sample size",
      max_n = "Maximum sample size"
    )
  ) +
  labs(
    x = NULL,
    y = "Sample size per group"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 15, colour = "black"),
    axis.title.x = element_text(size = 15, colour = "black"),
    
    axis.line = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    
    plot.margin = margin(10, 10, 10, 10)
  )  + scale_y_log10(
    breaks = c(5, 10, 20, 30, 50, 100, 200, 500, 1000),
    labels = c(5, 10, 20, 30, 50, 100, 200, 500, 1000)
  )


ggsave(
  filename = "figures/sample_sizes.eps",
  plot = pp_sample,
  device = cairo_ps,   
  width = 9,
  height = 7,
  units = "in"
)

######################################################################
##power plots

summary_max_n = summary(sample_size_max_min[sample_size_max_min$type 
                                            == "max_n",]$value)
summary_min_n = summary(sample_size_max_min[sample_size_max_min$type 
                                            == "min_n",]$value)

summary_max_n <- data.frame(
  Min    =  as.numeric(summary_max_n["Min."]),
  Median =  as.numeric(summary_max_n["Median"]),
  Max    = as.numeric(summary_max_n["Max."]))


summary_min_n <- data.frame(
  Min    =  as.numeric(summary_min_n["Min."]),
  Median =  as.numeric(summary_min_n["Median"]),
  Max    = as.numeric(summary_min_n["Max."]))

saveRDS(summary_max_n,
        file = paste0("results_other_data/summary_max_n.rds"))

saveRDS(summary_min_n,
        file = paste0("results_other_data/summary_min_n.rds"))
######################################################################
summary_datasets <- expand.grid(
  mean_level = c("Low mean count (-2.34)", "Median mean count (0.09)", 
                 "High mean count (2.20)"),
  foldchange  = c("Low (0.00)", "Median (0.07)", "High (0.66)"),
  n_level    = c("Small sample size (n=17)", "Median sample size (n=47)",
                 "Large sample size (n=798)"))

summary_datasets$logmean <- rep(as.numeric(summary_lmc_med), times = 9)
summary_datasets$abs_lfc <- rep(rep(as.numeric(summary_abslfc_med), each = 3), times = 3)
summary_datasets$sample_size <- rep(as.numeric(summary_max_n), each = 9)
summary_datasets$logsample_size <- log2(summary_datasets$sample_size)


power_df <- map_dfr(seq_along(gam_fit_list), function(i) {
  
  model <- gam_fit_list[[i]]$gam_mod
  
  pred <- predict(
    model,
    newdata = summary_datasets,
    type = "response"
  )
  
  cbind(
    Dataset = data_names[i],
    summary_datasets,
    power = pred
  )
})


power_df <- power_df |>
  mutate(Dataset_id = factor(as.numeric(factor(Dataset))))

pp_power <- ggplot(power_df, aes(x = foldchange, y = power, fill = foldchange)) +
  geom_boxplot(
    width = 0.65,
    alpha = 0.85,
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  geom_hline(
    yintercept = 0.8,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.7
  ) +
  ###
  geom_jitter(
    width = 0.15,
    height = 0,
    size = 1.8,
    alpha = 0.6,
    colour = "black"
  ) +
  #scale_y_log10() +
  scale_y_continuous(
    #trans = "log10"
    #breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1)
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(n_level ~ mean_level) +
  labs(
    x = expression("|" * log[2]("fold change") * "|"),
    y = "Predicted power" #(square-root scale)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(size = 11),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.title = element_text(size = 11),
    panel.spacing = unit(0.7, "lines")
  )

pp_power


ggsave(
  filename = "figures/power_plt.eps",
  plot = pp_power,
  device = cairo_ps,   
  width = 10,
  height = 8,
  units = "in"
)


######################################################################
summary_datasets2 <- expand.grid(
  mean_level = c("Low mean count (-2.34)", "Median mean count (-0.02)", 
                 "High mean count (2.20)"),
  foldchange  = c("Low (0.00)", "Median (0.05)", "High (0.66)"),
  power_level    = c("Power = 80%", "Power = 85%",
                 "Power = 90%"))

power_vec  <-  c(0.8, 0.85, 0.9)
summary_datasets2$logmean <- rep(as.numeric(summary_lmc_med), times = 9)
summary_datasets2$abs_lfc <- rep(rep(as.numeric(summary_abslfc_med), each = 3), times = 3)
summary_datasets2$target_power <- rep(as.numeric(power_vec), each = 9)


target_powers = c(0.8)


indx = setdiff(1:33, c(16,28,30))

res_list = list()
for(i in 1:30){
  model = gam_fit_list[[i]]$gam_mod
  res = uniroot_ss(target_power=0.8,logmean=max(as.numeric(summary_lmc_med)), 
                   abs_lfc=max(as.numeric(summary_abslfc_med)),
                   model,xmin = log2(10),
                     xmax = log2(5000),maxiter = 2000)
  
  res1 = uniroot_ss3(target_power=0.8,logmean=max(as.numeric(summary_lmc_med)), 
                     abs_lfc=max(as.numeric(summary_abslfc_med)),model,xmin = log2(10),
             xmax = log2(5000),maxiter = 2000)
  res_list[[i]] = res
  rr = c(res, res1$sample_size)
  print(rr)
}

############################################################
rrr = list()
for(i in 1:30){
  res=sample_size_grid(target_powers, vary_vals = max(as.numeric(summary_abslfc_med)), 
                       mod = gam_fit_list[[i]]$gam_mod,
                       vary =  "abs_lfc",
                       fixed_logmean = max(as.numeric(summary_lmc_med)),
                       fixed_abs_lfc = NULL,
                       xmin = log2(5),
                       xmax = log2(5000)) 
  rrr[[i]] = res$sample_size
  print(res$sample_size)
  
}


ddff = data.frame(Datasets = data_names,
                  est_ss = unlist(rrr))





summary_max_nn <- data.frame(
  stat = names(summary_max_n),
  n = as.numeric(summary_max_n[1, ])
)


ddff_plot <- ddff %>%
  mutate(
    est_ss_plot = ifelse(est_ss > 1500, 1500, est_ss),
    label = case_when(
      is.na(est_ss) ~ "Not reached",
      est_ss > 1500 ~ "> 1500",
      TRUE ~ round(est_ss, 0) |> as.character()
    ),
    Datasets = fct_reorder(Datasets, est_ss, .na_rm = FALSE)
  )


ss_plt <- ggplot(ddff_plot, aes(x = est_ss_plot, y = Datasets)) +
  geom_point(size = 3, na.rm = TRUE) +
  geom_vline(
    data = summary_max_nn,
    aes(xintercept = n, linetype = stat, color = stat),
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = label),
    hjust = -0.24,
    size = 3,
    na.rm = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 1700),
    breaks = c(17, 47, 100, 250, 500, 798, 1000, 1500),
    labels = c("17", "47", "100", "250", "500", "798", "1000", ">1500")
  ) +
  labs(
    x = "Estimated sample size per group",
    y = "Dataset",
    colour = "Observed group-size summary",
    linetype = "Observed group-size summary",
    #title = "Estimated sample size required to reach target power"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(size = 11),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.title = element_text(size = 11),
    panel.spacing = unit(0.7, "lines")
  ) 


ggsave(
  filename = "figures/ss_plt.eps",
  plot = ss_plt,
  device = cairo_ps,   
  width = 13,
  height = 5,
  units = "in"
)

############################################################
nsim = 5; notu = 1000
nsample_vec    =  c(30, 70, 150) #seq(10,150,20)
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
                         logsample_size = log2(sample_size))
  
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
if (!dir.exists("results_other_data/power_sim/")) {
  dir.create("results_other_data/power_sim/", recursive = TRUE)
}
saveRDS(power_dd_list, file = "results_other_data/power_sim/power_dd_list.rds")
saveRDS(pred_dd_list, file = "results_other_data/power_sim/pred_dd_list.rds")
saveRDS(logfoldchange_sim_list, file = "results_other_data/power_sim/logfoldchange_sim_list.rds")
saveRDS(logmean_sim_list, file = "results_other_data/power_sim/logmean_sim_list.rds")
################################################################
power_dd_list = readRDS("results_other_data/power_sim/power_dd_list.rds")
pred_dd_list = readRDS("results_other_data/power_sim/pred_dd_list.rds")
logfoldchange_sim_list = readRDS("results_other_data/power_sim/logfoldchange_sim_list.rds")
logmean_sim_list = readRDS("results_other_data/power_sim/logmean_sim_list.rds")
################################################################
subb = c("Blueberry", "glass_plastic_oberbeckmann", 
         "ob_ross",  "MALL", "Exercise")

sub_samples = c(30, 70, 150)
pred_df <- do.call(rbind, lapply(pred_dd_list[subb], function(df) {
  df[df$sample_size %in% sub_samples, ]
}))


power_df <- do.call(rbind, lapply(power_dd_list[subb], function(df) {
  df[df$sample_size %in% sub_samples, ]
}))


pred_df <- pred_df %>%
  mutate(
    dataset = recode(
      dataset,
      "Blueberry" = "blueberry soil",
      "glass_plastic_oberbeckmann" = "glass-plastic",
      "ob_ross" = "obesity",
      "MALL" = "mall",
      "Exercise" = "exercise"
    )
  )

power_df <- power_df %>%
  mutate(
    dataset = recode(
      dataset,
      "Blueberry" = "blueberry soil",
      "glass_plastic_oberbeckmann" = "glass-plastic",
      "ob_ross" = "obesity",
      "MALL" = "mall",
      "Exercise" = "exercise"
    )
  )
#############################################################################
cont_breaks <- c(0.2, 0.6,0.8) #seq(0.2, 0.8, 0.2)
ddplt  = plot_fun2(pred_df,
                   power_df,
                   cont_breaks, 
                   facet_labs)


ddplt

ggsave(
  filename = "figures/power_contour.eps",
  plot = ddplt,
  device = cairo_ps,   
  width = 11,
  height = 8,
  units = "in"
)



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





