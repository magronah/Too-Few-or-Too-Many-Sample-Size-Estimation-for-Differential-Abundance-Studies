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
library(ggbeeswarm)
#############################################################
##############################################################
source("functions.R")
source("param_vals.R")
path <- file.path(getwd(), "results_other_data")
file_path <-  "figures/"
##############################################################
logmean_df  <-  read_to_df("^logmeanEst_", "logmeanEst_")
logfc_df    <-  read_to_df("^logfoldchangeEst_", "logfoldchangeEst_")

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

lmc_dist <- ggplot(logmean_df, aes(Dataset, value)) +
  geom_violin(fill = "grey85",
              colour = "black",
              alpha = 0.5,
              trim = FALSE) +
  geom_boxplot(width = 0.1,
               fill = "white",
               outlier.size = 0.8) +
  labs(x = "Dataset",
       y = expression(log[2]("mean abundance"))) +
  my_theme(base_size = 12)


abslfc_dist <- ggplot(logfc_df, aes(Dataset, abs(value))) +
  geom_violin(fill = "grey85",
              colour = "black",
              alpha = 0.5,
              trim = FALSE) +
  geom_boxplot(width = 0.1,
               fill = "white",
               outlier.size = 0.8) +
  labs(x = "Dataset",
       y = expression("|" * log[2]("fold change") * "|")) +
  scale_y_continuous() +
  my_theme(base_size = 12)

############################################################
############################################################
ggsave(
  filename = paste0(file_path,"lmc_dist.eps"),
  plot = lmc_dist,
  device = cairo_ps,  
  width = 10,
  height = 8,
  units = "in"
)

ggsave(
  filename = paste0(file_path,"abslfc_dist.eps"), 
  plot = abslfc_dist,
  device = cairo_ps,  
  width = 10,
  height = 8,
  units = "in"
)

############################################################
############################################################
### Compute quantiles
logmean_quantiles <- logmean_df %>%
  group_by(Dataset) %>%
  summarise(
    q10 = quantile(value, 0.1, na.rm = TRUE),
    q50 = quantile(value, 0.5, na.rm = TRUE),
    q90 = quantile(value, 0.9, na.rm = TRUE),
    .groups = "drop"
  )


long_logmean_quant <- logmean_quantiles %>%
  pivot_longer(
    cols = c(q10, q50, q90),
    names_to = "quantile",
    values_to = "value")
############################################################
logfc_quantiles <- logfc_df %>%
  group_by(Dataset) %>%
  summarise(
    q10 = quantile(abs(value), 0.1, na.rm = TRUE),
    q50 = quantile(abs(value), 0.5, na.rm = TRUE),
    q90 = quantile(abs(value), 0.9, na.rm = TRUE),
    .groups = "drop"
  )

long_logfc_quant <- logfc_quantiles %>%
  pivot_longer(
    cols = c(q10, q50, q90),
    names_to = "quantile",
    values_to = "value")


pp_logmean <- plot_quantiles(
  df = long_logmean_quant,
  x_labels = c(
    q10 = "Low abundance\n(10th percentile)",
    q50 = "Medium abundance\n(50th percentile)",
    q90 = "High abundance\n(90th percentile)"
  ),
  y_lab = expression(log[2]*"(mean abundance)"),
  panel_label = "(A)"
)

pp_logfc <- plot_quantiles(
  df = long_logfc_quant,
  x_labels = c(
    q10 = "Small fold change \n(10th percentile)",
    q50 = "Medium fold change \n(50th percentile)",
    q90 = "Large fold change \n(90th percentile)"
  ),
  y_lab = expression(abs(log[2]*"(fold change)")),
  panel_label = "(B)"
)


summary_abslfc =  c(small   =  median(logfc_quantiles$q10), 
                    medium  =  median(logfc_quantiles$q50),
                    large   =  median(logfc_quantiles$q90))

summary_lmc =  c(low     =  median(logmean_quantiles$q10),
                 medium  =  median(logmean_quantiles$q50),
                 large   =  median(logmean_quantiles$q90))
############################################################
############################################################
ggsave(
  filename = paste0(file_path,"lmc_quant.eps"), 
  plot = pp_logmean,
  device = cairo_ps,  
  width = 10,
  height = 8,
  units = "in"
)

ggsave(
  filename = paste0(file_path,"abslfc_quant.eps"), 
  plot = pp_logfc,
  device = cairo_ps,  
  width = 10,
  height = 8,
  units = "in"
)

############################################################
########################################################
### plot sample sizes per group
control_n <- mapply(function(x, ref) {
  if (ref %in% names(x)) {
    return(x[ref])
  } else {
    return(NA)  
  }
}, sample_size_list, reference_names)

treatment_n <- mapply(function(x, ref) {
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


sample_size_df <- sample_size_df |>
  arrange(min_n) |>
  mutate(rank = row_number())
############################################################
############################################################
pp_sample_size <- ggplot(sample_size_df) +
  geom_segment(aes(x=rank, y=min_n, xend = rank, yend = max_n)) +
  geom_point(aes(x=rank, y = min_n)) +
  geom_point(aes(x=rank, y = max_n)) +
  labs(
    x = "Dataset (ordered by minimum group size)",
    y = "Group sample size",
  ) +
  my_theme() +
  scale_y_log10(
    breaks = c(5, 10, 20, 50, 100, 200, 500, 800)
  )

ggsave(
  filename = paste0(file_path,"sample_sizes.eps"),
  plot = pp_sample_size,
  device = cairo_ps,   
  width = 10,
  height = 8,
  units = "in"
)

############################################################
############################################################
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
  filename = paste0(file_path,"cont_treat_ss.eps"),
  plot = pp_ss,
  device = cairo_ps,   
  width = 10,
  height = 8,
  units = "in"
)

############################################################
############################################################
summary_sample_size  <- quantile(sample_size_df$max_n, 
                                probs = c(0.1, 0.5, 0.9))


summary_datasets <- expand.grid(
  mean_level = c("Low mean count (-1.62)", "Medium mean count (0.09)", 
                 "High mean count (3.34)"),
  foldchange  = c("Small (0.09)", "Medium (0.49)", "Large (1.32)"),
  n_level    = c("Small sample size (n=20)", "Medium sample size (n=45)",
                 "Large sample size (n=181)"))


summary_datasets$logmean <- rep(as.numeric(summary_lmc), times = 9)
summary_datasets$abs_lfc <- rep(rep(as.numeric(summary_abslfc), each = 3), 
                                times = 3)
summary_datasets$sample_size <- rep(as.numeric(summary_sample_size), each = 9)
summary_datasets$logsample_size <- log2(summary_datasets$sample_size)

############################################################
############################################################

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
  geom_jitter(
  width = 0.15,
    height = 0,
    size = 1.8,
    alpha = 0.6,
    colour = "black"
  ) +
  scale_y_continuous() +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(n_level ~ mean_level) +
  labs(
    x = expression("|" * log[2]("fold change") * "|"),
    y = "Predicted power",
    fill = expression("|" * log[2]("fold change") * "|")
    ) +
  my_theme(base_size = 18,
           grid = TRUE, 
           axis_title_size = 18)

############################################################
############################################################
ggsave(
  filename =  paste0(file_path,"power_plt.eps"),
  plot = pp_power,
  device = cairo_ps,   
  width = 15,
  height = 11,
  units = "in"
)

############################################################
############################################################
summary_datasets2 <- expand.grid(
  mean_level = c("Low mean count (-2.34)", "Median mean count (-0.02)", 
                 "High mean count (2.20)"),
  foldchange  = c("Low (0.00)", "Median (0.05)", "High (0.66)"),
  power_level    = c("Power = 80%", "Power = 85%",
                     "Power = 90%"))

power_vec  <-  c(0.8, 0.85, 0.9)
summary_datasets2$logmean <- rep(as.numeric(summary_lmc), times = 9)
summary_datasets2$abs_lfc <- rep(rep(as.numeric(summary_abslfc), each = 3), times = 3)
summary_datasets2$target_power <- rep(as.numeric(power_vec), each = 9)

############################################################
############################################################
target_powers = c(0.8)

res_list = list()

for(i in 1:length(data_names)){
  model = gam_fit_list[[i]]$gam_mod
  res = power.nb:::uniroot_ss(target_power=0.8,logmean=max(as.numeric(summary_lmc)), 
                   abs_lfc=max(as.numeric(summary_abslfc)),
                   model,xmin = log2(10),
                   xmax = log2(5000),maxiter = 2000)
  
  res_list[[i]] = res
}


############################################################
############################################################
ss_est  =  read_data(res_list,"sample_size_per_group")
ddf = data.frame(Dataset = data_names,
                  est_ss = unlist(ss_est))


ddff  <-  merge(ddf, sample_size_df, by  = "Dataset")
ddff  <-  ddff[complete.cases(ddff), ]

ddff_plot <- ddff %>%
  mutate(
    est_ss_plot = ifelse(est_ss > 1500, 1500, est_ss),
    label = case_when(
      is.na(est_ss) ~ "",
      est_ss > 1500 ~ "> 1500",
      TRUE ~ as.character(round(est_ss, 0))
    )
  ) %>%
  arrange(est_ss_plot) %>%
  mutate(Dataset = factor(Dataset, levels = Dataset))

summary_ss_dd <-  data.frame(n  =  as.numeric(round(summary_sample_size)),
                           stat =  paste0(names(round(summary_sample_size)),
                                          " ", "percentile"))


ss_plt <- ggplot(ddff_plot, aes(y = Dataset)) +
  geom_segment(
    aes(x = min_n, xend = max_n, yend = Dataset),
    linewidth = 1.6,
    colour = "grey65",
    alpha = 0.9,
    lineend = "round"
  ) +
  geom_point(
    aes(x = est_ss_plot),
    size = 3,
    colour = "black",
    na.rm = TRUE
  ) +
  geom_text(
    aes(x = est_ss_plot, label = label),
    hjust = -0.34,
    size = 4,
    na.rm = TRUE
  ) +
  geom_vline(
    data = summary_ss_dd,
    aes(xintercept = n, linetype = stat, color = stat),
    linewidth = 0.8
  ) +
  scale_x_continuous(
    limits = c(0, 1700),
    breaks = c(20, 45, 100, 181, 250, 500, 1000, 1500),
    labels = c("20", "45", "100", "181", "250", "500", "1000", ">1500")
  ) +
  labs(
    x = "Sample size per group",
    y = "Dataset",
    colour = "Group-size summary",
    linetype = "Group-size summary"
  ) +
  my_theme(base_size = 13, 
           x_text_angle = 0,
           grid = TRUE)


ggsave(
  filename = paste0(file_path,"ss_plt.eps"), 
  plot = ss_plt,
  device = cairo_ps,   
  width = 17,
  height = 9,
  units = "in"
)

############################################################
############################################################
nsim = 5; notu = 1000 
nsample_vec    =  seq(10,150,20)
grid_len =  20  
############################################################
power_dd_list  =  pred_dd_list = list()
logfoldchange_sim_list  =  logmean_sim_list  = list()

for(i in 1:length(data_names)){

  dispersion_param      =  dispersionFit_list[[i]]$param
  logmean_param         =  logmeanFit_list[[i]]$param
  logfoldchange_param   =  logfoldchangeFit_list[[i]]
  ############################################################
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
  ############################################################
  logfoldchange_list  =   read_data(countdata_sims_list,"logfoldchange_list")
  logmean_list        =   read_data(countdata_sims_list,"logmean_list")
  logfoldchange_sim_list[[i]]  =  logfoldchange_list
  logmean_sim_list[[i]]        =  logmean_list
  ############################################################
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
  ############################################################
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
  ############################################################
  sample_size         =   rep(nsample_vec,
                              times = sapply(lapply(logfoldchange_list, unlist),
                                             length))
  ############################################################
  pred_data = data.frame(logmean  =  unlist(logmean_list),
                         abs_lfc  =  abs(unlist(logfoldchange_list)),
                         sample_size   =  sample_size,
                         logsample_size = log2(sample_size))
  
  pred_data$pvalue_reject    =   factor(as.numeric(pval_reject))
  rownames(pred_data)        =   NULL
  ############################################################
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
############################################################
if (!dir.exists("results_other_data/power_sim/")) {
  dir.create("results_other_data/power_sim/", recursive = TRUE)
}
saveRDS(power_dd_list, file = "results_other_data/power_sim/power_dd_list.rds")
saveRDS(pred_dd_list, file = "results_other_data/power_sim/pred_dd_list.rds")
saveRDS(logfoldchange_sim_list, file = "results_other_data/power_sim/logfoldchange_sim_list.rds")
saveRDS(logmean_sim_list, file = "results_other_data/power_sim/logmean_sim_list.rds")
################################################################
# power_dd_list = readRDS("results_other_data/power_sim/power_dd_list.rds")
# pred_dd_list = readRDS("results_other_data/power_sim/pred_dd_list.rds")
# logfoldchange_sim_list = readRDS("results_other_data/power_sim/logfoldchange_sim_list.rds")
# logmean_sim_list = readRDS("results_other_data/power_sim/logmean_sim_list.rds")

################################################################
############################################################
subb = c("Blueberry", "glass_plastic_oberbeckmann", 
         "ob_ross",  "MALL", "Exercise")

sub_samples = c(30, 70, 150)
pred_df <- do.call(rbind, lapply(pred_dd_list[subb], function(df) {
  df[df$sample_size %in% sub_samples, ]
}))


power_df <- do.call(rbind, lapply(power_dd_list[subb], function(df) {
  df[df$sample_size %in% sub_samples, ]
}))

################################################################
############################################################
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
################################################################
############################################################
cont_breaks <- c(0.2, 0.6,0.8) 
ddplt  = plot_fun(pred_df,
                   power_df,
                   cont_breaks, 
                   facet_labs)


ddplt

ggsave(
  filename = paste0(file_path,"power_contour.eps"), 
  plot = ddplt,
  device = cairo_ps,   
  width = 11,
  height = 8,
  units = "in"
)
#############################################################################
##Contour plots for other datasets
## in sets of 5 and then 5 set

remaining_list  =  list()
dataset_names <- setdiff(names(sample_size_list), subb)
splits  <- split(dataset_names, ceiling(seq_along(dataset_names) / 5))
cont_breaks <- c(0.1,0.2, 0.6,0.8) 
plot_names  = paste0("cont_plot",1:length(splits),".eps")

for(i in 1:length(splits)){
  pred_df <- do.call(rbind, lapply(pred_dd_list[splits[[i]]], function(df) {
    df[df$sample_size %in% sub_samples, ]
  }))
  
  power_df <- do.call(rbind, lapply(power_dd_list[splits[[i]]], function(df) {
    df[df$sample_size %in% sub_samples, ]
  }))
  
  ddplt  = plot_fun(pred_df,
                     power_df,
                     cont_breaks, 
                     facet_labs)
  
  ggsave(
    filename = paste0(file_path,plot_names[i]),
    plot = ddplt,
    device = cairo_ps,   
    width = 12,
    height = 10,
    units = "in"
  )
  
  remaining_list[[i]] = ddplt
}

names(remaining_list) = plot_names

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





