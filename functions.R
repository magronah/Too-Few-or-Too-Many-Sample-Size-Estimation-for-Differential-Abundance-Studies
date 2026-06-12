# Common theme
my_theme <- function(
    base_size = 20,
    axis_title_size = base_size + 2,
    x_text_angle = 60,
    tick_length = 0.25,
    axis_line_width = 1,
    grid = FALSE
) {
  
  theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(
        angle = x_text_angle,
        hjust = 1,
        colour = "black"
      ),
      axis.text.y = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        size = axis_title_size,
        colour = "black"
      ),
      
      axis.ticks.length = unit(tick_length, "cm"),
      
      axis.line = element_line(
        linewidth = axis_line_width
      ),
      
      panel.grid = if (grid) {
        element_line(colour = "grey85")
      } else {
        element_blank()
      },
      
      panel.grid.minor = element_blank(),
  

      axis.ticks = element_line(
        linewidth = axis_line_width,
        colour = "black"
      )
    )
}
########################################################
plot_quantiles <- function(
    df,
    x_labels,
    y_lab,
    box_width = 0.45,
    box_fill = "grey85",
    box_colour = "black",
    whisker_linewidth = 1,
    box_linewidth = 1,
    point_size = 3,
    point_alpha = 1,
    point_colour = "black",
    point_cex = 1.7,
    point_fill = "white",
    x_text_size = 20,
    y_text_size = 20,
    base_size = 20,
    axis_title_size = base_size + 2,
    axis_line_width = 0.8,
    axis_tick_width = 0.8,
    tick_length = 0.2,
    margins = c(10, 10, 10, 10),
    panel_label = NULL,
    panel_x = -Inf,
    panel_y = Inf,
    panel_hjust = -0.6,
    panel_vjust = 1.4,
    panel_size = 6
) {
  
 p <- ggplot(df, aes(x = quantile, y = value)) +
    geom_boxplot(
      width = box_width,
      fill = box_fill,
      whisker.linewidth = whisker_linewidth,
      box.linewidth   =  box_linewidth,
      colour = box_colour
    ) +
    geom_beeswarm(
      size = point_size,
      alpha = point_alpha,
      colour = point_colour,
      cex = point_cex,
      fill = point_fill
    ) +
    scale_x_discrete(labels = x_labels) +
    labs(
      x = NULL,
      y = y_lab
    ) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      
      axis.text.x = element_text(
        size = x_text_size,
        colour = "black"
      ),
      axis.text.y = element_text(
        size = y_text_size,
        colour = "black"
      ),
      
      axis.title.x = element_text(
        size = axis_title_size,
        colour = "black"
      ),
      axis.title.y = element_text(
        size = axis_title_size,
        colour = "black"
      ),
      
      axis.line = element_line(
        linewidth = axis_line_width,
        colour = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = axis_tick_width,
        colour = "black"
      ),
      
      axis.ticks.length = unit(
        tick_length,
        "cm"
      ),
      
      plot.margin = margin(
        margins[1],
        margins[2],
        margins[3],
        margins[4]
      )
    )
 
 if (!is.null(panel_label)) {
   p <- p +
     annotate(
       "text",
       x = panel_x,
       y = panel_y,
       label = panel_label,
       hjust = panel_hjust,
       vjust = panel_vjust,
       fontface = "bold",
       size = panel_size
     )
 }
 
 p
 
}

########################################################
########################################################
# Helper to read + clean names + convert to df
read_to_df <- function(pattern, prefix_remove, abs_val = FALSE) {
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  names <- gsub(paste0("^", prefix_remove, "|\\.rds$"), "", basename(files))
  
  map_dfr(seq_along(files), function(i) {
    vals <- as.numeric(unlist(readRDS(files[i])))
    data.frame(Dataset = names[i], value = vals)
  })
}

########################################################
########################################################
read_to_list <- function(pattern, prefix_remove) {
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  
  out <- lapply(files, readRDS)
  
  names(out) <- gsub(
    paste0("^", prefix_remove, "|\\.rds$"),
    "",
    basename(files)
  )
  
  out
}

########################################################
########################################################
# reorder by median
reorder_by_median <- function(df, magnitude_only = FALSE) {
  
  if(magnitude_only == TRUE){
    ddf <- df %>% 
          group_by(Dataset) %>%
          mutate(med = median(abs(value))) %>%
          ungroup() %>%
          mutate(Dataset = reorder(Dataset, med))
  return(ddf)
  }else{
   ddf <- df %>%
      group_by(Dataset) %>%
      mutate(med = median(value)) %>%
      ungroup() %>%
      mutate(Dataset = reorder(Dataset, med))
   return(ddf)
   
  }
  
  
  
}

########################################################
########################################################
# ss_var_lmc <- function(target_powers,lmc_vals, mod){
#   dd_pow <- do.call(rbind, lapply(lmc_vals, function(lmc) {
#     
#     ss_vals <- sapply(target_powers, function(tp) {
#       ss_solver(
#         target_power = tp,
#         logmean = lmc,
#         abs_lfc = 1.5,
#         model = mod,
#         xmin = log2(5),
#         xmax = log2(500)
#       )
#     })
#     
#     data.frame(
#       sample_size = ss_vals,
#       power = target_powers,
#       logmean = lmc
#     )
#   }))
#   
# }
# #################################################################
facet_labs <- function(x) {
  paste0("n = ", x, " samples per group")
}

#####################################################################
plot_fun <- function(
    pred_df,
    power_df,
    cont_breaks,
    facet_labs,
    contour_args = list(
      linewidth = 1,
      colour = "blue"
    ),
    point_args = list(
      alpha = 0.4,
      colours = c("grey", "red"),
      labels = c("No", "Yes"),
      legend_name = "Significant Taxa"
    ),
    contour_label_args = list(
      digits = 2,
      size = 4,
      label_size = 0.25,
      fill = "white",
      colour = "black"
    ),
    hseg_args = list(
      x = 0.8,
      xend = 3.2,
      y = 0.0,
      yend = 0.0,
      arrow_length = 0.2,
      ends = "last",
      type = "closed",
      linewidth = 0.8,
      colour = "black"
    ),
    htext_args = list(
      x = 0.7,
      y = 0.33,
      label = "higher power",
      size = 4
    ),
    vseg_args = list(
      x = -2,
      xend = -2,
      y = 1.00,
      yend = 1.80,
      arrow_length = 0.2,
      ends = "last",
      type = "closed",
      linewidth = 0.8,
      colour = "black"
    ),
    vtext_args = list(
      x = -2.70,
      y = 2.00,
      label = "higher power",
      angle = 90,
      size = 4
    ),
    theme_args = list(
      base_size = 12,
      strip_text_size = 12,
      legend_title_size = 12,
      legend_text_size = 12,
      axis_text_x_size = 12,
      axis_text_y_size = 12
    ) 
) {
  
  gg_2dimc <- ggplot(pred_df) +
    aes(logmean, abs_lfc, group = sample_size) +
    
    ggrastr::rasterise(
      geom_point(
        aes(color = pvalue_reject),
        alpha = point_args$alpha
      )
    ) +
    scale_colour_manual(
      values = point_args$colours,
      labels = point_args$labels,
      name = point_args$legend_name
    ) +     
    ggplot2::geom_contour(
      data = power_df,
      aes(x = logmean, y = abs_lfc, z = .data$power),
      linewidth = contour_args$linewidth,
      breaks = cont_breaks,
      colour = contour_args$colour
    ) +
    
    metR::geom_label_contour(
      data = power_df,
      aes(
        x = logmean,
        y = abs_lfc,
        z = .data$power,
        label = sprintf(
          paste0("%.", contour_label_args$digits, "f"),
          after_stat(level)
        )
      ),
      breaks = cont_breaks,
      size = contour_label_args$size,
      label.size = contour_label_args$label_size,
      fill = contour_label_args$fill,
      colour = contour_label_args$colour
    ) +
    
    annotate(
      "segment",
      x = hseg_args$x,
      xend = hseg_args$xend,
      y = hseg_args$y,
      yend = hseg_args$yend,
      arrow = arrow(
        length = unit(hseg_args$arrow_length, "cm"),
        ends = hseg_args$ends,
        type = hseg_args$type
      ),
      linewidth = hseg_args$linewidth,
      colour = hseg_args$colour
    ) +
    
    annotate(
      "text",
      x = htext_args$x,
      y = htext_args$y,
      label = htext_args$label,
      size = htext_args$size
    ) +
    
    annotate(
      "segment",
      x = vseg_args$x,
      xend = vseg_args$xend,
      y = vseg_args$y,
      yend = vseg_args$yend,
      arrow = arrow(
        length = unit(vseg_args$arrow_length, "cm"),
        ends = vseg_args$ends,
        type = vseg_args$type
      ),
      linewidth = vseg_args$linewidth,
      colour = vseg_args$colour
    ) +
    
    annotate(
      "text",
      x = vtext_args$x,
      y = vtext_args$y,
      label = vtext_args$label,
      angle = vtext_args$angle,
      size = vtext_args$size
    ) +
    
    ggplot2::theme_bw(base_size = theme_args$base_size) +
    facet_grid(dataset ~ sample_size,
               labeller = labeller(sample_size = facet_labs),
               scales = "free") +
    theme(
      strip.text = element_text(size = theme_args$strip_text_size, color = "black"),
      legend.title = element_text(size = theme_args$legend_title_size, color = "black"),
      legend.text = element_text(size = theme_args$legend_text_size, color = "black"),
      axis.text.x = element_text(size = theme_args$axis_text_x_size, color = "black"),
      axis.text.y = element_text(size = theme_args$axis_text_y_size, color = "black")
    ) +
    xlab(TeX("$\\log_2$(mean counts)")) +
    ylab(TeX("|$\\log_2$(fold change)|")) 
    
  
  gg_2dimc
}
################################################################
###############################################################
plot_power_curve <- function(data,
                             group_var,
                             legend_label,
                             nn = 15,
                             okabe_ito = c(
                               "#000000", "#E69F00", "#56B4E9", "#009E73",
                               "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
                             ),
                             shape_values = c(16, 17, 15, 18, 8, 3, 7, 1),
                             actual_n = NULL) {
  
  group_var <- enquo(group_var)
  
  p <- ggplot(
    data,
    aes(
      x = sample_size,
      y = power,
      color = factor(!!group_var),
      shape = factor(!!group_var),
      group = factor(!!group_var)
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(
      yintercept = 0.8,
      linetype = "dashed",
      linewidth = 0.8,
      color = "grey40"
    ) +
    scale_color_manual(values = okabe_ito) +
    scale_shape_manual(values = shape_values) +
    scale_x_continuous(trans = "log2",
                       breaks = function(x) {
                         rng <- range(x, na.rm = TRUE)
                         2^seq(floor(log2(rng[1])), ceiling(log2(rng[2])))
                       },
                       labels = scales::label_number()) +
    labs(
      x = "Sample size per group",
      y = "Target Power",
      color = legend_label,
      shape = legend_label
    ) +
    theme_bw() +
    theme(
      legend.title = element_text(size = nn, color = "black"),
      legend.text = element_text(size = nn, color = "black"),
      axis.text.x = element_text(size = nn, color = "black"),
      axis.text.y = element_text(size = nn, color = "black"),
      axis.title.x = element_text(size = nn, color = "black"),
      axis.title.y = element_text(size = nn, color = "black"),
      strip.text = element_text(size = nn, color = "black"),
      panel.spacing = unit(0.4, "lines")
    ) +
    facet_wrap(~dataset, scales = "free_x")
  
  if (!is.null(actual_n)) {
    p <- p +
      geom_vline(
        data = actual_n,
        aes(xintercept = max_n),
        inherit.aes = FALSE,
        linetype = "dashed",
        linewidth = 0.8,
        color = "grey20"
      )
  }
  
  p
}
################################################################
###############################################################
# predict_power_at_n <- function(log_n, logmean, abs_lfc, model) {
#   
#   newdat <- data.frame(
#     logsample_size = log_n,
#     logmean = logmean,
#     abs_lfc = abs_lfc
#   )
#   
#   pred <- predict(
#     model,
#     newdata = newdat[rep(1, 2), ],
#     type = "response"
#   )
#   
#   pred[1]
# }


################################################################
###############################################################

# ss_solver <- function(target_power, logmean, abs_lfc, model,
#                       xmin = log2(5), xmax = log2(500)) {
#   
#   if(target_power == 0 || target_power == 1){
#     warning("statistical power for  0% or 100% is  unlikely")
#   }
#   out <- tryCatch({
#     
#     uniroot_ss(
#       target_power = target_power,
#       logmean = logmean,
#       abs_lfc = abs_lfc,
#       model = model,
#       xmin = xmin,
#       xmax = xmax
#     )
#     
#   }, error = function(e) {
#     
#     sample_size_ss_interp(
#       target_power = target_power,
#       logmean = logmean,
#       abs_lfc = abs_lfc,
#       model = model,
#       xmin = xmin,
#       xmax = xmax
#     )
#     
#   })
#   
#   out
# }

############################################################################
power_grid <- function(pred_data,model,grid_len){
  df <- with(pred_data,
             expand.grid(logmean = seq(min(logmean),
                                       max(logmean),
                                       length  = grid_len),
                         abs_lfc   =  seq(min(abs_lfc),
                                          max(abs_lfc),
                                          length  =  grid_len),
                         sample_size = nsample_vec))
  df$logsample_size = log2(df$sample_size)
  df$power   =   predict(model, newdata = df,type = "response")
  df$sample_size   <- factor(df$sample_size) 
  df
}
############################################################################
# sample_size_grid <- function(target_powers, vary_vals, mod,
#                              vary = c("abs_lfc", "logmean"),
#                              fixed_logmean = NULL,
#                              fixed_abs_lfc = NULL,
#                              xmin = log2(5),
#                              xmax = log2(500)) {
#   
#   vary <- match.arg(vary)
#   
#   do.call(rbind, lapply(vary_vals, function(val) {
#     
#     ss_vals <- sapply(target_powers, function(tp) {
#       
#       logmean_val <- if (vary == "logmean") val else fixed_logmean
#       abs_lfc_val <- if (vary == "abs_lfc") val else fixed_abs_lfc
#       
#       ss_solver(
#         target_power = tp,
#         logmean = logmean_val,
#         abs_lfc = abs_lfc_val,
#         model = mod,
#         xmin = xmin,
#         xmax = xmax
#       )
#     })
#     
#     out <- data.frame(
#       sample_size = ss_vals,
#       power = target_powers
#     )
#     
#     out[[vary]] <- val
#     out
#   }))
# }


# sample_size_est <- function(target_powers,lfc_vals, mod, logmean = 5, 
#                    ssmin = 5, ssmax = 500){
#   
#   ss <- do.call(rbind, lapply(lfc_vals, function(lfc) {
#     
#     ss_vals <- sapply(target_powers, function(tp) {
#       ss_solver(
#         target_power = tp,
#         logmean =  logmean,
#         abs_lfc = lfc,
#         model = mod,
#         xmin = log2(ssmin),
#         xmax = log2(ssmax)
#       )
#     })
#     
#     data.frame(
#       sample_size = ss_vals,
#       power = target_powers,
#       abs_lfc = lfc
#     )
#   }))
#   
#   ss
#   
# }




