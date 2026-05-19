# Common theme
my_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 12, colour = "black"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 16, colour = "black"),
    axis.title.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 16, colour = "black"),
    axis.ticks.length = unit(0.25, "cm"),   # increase tick length
    
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_line(linewidth = 1, colour = "black")
  )



# Helper to read + clean names + convert to df
read_to_df <- function(pattern, prefix_remove, abs_val = FALSE) {
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  names <- gsub(paste0("^", prefix_remove, "|\\.rds$"), "", basename(files))
  
  map_dfr(seq_along(files), function(i) {
    vals <- as.numeric(unlist(readRDS(files[i])))
    data.frame(Dataset = names[i], value = vals)
  })
}


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

# reorder by median
reorder_by_median <- function(df) {
  df %>%
    group_by(Dataset) %>%
    mutate(med = median(value)) %>%
    ungroup() %>%
    mutate(Dataset = reorder(Dataset, med))
}




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
#####################################################################
facet_labs <- function(x) {
  paste0("n = ", x, " samples per group")
}

plot_fun2 <- function(
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
###############################################################################
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
  
  # if (add_band) {
  # band_df <- data %>%
  #     group_by(dataset) %>%
  #     summarise(
  #       lower_95 = quantile(sample_size, 0.025, na.rm = TRUE),
  #       upper_95 = quantile(sample_size, 0.975, na.rm = TRUE))
  # 
  # p <- p +
  #   geom_rect(
  #     data = band_df,
  #     aes(xmin = lower_95, xmax = upper_95, ymin = -Inf, ymax = Inf),
  #     inherit.aes = FALSE,
  #     fill = "grey80",
  #     alpha = 0.3
  #   )
  # }
  
  p
}
###############################################################################
predict_power_at_n <- function(log_n, logmean, abs_lfc, model) {
  
  newdat <- data.frame(
    logsample_size = log_n,
    logmean = logmean,
    abs_lfc = abs_lfc
  )
  
  pred <- predict(
    model,
    newdata = newdat[rep(1, 2), ],
    type = "response"
  )
  
  pred[1]
}


uniroot_ss <- function(target_power, logmean, abs_lfc, model,
                       xmin = log2(10), xmax = log2(500),
                       maxiter = 2000) {
  
  f <- function(log_n) {
    predict_power_at_n(
      log_n = log_n,
      logmean = logmean,
      abs_lfc = abs_lfc,
      model = model
    ) - target_power
  }
  
  f_min <- f(xmin)
  f_max <- f(xmax)
  
  if (is.na(f_min) || is.na(f_max)) {
    return(NA_real_)
  }
  
  if (f_min >= 0) {
    return(2^xmin)
  }
  
  if (f_max < 0) {
    return(NA_real_)
  }
  
  root <- uniroot(
    f,
    interval = c(xmin, xmax),
    maxiter = maxiter, 
    extendInt = "yes"
   )$root
  
  2^root
}




uniroot_ss2 =  function(target_power,logmean, abs_lfc,model,xmin,xmax,
                        maxiter = 10000){

  root <- uniroot(function(ss) {

    data = data.frame(logsample_size = ss,
                      logmean = logmean,
                      abs_lfc = abs_lfc)

    pred = predict(model,
                   type = "response",
                   newdata = data[rep(1,2),])

    pred[[1]] - target_power
  },
  interval = c(xmin, xmax),
  extendInt = "yes",maxiter = maxiter)$root
  2^root
  
}



uniroot_ss3 <- function(target_power, logmean, abs_lfc, model,
                        xmin, xmax, maxiter = 10000,
                        max_report = 2000) {
  
  f <- function(ss) {
    
    data <- data.frame(
      logsample_size = ss,
      logmean = logmean,
      abs_lfc = abs_lfc
    )
    
    pred <- predict(
      model,
      type = "response",
      newdata = data[rep(1, 2), ]
    )
    
    pred[[1]] - target_power
  }
  
  f_min <- f(xmin)
  f_max <- f(xmax)
  
  pred_min <- f_min + target_power
  pred_max <- f_max + target_power
  
  if (is.na(f_min) || is.na(f_max)) {
    return(list(
      sample_size = NA,
      conclusion = "Sample size could not be estimated because the model returned NA predictions.",
      predicted_power_min_n = pred_min,
      predicted_power_max_n = pred_max
    ))
  }
  
  if (f_min >= 0) {
    ss_est <- 2^xmin
    
    return(list(
      sample_size = ss_est,
      conclusion = "Target power is already achieved at the minimum sample size.",
      predicted_power_min_n = pred_min,
      predicted_power_max_n = pred_max
    ))
  }
  
  if (f_max < 0) {
    return(list(
      sample_size = NA,
      conclusion = paste0(
        "No sign change was found. The predicted power does not reach the target power ",
        "within the specified sample-size interval."
      ),
      predicted_power_min_n = pred_min,
      predicted_power_max_n = pred_max
    ))
  }
  
  root <- tryCatch(
    uniroot(
      f,
      interval = c(xmin, xmax),
      maxiter = maxiter
    )$root,
    error = function(e) NA_real_
  )
  
  if (is.na(root)) {
    return(list(
      sample_size = NA,
      conclusion = "Root finding failed even though the endpoint checks suggested a solution.",
      predicted_power_min_n = pred_min,
      predicted_power_max_n = pred_max
    ))
  }
  
  ss_est <- 2^root
  
  if (ss_est > max_report) {
    return(list(
      sample_size = paste0("> ", max_report),
      conclusion = paste0(
        "The estimated sample size exceeds ",
        max_report,
        "."
      ),
      predicted_power_min_n = pred_min,
      predicted_power_max_n = pred_max
    ))
  }
  
  list(
    sample_size = ss_est,
    conclusion = "Success.",
    predicted_power_min_n = pred_min,
    predicted_power_max_n = pred_max
  )
}

###############################################################################
sample_size_ss_interp <- function(target_power, logmean, abs_lfc, model,
                                  xmin = log2(5), xmax = log2(500),
                                  ngrid = 1000) {
  
  ss_grid <- seq(xmin, xmax, length.out = ngrid)
  
  nd <- data.frame(
    logsample_size = ss_grid,
    logmean = logmean,
    abs_lfc = abs_lfc
  )
  
  pred <- predict(model, type = "response", newdata = nd)
  
  idx <- which(pred >= target_power)[1]
  
  if (is.na(idx)) return(NA_real_)     # target never reached
  if (idx == 1) return(2^ss_grid[1])   # already achieved at minimum
  
  x0 <- ss_grid[idx - 1]
  x1 <- ss_grid[idx]
  y0 <- pred[idx - 1]
  y1 <- pred[idx]
  
  # linear interpolation on log2(sample size) scale
  x_star <- x0 + (target_power - y0) * (x1 - x0) / (y1 - y0)
  
  2^x_star
}

###############################################################################
ss_solver <- function(target_power, logmean, abs_lfc, model,
                      xmin = log2(5), xmax = log2(500)) {
  
  if(target_power == 0 || target_power == 1){
    warning("statistical power for  0% or 100% is  unlikely")
  }
  out <- tryCatch({
    
    uniroot_ss(
      target_power = target_power,
      logmean = logmean,
      abs_lfc = abs_lfc,
      model = model,
      xmin = xmin,
      xmax = xmax
    )
    
  }, error = function(e) {
    
    sample_size_ss_interp(
      target_power = target_power,
      logmean = logmean,
      abs_lfc = abs_lfc,
      model = model,
      xmin = xmin,
      xmax = xmax
    )
    
  })
  
  out
}
##############################################################################

# plot_fun <- function(dd_sub,pow_sub,cont_breaks,facet_labs,
#                      
#     contour_args = list(linewidth = 1,colour = "blue"),
#     
#     point_args   =  list(alpha = 0.4,colours = c("grey", "red"),
#                          labels = c("No", "Yes"),legend_name = "Significant Taxa"),
#     
#     contour_label_args = list(digits = 2, size = 5, label_size = 0.25,
#                               fill = "white",colour = "black"),
#     hseg_args = list(
#       x = 0.8,
#       xend = 3.2,
#       y = 0.20,
#       yend = 0.20,
#       arrow_length = 0.2,
#       ends = "last",
#       type = "closed",
#       linewidth = 0.8,
#       colour = "black"
#     ),
#     htext_args = list(
#       x = 2.0,
#       y = 0.33,
#       label = "higher power",
#       size = 4.5
#     ),
#     vseg_args = list(
#       x = -2,
#       xend = -2,
#       y = 0.55,
#       yend = 1.80,
#       arrow_length = 0.2,
#       ends = "last",
#       type = "closed",
#       linewidth = 0.8,
#       colour = "black"
#     ),
#     vtext_args = list(
#       x = -2.70,
#       y = 3.20,
#       label = "higher power",
#       angle = 90,
#       size = 4.5
#     ),
#     theme_args = list(
#       base_size = 16,
#       strip_text_size = 15,
#       legend_title_size = 15,
#       legend_text_size = 13,
#       axis_text_x_size = 16,
#       axis_text_y_size = 16
#     )
# ) {
#   
#   
#   gg_2dimc <- ggplot(dd_sub) +
#     aes(logmean, abs_lfc, group = sample_size) +
#     
#     ggrastr::rasterise(
#       geom_point(
#         aes(color = pvalue_reject),
#         alpha = point_args$alpha
#       )
#     ) +
#     
#     xlab(TeX("$\\log_2$(mean counts)")) +
#     ylab(TeX("|$\\log_2$(fold change)|")) +
#     
#     scale_colour_manual(
#       values = point_args$colours,
#       labels = point_args$labels,
#       name = point_args$legend_name
#     ) +
#     
#     ggplot2::geom_contour(
#       data = pow_sub,
#       aes(x = logmean, y = abs_lfc, z = .data$power),
#       linewidth = contour_args$linewidth,
#       breaks = cont_breaks,
#       colour = contour_args$colour
#     ) +
#     
#     metR::geom_label_contour(
#       data = pow_sub,
#       aes(
#         x = logmean,
#         y = abs_lfc,
#         z = .data$power,
#         label = sprintf(
#           paste0("%.", contour_label_args$digits, "f"),
#           after_stat(level)
#         )
#       ),
#       breaks = cont_breaks,
#       size = contour_label_args$size,
#       label.size = contour_label_args$label_size,
#       fill = contour_label_args$fill,
#       colour = contour_label_args$colour
#     ) +
#     
#     annotate(
#       "segment",
#       x = hseg_args$x,
#       xend = hseg_args$xend,
#       y = hseg_args$y,
#       yend = hseg_args$yend,
#       arrow = arrow(
#         length = unit(hseg_args$arrow_length, "cm"),
#         ends = hseg_args$ends,
#         type = hseg_args$type
#       ),
#       linewidth = hseg_args$linewidth,
#       colour = hseg_args$colour
#     ) +
#     
#     annotate(
#       "text",
#       x = htext_args$x,
#       y = htext_args$y,
#       label = htext_args$label,
#       size = htext_args$size
#     ) +
#     
#     annotate(
#       "segment",
#       x = vseg_args$x,
#       xend = vseg_args$xend,
#       y = vseg_args$y,
#       yend = vseg_args$yend,
#       arrow = arrow(
#         length = unit(vseg_args$arrow_length, "cm"),
#         ends = vseg_args$ends,
#         type = vseg_args$type
#       ),
#       linewidth = vseg_args$linewidth,
#       colour = vseg_args$colour
#     ) +
#     
#     annotate(
#       "text",
#       x = vtext_args$x,
#       y = vtext_args$y,
#       label = vtext_args$label,
#       angle = vtext_args$angle,
#       size = vtext_args$size
#     ) +
#     
#     ggplot2::theme_bw(base_size = theme_args$base_size) +
#     
#     facet_wrap(
#       ~sample_size,
#       labeller = labeller(sample_size = facet_labs)
#     ) +
#     
#     theme(
#       strip.text = element_text(size = theme_args$strip_text_size, color = "black"),
#       legend.title = element_text(size = theme_args$legend_title_size, color = "black"),
#       legend.text = element_text(size = theme_args$legend_text_size, color = "black"),
#       axis.text.x = element_text(size = theme_args$axis_text_x_size, color = "black"),
#       axis.text.y = element_text(size = theme_args$axis_text_y_size, color = "black")
#     )
#   
#   gg_2dimc
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
sample_size_grid <- function(target_powers, vary_vals, mod,
                             vary = c("abs_lfc", "logmean"),
                             fixed_logmean = NULL,
                             fixed_abs_lfc = NULL,
                             xmin = log2(5),
                             xmax = log2(500)) {
  
  vary <- match.arg(vary)
  
  do.call(rbind, lapply(vary_vals, function(val) {
    
    ss_vals <- sapply(target_powers, function(tp) {
      
      logmean_val <- if (vary == "logmean") val else fixed_logmean
      abs_lfc_val <- if (vary == "abs_lfc") val else fixed_abs_lfc
      
      ss_solver(
        target_power = tp,
        logmean = logmean_val,
        abs_lfc = abs_lfc_val,
        model = mod,
        xmin = xmin,
        xmax = xmax
      )
    })
    
    out <- data.frame(
      sample_size = ss_vals,
      power = target_powers
    )
    
    out[[vary]] <- val
    out
  }))
}


sample_size_est <- function(target_powers,lfc_vals, mod, logmean = 5, 
                   ssmin = 5, ssmax = 500){
  
  ss <- do.call(rbind, lapply(lfc_vals, function(lfc) {
    
    ss_vals <- sapply(target_powers, function(tp) {
      ss_solver(
        target_power = tp,
        logmean =  logmean,
        abs_lfc = lfc,
        model = mod,
        xmin = log2(ssmin),
        xmax = log2(ssmax)
      )
    })
    
    data.frame(
      sample_size = ss_vals,
      power = target_powers,
      abs_lfc = lfc
    )
  }))
  
  ss
  
}




