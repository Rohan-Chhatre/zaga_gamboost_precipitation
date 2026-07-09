library(data.table)
library(ggplot2)
library(lubridate)

files <- list(
  Persistence = "/Users/rohan/Documents/Rohan/Precip_project/github/output/persistent_forecasts/rolling_validation/meriden_persistent_only/rolling_predictions.csv",
  ARIMA = "/Users/rohan/Documents/Rohan/Precip_project/github/output/arima/rolling_validation/meriden_arima/rolling_predictions.csv",
  RF = "/Users/rohan/Documents/Rohan/Precip_project/github/output/random_forest_new/rolling_recursive_validation/meriden_model1/rolling_predictions.csv",
  LSTM = '/Users/rohan/Documents/Rohan/Precip_project/github/output/lstm/rolling_recursive_validation_2/meriden_model1/rolling_predictions.csv',
  Model1 = "/Users/rohan/Documents/Rohan/Precip_project/github/output/gamboostlss/rolling_recursive_validation/meriden_model1/rolling_predictions.csv",
  Model2 = "/Users/rohan/Documents/Rohan/Precip_project/github/output/gamboostlss/rolling_recursive_validation/meriden_model2/rolling_predictions.csv",
  Model3 = "/Users/rohan/Documents/Rohan/Precip_project/github/output/gamboostlss/rolling_recursive_validation/meriden_allstations_model/rolling_predictions.csv"
)

outdir <- "/Users/rohan/Documents/Rohan/Precip_project/github/output/event_comparison_5events"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pred_col <- function(dt) {
  if ("zaga" %in% names(dt)) return("zaga")
  if ("arima" %in% names(dt)) return("arima")
  if ("arimax" %in% names(dt)) return("arimax")
  if ("rf" %in% names(dt)) return("rf")
  if ("persistence" %in% names(dt)) return("persistence")
  #if ("naive" %in% names(dt)) return("naive")
  if ("lstm" %in% names(dt)) return("lstm")
  stop("No prediction column found.")
}

rmse <- function(actual, pred) {
  sqrt(mean((actual - pred)^2, na.rm = TRUE))
}

mae <- function(actual, pred) {
  mean(abs(actual - pred), na.rm = TRUE)
}

# ------------------------------------------------------------
# Event windows
# You can adjust these if your station-level rainfall peak differs.
# preferred_origin is chosen BEFORE the event starts.
# ------------------------------------------------------------

events <- data.table(
  event = c(
    "December 2020 Nor'easter",
    "Tropical Storm Elsa",
    "Hurricane Henri",
    "Hurricane Ida Remnants",
    "Summer Convective Storm"
  ),
  event_start = as.POSIXct(c(
    "2020-12-05 00:00:00",
    "2021-07-08 00:00:00",
    "2021-08-21 00:00:00",
    "2021-09-01 00:00:00",
    "2021-07-17 00:00:00"
  ), tz = "UTC"),
  event_end = as.POSIXct(c(
    "2020-12-07 00:00:00",
    "2021-07-10 00:00:00",
    "2021-08-24 00:00:00",
    "2021-09-03 00:00:00",
    "2021-07-19 00:00:00"
  ), tz = "UTC"),
  preferred_origin = as.POSIXct(c(
    "2020-12-04 12:00:00",
    "2021-07-07 12:00:00",
    "2021-08-20 12:00:00",
    "2021-08-31 12:00:00",
    "2021-07-16 12:00:00"
  ), tz = "UTC")
)

all_event_metrics <- list()
all_event_predictions <- list()

for (e in seq_len(nrow(events))) {
  
  ev_name <- events$event[e]
  ev_start <- events$event_start[e]
  ev_end <- events$event_end[e]
  target_origin <- events$preferred_origin[e]
  
  plot_list <- list()
  metric_list <- list()
  
  for (mod in names(files)) {
    
    dt <- fread(files[[mod]])
    
    dt[, DATE := as.POSIXct(DATE, tz = "UTC")]
    dt[, origin := as.POSIXct(origin, tz = "UTC")]
    
    pcol <- pred_col(dt)
    
    # Choose a forecast origin that starts before the event
    # and whose 120-hour horizon covers the whole event.
    valid_origins <- unique(
      dt[
        origin <= ev_start &
          origin + hours(120) >= ev_end,
        origin
      ]
    )
    
    if (length(valid_origins) == 0) {
      warning(paste("No valid origin found for", ev_name, "and model", mod))
      next
    }
    
    chosen_origin <- valid_origins[
      which.min(abs(as.numeric(difftime(valid_origins, target_origin, units = "hours"))))
    ]
    
    event_dt <- dt[
      origin == chosen_origin &
        DATE >= ev_start &
        DATE <= ev_end
    ]
    
    if (nrow(event_dt) == 0) {
      warning(paste("No event rows found for", ev_name, "and model", mod))
      next
    }
    
    event_dt[, model := mod]
    event_dt[, pred := get(pcol)]
    event_dt[, model := fcase(
      mod == "Persistence", "Persistence",
      mod == "ARIMA",       "ARIMA",
      mod == "RF",          "RF",
      mod == "LSTM",        "LSTM",
      mod == "Model1",      "Model 1",
      mod == "Model2",      "Model 2",
      mod == "Model3",      "Model 3"
    )]
    metric_list[[mod]] <- data.table(
      event = ev_name,
      model = mod,
      origin = chosen_origin,
      n = nrow(event_dt),
      rmse = rmse(event_dt$actual, event_dt$pred),
      mae = mae(event_dt$actual, event_dt$pred),
      max_observed = max(event_dt$actual, na.rm = TRUE),
      max_predicted = max(event_dt$pred, na.rm = TRUE)
    )
    
    plot_list[[mod]] <- event_dt[, .(
      event = ev_name,
      origin = chosen_origin,
      DATE,
      actual,
      pred,
      lo95 = if ("lo95" %in% names(event_dt)) lo95 else NA_real_,
      hi95 = if ("hi95" %in% names(event_dt)) hi95 else NA_real_,
      model
    )]
  }
  
  plot_dt <- rbindlist(plot_list, fill = TRUE)
  metrics_dt <- rbindlist(metric_list, fill = TRUE)
  
  all_event_metrics[[ev_name]] <- metrics_dt
  all_event_predictions[[ev_name]] <- plot_dt
  
  observed_dt <- unique(plot_dt[, .(DATE, actual)])
  
  observed_long <- observed_dt[, .(
    DATE,
    series = "Observed",
    value = actual
  )]
  
  pred_long <- plot_dt[, .(
    DATE,
    series = model,
    value = pred
  )]
  
  long_dt <- rbind(observed_long, pred_long, fill = TRUE)
  
  long_dt[, series := factor(
    series,
    levels = c("Observed", "Persistence", "ARIMA", "RF", "LSTM", "Model 1", "Model 2", "Model 3")
  )]
  
  model_cols <- c(
    "Persistence" = "#E69F00",
    "ARIMA"       = "#56B4E9",
    "RF"          = "#009E73",
    "LSTM"        = "#6A3D9A",
    "Model 1"     = "#0072B2",
    "Model 2"     = "#D55E00",
    "Model 3"     = "#CC79A7"
  )
  
  ribbon_dt <- plot_dt[
    !is.na(lo95) & !is.na(hi95)
  ]
  
  p1 <- ggplot() +
    
    ## 95% prediction interval bands
    geom_ribbon(
      data = ribbon_dt,
      aes(
        x = DATE,
        ymin = lo95,
        ymax = hi95,
        fill = model
      ),
      alpha = 0.16,
      colour = NA
    ) +
    
    ## Observed
    geom_line(
      data = subset(long_dt, series == "Observed"),
      aes(DATE, value),
      colour = "black",
      linewidth = 1.4
    ) +
    
    ## Forecasts
    geom_line(
      data = subset(long_dt, series != "Observed"),
      aes(DATE, value, colour = series),
      linewidth = 0.9
    ) +
    
    scale_colour_manual(values = model_cols) +
    scale_fill_manual(values = model_cols) +
    
    labs(
      title = paste0(ev_name, ": observed and forecast precipitation"),
      subtitle = paste0("Event window: ", ev_start, " to ", ev_end, " UTC"),
      x = "Date",
      y = "Hourly precipitation (mm)",
      colour = "Forecast model",
      fill = "95% interval"
    ) +
    
    theme_bw(base_size = 14) +
    
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
  print(p1)
  
  ggsave(
    filename = file.path(
      outdir,
      paste0(gsub("[^A-Za-z0-9]", "_", ev_name), "_forecast_plot.png")
    ),
    plot = p1,
    width = 10,
    height = 5.5,
    dpi = 300
  )
  
  p2 <- ggplot(metrics_dt, aes(x = reorder(model, rmse), y = rmse, fill = model)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    labs(
      title = paste0(ev_name, ": RMSE by model"),
      x = "Model",
      y = "RMSE"
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  print(p2)
  
  ggsave(
    filename = file.path(
      outdir,
      paste0(gsub("[^A-Za-z0-9]", "_", ev_name), "_rmse_barplot.png")
    ),
    plot = p2,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  fwrite(
    metrics_dt,
    file.path(
      outdir,
      paste0(gsub("[^A-Za-z0-9]", "_", ev_name), "_rmse_metrics.csv")
    )
  )
}

event_metrics <- rbindlist(all_event_metrics, fill = TRUE)
event_predictions <- rbindlist(all_event_predictions, fill = TRUE)

fwrite(
  event_metrics,
  file.path(outdir, "event_level_rmse_metrics.csv")
)

fwrite(
  event_predictions,
  file.path(outdir, "event_level_predictions.csv")
)

event_metrics

