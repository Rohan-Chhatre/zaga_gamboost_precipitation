library(data.table)
library(ggplot2)

files <- list(
  Persistence = "/Users/rohan/Documents/Rohan/Precip_project/github/output/persistent_forecasts/rolling_validation/meriden_persistent_only/rolling_predictions.csv",
  ARIMA = "/Users/rohan/Documents/Rohan/Precip_project/github/output/arima/rolling_validation/meriden_arima/rolling_predictions.csv",
  RF = "/Users/rohan/Documents/Rohan/Precip_project/github/output/random_forest_new/rolling_recursive_validation/meriden_model1/rolling_predictions.csv",
  LSTM = '/Users/rohan/Documents/Rohan/Precip_project/github/output/lstm/rolling_recursive_validation_2/meriden_model1/rolling_predictions.csv',
  Model1 = "/Users/rohan/Documents/Rohan/Precip_project/github/output/gamboostlss/rolling_recursive_validation/meriden_model1/rolling_predictions.csv",
  Model2 = "/Users/rohan/Documents/Rohan/Precip_project/github/output/gamboostlss/rolling_recursive_validation/meriden_model2/rolling_predictions.csv",
  Model3 = "/Users/rohan/Documents/Rohan/Precip_project/github/output/gamboostlss/rolling_recursive_validation/meriden_allstations_model/rolling_predictions.csv"
)

model_cols <- c(
  "Persistence" = "#E69F00",
  "ARIMA"       = "#56B4E9",
  "RF"          = "#009E73",
  "LSTM"        = "#6A3D9A",
  "Model 1"     = "#0072B2",
  "Model 2"     = "#D55E00",
  "Model 3"     = "#CC79A7"
)
pred_col <- function(dt) {
  if ("zaga" %in% names(dt)) return("zaga")
  if ("arimax" %in% names(dt)) return("arimax")
  if ("arima" %in% names(dt)) return("arima")
  if ("rf" %in% names(dt)) return("rf")
  if ("prediction" %in% names(dt)) return("prediction")
  if ("persistence" %in% names(dt)) return("persistence")
  #if ("naive" %in% names(dt)) return("naive")
  if ("lstm" %in% names(dt)) return("lstm")
  stop("No prediction column found.")
}

all_lead_metrics <- rbindlist(
  lapply(names(files), function(m) {
    
    dt <- fread(files[[m]])
    pcol <- pred_col(dt)
    
    setorder(dt, origin, DATE)
    
    dt[, lead_hour := seq_len(.N), by = origin]
    
    dt[
      ,
      .(
        rmse = sqrt(mean((actual - get(pcol))^2, na.rm = TRUE)),
        mae = mean(abs(actual - get(pcol)), na.rm = TRUE),
        mean_error = mean(actual - get(pcol), na.rm = TRUE),
        error_sd = sd(actual - get(pcol), na.rm = TRUE),
        n = .N
      ),
      by = lead_hour
    ][, model := m]
  }),
  fill = TRUE
)

all_lead_metrics[, model := factor(
  model,
  levels = c("Persistence", "ARIMA", "RF", "LSTM", "Model1", "Model2", "Model3"),
  labels = c("Persistence", "ARIMA", "RF", "LSTM", "Model 1", "Model 2", "Model 3")
)]

################################################################################

library(data.table)
library(ggplot2)

# ============================================================
# Prepare origin-lead squared errors
# ============================================================

all_origin_lead_errors <- rbindlist(
  lapply(names(files), function(m) {
    
    dt <- fread(files[[m]])
    pcol <- pred_col(dt)
    
    setorder(dt, origin, DATE)
    dt[, lead_hour := seq_len(.N), by = origin]
    
    dt[
      ,
      .(
        sq_error = (actual - get(pcol))^2
      ),
      by = .(origin, lead_hour)
    ][, model := m]
  }),
  fill = TRUE
)

all_origin_lead_errors[, model := factor(
  model,
  levels = c("Persistence", "ARIMA", "RF", "LSTM", "Model1", "Model2", "Model3"),
  labels = c("Persistence", "ARIMA", "RF", "LSTM", "Model 1", "Model 2", "Model 3")
)]

# ============================================================
# Mean RMSE + 5th–95th percentile ribbon
# ============================================================
model_cols <- c(
  "Persistence" = "#E69F00",
  "ARIMA"       = "#56B4E9",
  "RF"          = "#009E73",
  "LSTM"        = "#6A3D9A",
  "Model 1"     = "#0072B2",
  "Model 2"     = "#D55E00",
  "Model 3"     = "#CC79A7"
)

# RMSE by forecast lead time
p_rmse <- ggplot(all_lead_metrics, aes(x = lead_hour, y = rmse, color = model)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = model_cols) +
  labs(
    x = "Forecast lead time (hours)",
    y = "RMSE (mm)",
    color = "Model",
    title = "Overall RMSE vs lead time"
  ) +
  theme_bw(base_size = 14) +
  
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p_rmse)

ggsave(
  filename = "/Users/rohan/Documents/Rohan/Precip_project/github/output/figure5_rmse_by_lead_time.png",
  plot = p_rmse,
  width = 8,
  height = 5,
  dpi = 300
)


#####

library(data.table)
library(ggplot2)
library(lubridate)

# ============================================================
# Helpers
# ============================================================

get_season <- function(d) {
  m <- lubridate::month(d)
  fifelse(m %in% c(6, 7, 8, 9, 10, 11), "Hurricane", "Calm")
}

bucket_fun <- function(x) {
  fifelse(
    x == 0, "Dry",
    fifelse(
      x <= 2, "Light",
      fifelse(x <= 10, "Moderate", "Heavy")
    )
  )
}

model_cols <- c(
  "Persistence" = "#E69F00",
  "ARIMA"       = "#56B4E9",
  "RF"          = "#009E73",
  "LSTM"        = "#6A3D9A",
  "Model 1"     = "#0072B2",
  "Model 2"     = "#D55E00",
  "Model 3"     = "#CC79A7"
)

# ============================================================
# Build combined prediction-error data
# ============================================================

all_pred_errors <- rbindlist(
  lapply(names(files), function(m) {
    
    dt <- fread(files[[m]])
    pcol <- pred_col(dt)
    
    dt[, DATE := as.POSIXct(DATE, tz = "UTC")]
    dt[, origin := as.POSIXct(origin, tz = "UTC")]
    
    setorder(dt, origin, DATE)
    dt[, lead_hour := seq_len(.N), by = origin]
    
    dt[, error := actual - get(pcol)]
    dt[, season := get_season(DATE)]
    dt[, rainfall_bucket := bucket_fun(actual)]
    
    dt[, .(
      model = m,
      origin,
      DATE,
      lead_hour,
      actual,
      error,
      season,
      rainfall_bucket
    )]
  }),
  fill = TRUE
)

all_pred_errors[, model := factor(
  model,
  levels = c("Persistence", "ARIMA", "RF", "LSTM", "Model1", "Model2", "Model3"),
  labels = c("Persistence", "ARIMA", "RF", "LSTM", "Model 1", "Model 2", "Model 3")
)]

all_pred_errors[, season := factor(
  season,
  levels = c("Calm", "Hurricane")
)]

all_pred_errors[, rainfall_bucket := factor(
  rainfall_bucket,
  levels = c("Dry", "Light", "Moderate", "Heavy")
)]

# ============================================================
# 2) RMSE by lead time and season
# ============================================================

season_lead_rmse <- all_pred_errors[
  ,
  .(
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    n = .N
  ),
  by = .(season, model, lead_hour)
]

p_rmse_season <- ggplot(
  season_lead_rmse,
  aes(x = lead_hour, y = rmse, color = model)
) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = model_cols) +
  facet_wrap(~ season, ncol = 2, scales = "free_y") +
  labs(
    title = "RMSE Across Forecast Lead Times by Season",
    x = "Forecast lead time (hours)",
    y = "RMSE (mm)",
    color = "Model"
  ) +
  theme_bw(base_size = 14) +
  
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )


print(p_rmse_season)

ggsave(
  filename = "/Users/rohan/Documents/Rohan/Precip_project/github/output/rmse_by_lead_time_season.png",
  plot = p_rmse_season,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# 3) RMSE by lead time and actual rainfall bucket
# ============================================================

bucket_lead_rmse <- all_pred_errors[
  ,
  .(
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    n = .N
  ),
  by = .(rainfall_bucket, model, lead_hour)
]

p_rmse_bucket <- ggplot(
  bucket_lead_rmse,
  aes(x = lead_hour, y = rmse, color = model)
) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = model_cols) +
  facet_wrap(~ rainfall_bucket, ncol = 2, scales = "free_y", nrow=2) +
  labs(
    title = "RMSE Across Forecast Lead Times by rainfall bucket summary",
    x = "Forecast lead time (hours)",
    y = "RMSE (mm)",
    color = "Model"
  ) +
  theme_bw(base_size = 14) +
  
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p_rmse_bucket)

ggsave(
  filename = "/Users/rohan/Documents/Rohan/Precip_project/github/output/rmse_by_lead_time_rainfall_bucket.png",
  plot = p_rmse_bucket,
  width = 9,
  height = 10,
  dpi = 300
)

