# ============================================================
# ROLLING-WINDOW VALIDATION
# ZAGA gamboostLSS
# Model 1: m precipitation ~ m precipitation lags only
#
# Validation design:
# Training window  = previous 5 years
# Forecast horizon = 120 hours
# Step size        = 1 hours
# Validation start = 2021-01-01 00:00:00 UTC
# Validation end   = 2024-01-01 00:00:00 UTC
#
# Expected model fits:
# approx 8000 fits
# ============================================================

# -----------------------------
# Personal library + packages
# -----------------------------
#userlib <- "/Users/rohan/Documents/Rohan/"
#dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
#.libPaths(c(userlib, .libPaths()))

req <- c(
  "data.table",
  "lubridate",
  "gamboostLSS",
  "mboost",
  "gamlss.dist",
  "foreach",
  "doParallel"
)

#for (p in req) {
#  if (!requireNamespace(p, quietly = TRUE)) {
#    install.packages(p, lib = userlib, repos = "https://cloud.r-project.org")
#  }
#}

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(gamboostLSS)
  library(mboost)
  library(gamlss.dist)
  library(foreach)
  library(doParallel)
})

# -----------------------------
# Helper functions
# -----------------------------
numify <- function(x) {
  if (is.character(x)) x[x == ""] <- NA_character_
  as.numeric(x)
}

rmse <- function(actual, pred) {
  sqrt(mean((actual - pred)^2, na.rm = TRUE))
}

mae <- function(actual, pred) {
  mean(abs(actual - pred), na.rm = TRUE)
}

bias <- function(actual, pred) {
  mean(pred - actual, na.rm = TRUE)
}

smape <- function(actual, pred) {
  
  denom <- (abs(actual) + abs(pred))/2
  
  idx <- denom > 0
  
  if(sum(idx)==0) return(NA_real_)
  
  100 * mean(
    abs(actual[idx]-pred[idx]) / denom[idx],
    na.rm=TRUE
  )
  
}

get_season <- function(d) {
  m <- lubridate::month(d)
  
  data.table::fifelse(
    m %in% c(6, 7, 8, 9, 10, 11),
    "Hurricane",
    "Calm"
  )
}

rmse_peaks <- function(actual, pred, q = 0.95) {
  thr <- as.numeric(quantile(actual, probs = q, na.rm = TRUE))
  idx <- which(actual >= thr)
  
  if (length(idx) < 5) {
    return(list(threshold = thr, n = length(idx), rmse = NA_real_))
  }
  
  list(
    threshold = thr,
    n = length(idx),
    rmse = rmse(actual[idx], pred[idx])
  )
}

qmape <- function(actual, pred) {
  denom <- sum(actual, na.rm = TRUE)
  if (is.na(denom) || denom == 0) return(NA_real_)
  sum(abs(actual - pred), na.rm = TRUE) / denom
}

mase <- function(actual, pred, naive_actual, naive_pred) {
  model_mae <- mae(actual, pred)
  naive_mae <- mae(naive_actual, naive_pred)
  
  if (is.na(naive_mae) || naive_mae == 0) return(NA_real_)
  model_mae / naive_mae
}

bucket_fun <- function(x){
  
  fifelse(
    x == 0,
    "Dry",
    fifelse(
      x <= 2,
      "Light",
      fifelse(
        x <= 10,
        "Moderate",
        "Heavy"
      )
    )
  )
  
}
# -----------------------------
# Paths
# -----------------------------
path_G <- "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/concatenated_files/USW00054788_merged_imputed_2004_2024.csv"

outdir <- "/Users/rohan/Documents/Rohan/Precip_project/github/output/gamboostlss/rolling_validation/meriden_model1"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
persist_outdir <-
  "/Users/rohan/Documents/Rohan/Precip_project/github/output/persistent_forecasts/rolling_validation/meriden_persistent"

dir.create(
  persist_outdir,
  recursive = TRUE,
  showWarnings = FALSE
)
# -----------------------------
# Read Groton data
# -----------------------------
cols_needed <- c(
  "STATION",
  "DATE",
  "HourlyPrecipitation"
)

DT <- fread(path_G, select = cols_needed)

DT[, DATE := as.POSIXct(DATE, tz = "UTC")]
setorder(DT, DATE)

DT <- DT[
  DATE >= as.POSIXct("2004-01-01 00:00:00", tz = "UTC") &
    DATE <= as.POSIXct("2024-10-09 00:00:00", tz = "UTC")
]

# -----------------------------
# Build  features
# -----------------------------
L <- 6

DT[, rain_M := numify(HourlyPrecipitation)]
DT[, doy := yday(DATE)]

for (i in 1:L) {
  DT[, paste0("rain_M_lag", i) := shift(rain_M, n = i, type = "lag")]
  DT[, paste0("rain_M_lag", i) := log1p(get(paste0("rain_M_lag", i)))]
}

vars <- c(
  "DATE",
  "doy",
  "rain_M",
  paste0("rain_M_lag", 1:L)
)

DT <- DT[, ..vars]
DT <- na.omit(DT)

# -----------------------------
# Rolling validation settings
# -----------------------------
train_years <- 4
horizon_hours <- 120
step_hours <- 1 #120

validation_start <- as.POSIXct("2020-12-01 00:00:00", tz = "UTC")
validation_end   <- as.POSIXct("2021-12-01 00:00:00", tz = "UTC") #as.POSIXct("2024-01-01 00:00:00", tz = "UTC")

forecast_origins <- seq(
  from = validation_start,
  to = validation_end - hours(horizon_hours),
  by = paste(step_hours, "hours")
)

cat("Number of forecast origins:", length(forecast_origins), "\n")
cat("Number of  model fits:", length(forecast_origins), "\n")

# -----------------------------
#  formula
# -----------------------------
df_bbs <- 3

bbs_term <- function(v) {
  paste0("bbs(", v, ", df = ", df_bbs, ")")
}

bols_term <- function(v) {
  paste0("bols(", v, ", intercept = TRUE)")
}

rhs_mu <- paste(
  c(
    bbs_term("doy"),
    bols_term("rain_M_lag1"),
    bbs_term("rain_M_lag2"),
    bbs_term("rain_M_lag3"),
    bbs_term("rain_M_lag4"),
    bbs_term("rain_M_lag5"),
    bbs_term("rain_M_lag6")
  ),
  collapse = " + "
)

rhs_sigma <- paste(
  c(
    bbs_term("doy"),
    bols_term("rain_M_lag1"),
    bbs_term("rain_M_lag2")
  ),
  collapse = " + "
)

rhs_nu <- paste(
  c(
    bbs_term("doy"),
    bbs_term("rain_M_lag1")
  ),
  collapse = " + "
)

form_list <- list(
  mu    = as.formula(paste("rain_M ~", rhs_mu)),
  sigma = as.formula(paste("rain_M ~", rhs_sigma)),
  nu    = as.formula(paste("rain_M ~", rhs_nu))
)

# -----------------------------
# Storage objects
# -----------------------------
all_metrics <- list()
all_predictions <- list()
all_runtime <- list()
all_errors <- list()
#forecast_origins = forecast_origins[1]

# -----------------------------
# Rolling validation loop
# -----------------------------
library(doSNOW)

t0_total <- Sys.time()

ncores <- 10

cl <- parallel::makeCluster(ncores)
doSNOW::registerDoSNOW(cl)

cat("Using", ncores, "cores\n")

pb <- txtProgressBar(max = length(forecast_origins), style = 3)

progress <- function(n) {
  setTxtProgressBar(pb, n)
}

opts <- list(progress = progress)

results <- foreach(
  k = seq_along(forecast_origins),
  .options.snow = opts,
  .packages = c(
    "data.table",
    "lubridate",
    "gamboostLSS",
    "mboost",
    "gamlss.dist"
  ),
  .export = c(
    "DT",
    "forecast_origins",
    "train_years",
    "horizon_hours",
    "form_list",
    "rmse",
    "mae",
    "bias",
    "smape",
    "qmape",
    "mase",
    "rmse_peaks",
    "get_season",
    "bucket_fun"
  ),
  .errorhandling = "pass"
) %dopar% {
  
  origin <- forecast_origins[k]
  
  train_start <- origin %m-% years(train_years)
  train_end   <- origin
  test_start  <- origin + hours(1)
  test_end    <- origin + hours(horizon_hours)
  
  t0_fit <- Sys.time()
  
  tryCatch({
    
    train_dt <- DT[
      DATE >= train_start &
        DATE <= train_end
    ]
    
    test_dt <- DT[
      DATE >= test_start &
        DATE <= test_end
    ]
    
    train_dt <- na.omit(train_dt)
    test_dt  <- na.omit(test_dt)
    
    if (nrow(test_dt) < horizon_hours) {
      stop("Test set has fewer than 120 observations.")
    }
    
    dtr <- as.data.frame(train_dt[, !"DATE"])
    dte <- as.data.frame(test_dt[, !"DATE"])
    
    Y_tr <- train_dt$rain_M
    Y_te <- test_dt$rain_M
    
    # Naive persistence forecast: y_hat_t = y_{t-1}
    naive_te <- expm1(test_dt$rain_M_lag1)
    
    persist_k <- data.table(
      origin_id = k,
      origin = origin,
      DATE = test_dt$DATE,
      actual = Y_te,
      persistence = naive_te
    )
    
    # Training denominator for MASE
    naive_tr_actual <- train_dt$rain_M
    naive_tr_pred   <- expm1(train_dt$rain_M_lag1)
    
    gc()
    
    fit.gb <- gamboostLSS(
      formula  = form_list,
      data     = dtr,
      families = as.families("ZAGA"),
      control  = boost_control(mstop = 300, nu = 0.01)
    )
    
    mu_te <- as.numeric(
      predict(fit.gb, newdata = dte, parameter = "mu", type = "response")
    )
    
    sigma_te <- as.numeric(
      predict(fit.gb, newdata = dte, parameter = "sigma", type = "response")
    )
    
    nu_te <- as.numeric(
      predict(fit.gb, newdata = dte, parameter = "nu", type = "response")
    )
    
    pred_te <- (1 - nu_te) * mu_te
    
    lo80_te <- gamlss.dist::qZAGA(
      0.10, mu = mu_te, sigma = sigma_te, nu = nu_te
    )
    
    hi80_te <- gamlss.dist::qZAGA(
      0.90, mu = mu_te, sigma = sigma_te, nu = nu_te
    )
    
    lo95_te <- gamlss.dist::qZAGA(
      0.025, mu = mu_te, sigma = sigma_te, nu = nu_te
    )
    
    hi95_te <- gamlss.dist::qZAGA(
      0.975, mu = mu_te, sigma = sigma_te, nu = nu_te
    )
    
    peak90 <- rmse_peaks(Y_te, pred_te, q = 0.90)
    peak95 <- rmse_peaks(Y_te, pred_te, q = 0.95)
    peak99 <- rmse_peaks(Y_te, pred_te, q = 0.99)
    
    metrics_k <- data.table(
      origin_id = k,
      origin = origin,
      train_start = train_start,
      train_end = train_end,
      test_start = test_start,
      test_end = test_end,
      season = get_season(test_start),
      n_train = nrow(train_dt),
      n_test = nrow(test_dt),
      
      rmse = rmse(Y_te, pred_te),
      mae = mae(Y_te, pred_te),
      bias = bias(Y_te, pred_te),
      mase = mase(
        actual = Y_te,
        pred = pred_te,
        naive_actual = naive_tr_actual,
        naive_pred = naive_tr_pred
      ),
      smape = smape(Y_te, pred_te),
      naive_rmse = rmse(Y_te, naive_te),
      naive_mae = mae(Y_te, naive_te),
      
      peak90_threshold = peak90$threshold,
      peak90_n = peak90$n,
      peak90_rmse = peak90$rmse,
      
      peak95_threshold = peak95$threshold,
      peak95_n = peak95$n,
      peak95_rmse = peak95$rmse,
      
      peak99_threshold = peak99$threshold,
      peak99_n = peak99$n,
      peak99_rmse = peak99$rmse
    )
    
    pred_k <- data.table(
      origin_id = k,
      origin = origin,
      DATE = test_dt$DATE,
      season = get_season(test_dt$DATE),
      actual = Y_te,
      mu = mu_te,
      sigma = sigma_te,
      nu = nu_te,
      zaga = pred_te,
      naive = naive_te,   # ADD THIS
      lo80 = lo80_te,
      hi80 = hi80_te,
      lo95 = lo95_te,
      hi95 = hi95_te
    )
    
    # Bucket based on actual observed precipitation
    pred_k[, rainfall_bucket := bucket_fun(actual)]
    
    t1_fit <- Sys.time()
    
    runtime_k <- data.table(
      origin_id = k,
      origin = origin,
      status = "success",
      runtime_seconds = as.numeric(difftime(t1_fit, t0_fit, units = "secs"))
    )
    
    list(
      metrics = metrics_k,
      predictions = pred_k,
      runtime = runtime_k,
      error = NULL,
      persistence = persist_k
    )
    
  }, error = function(e) {
    
    t1_fit <- Sys.time()
    
    runtime_k <- data.table(
      origin_id = k,
      origin = origin,
      status = "failed",
      runtime_seconds = as.numeric(difftime(t1_fit, t0_fit, units = "secs"))
    )
    
    error_k <- data.table(
      origin_id = k,
      origin = origin,
      train_start = train_start,
      train_end = train_end,
      test_start = test_start,
      test_end = test_end,
      error_message = conditionMessage(e)
    )
    
    list(
      metrics = NULL,
      predictions = NULL,
      runtime = runtime_k,
      error = error_k,
      persistence = NULL
    )
  })
}
close(pb)
parallel::stopCluster(cl)
# -----------------------------
# Collect parallel outputs
# -----------------------------
all_metrics     <- lapply(results, `[[`, "metrics")
all_predictions <- lapply(results, `[[`, "predictions")
all_runtime     <- lapply(results, `[[`, "runtime")
all_errors      <- lapply(results, `[[`, "error")
all_persistence <- lapply(results, `[[`, "persistence")

# -----------------------------
# Final outputs
# -----------------------------
metrics <- rbindlist(all_metrics, fill = TRUE)
predictions <- rbindlist(all_predictions, fill = TRUE)
runtime <- rbindlist(all_runtime, fill = TRUE)
persistence <- rbindlist(all_persistence, fill = TRUE)

fwrite(
  persistence,
  file.path(persist_outdir, "persistent_predictions.csv")
)
fwrite(metrics, file.path(outdir, "rolling_metrics.csv"))
fwrite(predictions, file.path(outdir, "rolling_predictions.csv"))
fwrite(runtime, file.path(outdir, "runtime.csv"))

if (length(all_errors) > 0) {
  errors <- rbindlist(all_errors, fill = TRUE)
  fwrite(errors, file.path(outdir, "errors.csv"))
}

# -----------------------------
# Seasonal summary
# -----------------------------
seasonal_summary <- metrics[
  ,
  .(
    n_windows = .N,
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    median_rmse = median(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    mean_peak95_rmse = mean(peak95_rmse, na.rm = TRUE),
    #mean_qmape = mean(qmape, na.rm = TRUE),
    mean_smape = mean(smape, na.rm = TRUE),
    mean_mase = mean(mase, na.rm = TRUE)
  ),
  by = season
]

fwrite(
  seasonal_summary,
  file.path(outdir, "seasonal_summary.csv")
)

# -----------------------------
# Overall summary
# -----------------------------
overall_summary <- metrics[
  ,
  .(
    n_windows = .N,
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    median_rmse = median(rmse, na.rm = TRUE),
    min_rmse = min(rmse, na.rm = TRUE),
    max_rmse = max(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    mean_peak90_rmse = mean(peak90_rmse, na.rm = TRUE),
    mean_peak95_rmse = mean(peak95_rmse, na.rm = TRUE),
    mean_peak99_rmse = mean(peak99_rmse, na.rm = TRUE),
    #mean_qmape = mean(qmape, na.rm = TRUE),
    #sd_qmape = sd(qmape, na.rm = TRUE),
    #median_qmape = median(qmape, na.rm = TRUE),
    mean_smape = mean(smape, na.rm = TRUE),
    sd_smape = sd(smape, na.rm = TRUE),
    median_smape = median(smape, na.rm = TRUE),
    mean_mase = mean(mase, na.rm = TRUE),
    sd_mase = sd(mase, na.rm = TRUE),
    median_mase = median(mase, na.rm = TRUE)
  )
]

fwrite(
  overall_summary,
  file.path(outdir, "overall_summary.csv")
)

bucket_summary <- predictions[
  ,
  .(
    n = .N,
    rmse = rmse(actual, zaga),
    mae = mae(actual, zaga),
    bias = bias(actual, zaga),
    smape = smape(actual, zaga),
    qmape = qmape(actual, zaga),
    
    naive_mae = mae(actual, naive),
    
    mase = mae(actual, zaga) /
      mae(actual, naive)
  ),
  by = rainfall_bucket
]

fwrite(
  bucket_summary,
  file.path(outdir, "bucket_summary.csv")
)

# -----------------------------
# Total runtime
# -----------------------------
t1_total <- Sys.time()

total_runtime <- data.table(
  total_runtime_seconds = as.numeric(difftime(t1_total, t0_total, units = "secs")),
  total_runtime_hours = as.numeric(difftime(t1_total, t0_total, units = "hours")),
  n_forecast_origins = length(forecast_origins),
  n_success = nrow(runtime[status == "success"]),
  n_failed = nrow(runtime[status == "failed"])
)

fwrite(
  total_runtime,
  file.path(outdir, "total_runtime.csv")
)

cat("\nDONE.\n")
cat("Outputs saved in:\n")
cat(normalizePath(outdir), "\n")
cat("\nNumber of forecast origins:", length(forecast_origins), "\n")
cat("Successful fits:", nrow(runtime[status == "success"]), "\n")
cat("Failed fits:", nrow(runtime[status == "failed"]), "\n")
