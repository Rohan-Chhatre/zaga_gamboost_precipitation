# ============================================================
# ROLLING-WINDOW VALIDATION
# ZAGA gamboostLSS + NARR imputation
# MODEL: Meriden ~ Meriden + Bradley
# ============================================================

#userlib <- "/home/fbs24003/Rlibs/r452"
#dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
#.libPaths(c(userlib, .libPaths()))

req <- c(
  "data.table",
  "lubridate",
  "gamboostLSS",
  "mboost",
  "gamlss.dist",
  "zoo",
  "foreach",
  "doSNOW"
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
  library(zoo)
  library(foreach)
  library(doSNOW)
})

# -----------------------------
# Helpers
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
  if (sum(idx) == 0) return(NA_real_)
  100 * mean(abs(actual[idx] - pred[idx]) / denom[idx], na.rm = TRUE)
}

get_season <- function(d) {
  m <- lubridate::month(d)
  data.table::fifelse(
    m %in% c(6, 7, 8, 9, 10, 11),
    "Hurricane",
    "Calm"
  )
}

bucket_fun <- function(x) {
  data.table::fifelse(
    x == 0,
    "Dry",
    data.table::fifelse(
      x <= 2,
      "Light",
      data.table::fifelse(
        x <= 10,
        "Moderate",
        "Heavy"
      )
    )
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
# -----------------------------
# NARR helper
# -----------------------------
narr_to_hourly <- function(narr_path) {
  
  narr <- fread(narr_path)
  
  if ("time" %in% names(narr)) {
    narr[, time := as.POSIXct(time, tz = "UTC")]
    setnames(narr, "time", "TIME_3H")
  } else if ("DATE" %in% names(narr)) {
    narr[, DATE := as.POSIXct(DATE, tz = "UTC")]
    setnames(narr, "DATE", "TIME_3H")
  } else {
    stop("NARR file has no 'time' or 'DATE' column.")
  }
  
  if (all(c("uwnd", "vwnd") %in% names(narr))) {
    narr[, `:=`(u_3h = numify(uwnd), v_3h = numify(vwnd))]
  } else if (all(c("u_3h", "v_3h") %in% names(narr))) {
    narr[, `:=`(u_3h = numify(u_3h), v_3h = numify(v_3h))]
  } else {
    uv_cols <- tail(names(narr), 2)
    narr[, `:=`(
      u_3h = numify(get(uv_cols[1])),
      v_3h = numify(get(uv_cols[2]))
    )]
  }
  
  if ("dpt" %in% names(narr)) {
    narr[, dpt_3h := numify(dpt)]
  } else if ("dpt_narr" %in% names(narr)) {
    narr[, dpt_3h := numify(dpt_narr)]
  } else {
    narr[, dpt_3h := NA_real_]
  }
  
  if ("rhum" %in% names(narr)) {
    narr[, rh_3h := numify(rhum)]
  } else if ("RH" %in% names(narr)) {
    narr[, rh_3h := numify(RH)]
  } else if ("rh_3h" %in% names(narr)) {
    narr[, rh_3h := numify(rh_3h)]
  } else {
    narr[, rh_3h := NA_real_]
  }
  
  if ("prmsl" %in% names(narr)) {
    narr[, slp_3h := numify(prmsl)]
  } else if ("mslp" %in% names(narr)) {
    narr[, slp_3h := numify(mslp)]
  } else if ("slp_3h" %in% names(narr)) {
    narr[, slp_3h := numify(slp_3h)]
  } else {
    narr[, slp_3h := NA_real_]
  }
  
  narr <- narr[!is.na(TIME_3H)]
  setorder(narr, TIME_3H)
  
  if (!all(is.na(narr$dpt_3h))) {
    if (median(narr$dpt_3h, na.rm = TRUE) > 100) {
      narr[, dpt_3h := dpt_3h - 273.15]
    }
  }
  
  if (!all(is.na(narr$slp_3h))) {
    if (median(narr$slp_3h, na.rm = TRUE) > 2000) {
      narr[, slp_3h := slp_3h / 100]
    }
  }
  
  hour_grid <- data.table(
    TIME_3H = seq(min(narr$TIME_3H), max(narr$TIME_3H), by = "1 hour")
  )
  
  narr_h <- narr[hour_grid, on = "TIME_3H"]
  
  narr_h[, u   := zoo::na.approx(u_3h,   x = as.numeric(TIME_3H), na.rm = FALSE)]
  narr_h[, v   := zoo::na.approx(v_3h,   x = as.numeric(TIME_3H), na.rm = FALSE)]
  narr_h[, dew := zoo::na.approx(dpt_3h, x = as.numeric(TIME_3H), na.rm = FALSE)]
  narr_h[, rh  := zoo::na.approx(rh_3h,  x = as.numeric(TIME_3H), na.rm = FALSE)]
  narr_h[, slp := zoo::na.approx(slp_3h, x = as.numeric(TIME_3H), na.rm = FALSE)]
  
  narr_h[, .(DATE = TIME_3H, u, v, dew, rh, slp)]
}

# -----------------------------
# Paths
# -----------------------------
path_M <- "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/concatenated_files/merged_imputed_final/USW00054788_merged_imputed_2004_2024.csv"
narr_M <- "/Users/rohan/Documents/Rohan/Precip_project/narr_reanalysis_data/USW00054788.csv"

path_B <- "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/concatenated_files/merged_imputed_final/USW00014740_merged_imputed_2004_2024.csv"
narr_B <- "/Users/rohan/Documents/Rohan/Precip_project/narr_reanalysis_data/USW00014740.csv"

outdir <- "/Users/rohan/Documents/Rohan/Precip_project/github/output/gamboostlss/rolling_validation/meriden_bradley_model3"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
persist_outdir <-
  "/Users/rohan/Documents/Rohan/Precip_project/github/output/persistent_forecasts/rolling_validation/meriden_bradley_model3_persistent"

dir.create(persist_outdir, recursive = TRUE, showWarnings = FALSE)
# -----------------------------
# Read data
# -----------------------------
cols_needed <- c(
  "STATION", "DATE", "LATITUDE", "LONGITUDE", "ELEVATION",
  "HourlyPrecipitation", "HourlyRelativeHumidity",
  "HourlyDewPointTemperature", "HourlySeaLevelPressure",
  "HourlyWindSpeed", "HourlyWindDirection"
)

datM <- fread(path_M, select = cols_needed)
datB <- fread(path_B, select = cols_needed)

datM[, DATE := as.POSIXct(DATE, tz = "UTC")]
datB[, DATE := as.POSIXct(DATE, tz = "UTC")]

datM <- datM[
  DATE >= as.POSIXct("2004-01-01 00:00:00", tz = "UTC") &
    DATE <= as.POSIXct("2024-10-09 00:00:00", tz = "UTC")
]

datB <- datB[
  DATE >= as.POSIXct("2004-01-01 00:00:00", tz = "UTC") &
    DATE <= as.POSIXct("2024-10-09 00:00:00", tz = "UTC")
]

setorder(datM, DATE)
setorder(datB, DATE)

DT <- merge(
  datM,
  datB,
  by = "DATE",
  suffixes = c("_M", "_B")
)

setorder(DT, DATE)

# -----------------------------
# Feature construction
# -----------------------------
L <- 6
L_sigma <- 2

DT[, doy := yday(DATE)]

# Meriden
DT[, rain_M := numify(HourlyPrecipitation_M)]
DT[, RH_M   := numify(HourlyRelativeHumidity_M)]
DT[, dew_M  := numify(HourlyDewPointTemperature_M)]
DT[, slp_M  := numify(HourlySeaLevelPressure_M)]

DT[, S_h_M     := numify(HourlyWindSpeed_M)]
DT[, theta_h_M := numify(HourlyWindDirection_M)]

DT[, rad_M := ((630 - theta_h_M) %% 360) * pi / 180]
DT[, WindE_M := S_h_M * cos(rad_M)]
DT[, WindN_M := S_h_M * sin(rad_M)]

# Bradley
DT[, rain_B := numify(HourlyPrecipitation_B)]
DT[, RH_B   := numify(HourlyRelativeHumidity_B)]
DT[, dew_B  := numify(HourlyDewPointTemperature_B)]
DT[, slp_B  := numify(HourlySeaLevelPressure_B)]

DT[, S_h_B     := numify(HourlyWindSpeed_B)]
DT[, theta_h_B := numify(HourlyWindDirection_B)]

DT[, rad_B := ((630 - theta_h_B) %% 360) * pi / 180]
DT[, WindE_B := S_h_B * cos(rad_B)]
DT[, WindN_B := S_h_B * sin(rad_B)]

# -----------------------------
# NARR imputation: Meriden
# -----------------------------
narr_h_M <- narr_to_hourly(narr_M)
DT <- DT[narr_h_M, on = "DATE"]

DT[is.na(WindE_M) & !is.na(u),   WindE_M := u]
DT[is.na(WindN_M) & !is.na(v),   WindN_M := v]
DT[is.na(dew_M)   & !is.na(dew), dew_M   := dew]
DT[is.na(RH_M)    & !is.na(rh),  RH_M    := rh]
DT[is.na(slp_M)   & !is.na(slp), slp_M   := slp]

DT[, c("u", "v", "dew", "rh", "slp") := NULL]

# -----------------------------
# NARR imputation: Bradley
# -----------------------------
narr_h_B <- narr_to_hourly(narr_B)

setnames(
  narr_h_B,
  c("u", "v", "dew", "rh", "slp"),
  c("u_B", "v_B", "dew_B_narr", "rh_B_narr", "slp_B_narr")
)

DT <- DT[narr_h_B, on = "DATE"]

DT[is.na(WindE_B) & !is.na(u_B),        WindE_B := u_B]
DT[is.na(WindN_B) & !is.na(v_B),        WindN_B := v_B]
DT[is.na(dew_B)   & !is.na(dew_B_narr), dew_B   := dew_B_narr]
DT[is.na(RH_B)    & !is.na(rh_B_narr),  RH_B    := rh_B_narr]
DT[is.na(slp_B)   & !is.na(slp_B_narr), slp_B   := slp_B_narr]

DT[, c("u_B", "v_B", "dew_B_narr", "rh_B_narr", "slp_B_narr") := NULL]

# Interactions
DT[, dew_WindE_M := dew_M * WindE_M]
DT[, dew_WindN_M := dew_M * WindN_M]

DT[, dew_WindE_B := dew_B * WindE_B]
DT[, dew_WindN_B := dew_B * WindN_B]

# Lags
for (i in 1:L) {
  
  DT[, paste0("rain_M_lag", i) := shift(rain_M, n = i, type = "lag")]
  DT[, paste0("rain_B_lag", i) := shift(rain_B, n = i, type = "lag")]
  
  DT[, paste0("RH_M_lag", i) := shift(RH_M, n = i, type = "lag")]
  DT[, paste0("dew_M_lag", i) := shift(dew_M, n = i, type = "lag")]
  DT[, paste0("slp_M_lag", i) := shift(slp_M, n = i, type = "lag")]
  DT[, paste0("WindE_M_lag", i) := shift(WindE_M, n = i, type = "lag")]
  DT[, paste0("WindN_M_lag", i) := shift(WindN_M, n = i, type = "lag")]
  DT[, paste0("dew_WindE_M_lag", i) := shift(dew_WindE_M, n = i, type = "lag")]
  DT[, paste0("dew_WindN_M_lag", i) := shift(dew_WindN_M, n = i, type = "lag")]
  
  DT[, paste0("RH_B_lag", i) := shift(RH_B, n = i, type = "lag")]
  DT[, paste0("dew_B_lag", i) := shift(dew_B, n = i, type = "lag")]
  DT[, paste0("slp_B_lag", i) := shift(slp_B, n = i, type = "lag")]
  DT[, paste0("WindE_B_lag", i) := shift(WindE_B, n = i, type = "lag")]
  DT[, paste0("WindN_B_lag", i) := shift(WindN_B, n = i, type = "lag")]
  DT[, paste0("dew_WindE_B_lag", i) := shift(dew_WindE_B, n = i, type = "lag")]
  DT[, paste0("dew_WindN_B_lag", i) := shift(dew_WindN_B, n = i, type = "lag")]
}

for (i in 1:L) {
  DT[, paste0("rain_M_lag", i) := log1p(get(paste0("rain_M_lag", i)))]
  DT[, paste0("rain_B_lag", i) := log1p(get(paste0("rain_B_lag", i)))]
}

DT[, rain_B_log1p := log1p(rain_B)]

vars <- c(
  "DATE", "doy", "rain_M",
  
  paste0("rain_M_lag", 1:L),
  "RH_M", paste0("RH_M_lag", 1:L),
  "dew_M", paste0("dew_M_lag", 1:L),
  "slp_M", paste0("slp_M_lag", 1:L),
  "WindE_M", paste0("WindE_M_lag", 1:L),
  "WindN_M", paste0("WindN_M_lag", 1:L),
  "dew_WindE_M", paste0("dew_WindE_M_lag", 1:L),
  "dew_WindN_M", paste0("dew_WindN_M_lag", 1:L),
  
  "rain_B_log1p",
  paste0("rain_B_lag", 1:L),
  "RH_B", paste0("RH_B_lag", 1:L),
  "dew_B", paste0("dew_B_lag", 1:L),
  "slp_B", paste0("slp_B_lag", 1:L),
  "WindE_B", paste0("WindE_B_lag", 1:L),
  "WindN_B", paste0("WindN_B_lag", 1:L),
  "dew_WindE_B", paste0("dew_WindE_B_lag", 1:L),
  "dew_WindN_B", paste0("dew_WindN_B_lag", 1:L)
)

DT <- DT[, ..vars]
dim(DT)
DT <- na.omit(DT)
dim(DT)
# -----------------------------
# Rolling validation settings
# -----------------------------
train_years <- 4
horizon_hours <- 120
step_hours <- 1

validation_start <- as.POSIXct("2020-12-01 00:00:00", tz = "UTC")
validation_end   <- as.POSIXct("2021-12-01 00:00:00", tz = "UTC")

forecast_origins <- seq(
  from = validation_start,
  to = validation_end - hours(horizon_hours),
  by = paste(step_hours, "hours")
)

cat("Number of forecast origins:", length(forecast_origins), "\n")

# -----------------------------
# Formula
# -----------------------------
df_bbs <- 3

bbs_term <- function(v) {
  paste0("bbs(", v, ", df = ", df_bbs, ")")
}

bols_term <- function(v) {
  paste0("bols(", v, ", intercept = TRUE)")
}

meriden_exo <- c(
  "RH_M", "dew_M", "slp_M",
  "WindE_M", "WindN_M",
  "dew_WindE_M", "dew_WindN_M"
)

bradley_exo <- c(
  "RH_B", "dew_B", "slp_B",
  "WindE_B", "WindN_B",
  "dew_WindE_B", "dew_WindN_B"
)

rhs_mu <- paste(
  c(
    bbs_term("doy"),
    
    vapply(meriden_exo, bbs_term, character(1)),
    vapply(bradley_exo, bbs_term, character(1)),
    bbs_term("rain_B_log1p"),
    
    bols_term("rain_M_lag1"),
    bbs_term("rain_M_lag2"),
    bbs_term("rain_M_lag3"),
    bbs_term("rain_M_lag4"),
    bbs_term("rain_M_lag5"),
    bbs_term("rain_M_lag6"),
    
    bols_term("rain_B_lag1"),
    bbs_term("rain_B_lag2"),
    bbs_term("rain_B_lag3"),
    bbs_term("rain_B_lag4"),
    bbs_term("rain_B_lag5"),
    bbs_term("rain_B_lag6"),
    
    unlist(lapply(1:L, function(h) {
      vapply(paste0(meriden_exo, "_lag", h), bbs_term, character(1))
    })),
    
    unlist(lapply(1:L, function(h) {
      vapply(paste0(bradley_exo, "_lag", h), bbs_term, character(1))
    }))
  ),
  collapse = " + "
)

rhs_sigma <- paste(
  c(
    bbs_term("doy"),
    bbs_term("rain_B_log1p"),
    
    bols_term("rain_M_lag1"),
    bbs_term("rain_M_lag2"),
    
    bbs_term("rain_B_lag1"),
    bbs_term("rain_B_lag2"),
    
    unlist(lapply(1:L_sigma, function(h) {
      vapply(paste0(meriden_exo, "_lag", h), bbs_term, character(1))
    })),
    
    unlist(lapply(1:L_sigma, function(h) {
      vapply(paste0(bradley_exo, "_lag", h), bbs_term, character(1))
    }))
  ),
  collapse = " + "
)

rhs_nu <- paste(
  c(
    bbs_term("doy"),
    bbs_term("rain_M_lag1"),
    bbs_term("rain_B_lag1")
  ),
  collapse = " + "
)

form_list <- list(
  mu    = as.formula(paste("rain_M ~", rhs_mu)),
  sigma = as.formula(paste("rain_M ~", rhs_sigma)),
  nu    = as.formula(paste("rain_M ~", rhs_nu))
)

# -----------------------------
# Storage
# -----------------------------
all_metrics <- list()
all_predictions <- list()
all_runtime <- list()
all_errors <- list()

# -----------------------------
# Rolling validation loop - PARALLEL VERSION
# -----------------------------
library(doSNOW)

t0_total <- Sys.time()

ncores <- detectCores()-2

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
    
    naive_te <- expm1(test_dt$rain_M_lag1)
    
    persist_k <- data.table(
      origin_id = k,
      origin = origin,
      DATE = test_dt$DATE,
      actual = Y_te,
      persistence = naive_te
    )
    
    naive_tr_actual <- train_dt$rain_M
    naive_tr_pred   <- expm1(train_dt$rain_M_lag1)
    
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
      qmape = qmape(Y_te, pred_te),
      smape = smape(Y_te, pred_te),
      mase = mase(
        actual = Y_te,
        pred = pred_te,
        naive_actual = naive_tr_actual,
        naive_pred = naive_tr_pred
      ),
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
      naive = naive_te,
      lo80 = lo80_te,
      hi80 = hi80_te,
      lo95 = lo95_te,
      hi95 = hi95_te
    )
    
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

metrics <- rbindlist(all_metrics, fill = TRUE)
predictions <- rbindlist(all_predictions, fill = TRUE)
runtime <- rbindlist(all_runtime, fill = TRUE)
persistence <- rbindlist(all_persistence, fill = TRUE)

# -----------------------------
# Save main outputs
# -----------------------------
fwrite(metrics, file.path(outdir, "rolling_metrics.csv"))
fwrite(predictions, file.path(outdir, "rolling_predictions.csv"))
fwrite(runtime, file.path(outdir, "runtime.csv"))

fwrite(
  persistence,
  file.path(persist_outdir, "persistent_predictions.csv")
)

if (length(all_errors) > 0) {
  errors <- rbindlist(all_errors, fill = TRUE)
  
  if (nrow(errors) > 0) {
    fwrite(errors, file.path(outdir, "errors.csv"))
  }
}

# -----------------------------
# Hurricane / Calm summary
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
    mean_qmape = mean(qmape, na.rm = TRUE),
    mean_smape = mean(smape, na.rm = TRUE),
    mean_mase = mean(mase, na.rm = TRUE),
    mean_peak90_rmse = mean(peak90_rmse, na.rm = TRUE),
    mean_peak95_rmse = mean(peak95_rmse, na.rm = TRUE),
    mean_peak99_rmse = mean(peak99_rmse, na.rm = TRUE)
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
    
    mean_qmape = mean(qmape, na.rm = TRUE),
    sd_qmape = sd(qmape, na.rm = TRUE),
    median_qmape = median(qmape, na.rm = TRUE),
    
    mean_smape = mean(smape, na.rm = TRUE),
    sd_smape = sd(smape, na.rm = TRUE),
    median_smape = median(smape, na.rm = TRUE),
    
    mean_mase = mean(mase, na.rm = TRUE),
    sd_mase = sd(mase, na.rm = TRUE),
    median_mase = median(mase, na.rm = TRUE),
    
    mean_peak90_rmse = mean(peak90_rmse, na.rm = TRUE),
    mean_peak95_rmse = mean(peak95_rmse, na.rm = TRUE),
    mean_peak99_rmse = mean(peak99_rmse, na.rm = TRUE)
  )
]

fwrite(
  overall_summary,
  file.path(outdir, "overall_summary.csv")
)

# -----------------------------
# Bucket-wise summary by observed rainfall
# -----------------------------
bucket_summary <- predictions[
  ,
  .(
    n = .N,
    rmse = rmse(actual, zaga),
    mae = mae(actual, zaga),
    bias = bias(actual, zaga),
    qmape = qmape(actual, zaga),
    smape = smape(actual, zaga),
    naive_mae = mae(actual, naive),
    mase = mae(actual, zaga) / mae(actual, naive)
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