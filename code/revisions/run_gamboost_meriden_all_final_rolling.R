# ============================================================
# ROLLING-WINDOW VALIDATION
# ZAGA gamboostLSS + NARR imputation
# MODEL: Meriden ~ Meriden + ALL OTHER CT STATIONS
# ============================================================

userlib <- "/home/fbs24003/Rlibs/r452"
dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(userlib, .libPaths()))

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

for (p in req) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, lib = userlib, repos = "https://cloud.r-project.org")
  }
}

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

# ============================================================
# Helpers
# ============================================================

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

qmape <- function(actual, pred) {
  denom <- sum(actual, na.rm = TRUE)
  if (is.na(denom) || denom == 0) return(NA_real_)
  sum(abs(actual - pred), na.rm = TRUE) / denom
}

smape <- function(actual, pred) {
  denom <- (abs(actual) + abs(pred)) / 2
  idx <- denom > 0
  if (sum(idx) == 0) return(NA_real_)
  100 * mean(abs(actual[idx] - pred[idx]) / denom[idx], na.rm = TRUE)
}

mase <- function(actual, pred, naive_actual, naive_pred) {
  model_mae <- mae(actual, pred)
  naive_mae <- mae(naive_actual, naive_pred)
  if (is.na(naive_mae) || naive_mae == 0) return(NA_real_)
  model_mae / naive_mae
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

# ============================================================
# NARR helper
# ============================================================

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

# ============================================================
# Paths and stations
# ============================================================

base_data_dir <- "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/concatenated_files/merged_imputed_final"
base_narr_dir <- "/Users/rohan/Documents/Rohan/Precip_project/narr_reanalysis_data"

station_map <- data.table(
  code = c("M", "B", "G", "N", "D", "W", "S"),
  id   = c(
    "USW00054788",  # Meriden
    "USW00014740",  # Bradley
    "USW00014707",  # Groton
    "USW00014758",  # New Haven
    "USW00054734",  # Danbury
    "USW00054767",  # Windham
    "USW00094702"   # Sikorsky
  )
)

target_code <- "M"
other_codes <- setdiff(station_map$code, target_code)

outdir <- "/Users/rohan/Documents/Rohan/Precip_project/github/output/gamboostlss/rolling_validation/meriden_allstations_model"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

persist_outdir <- "/Users/rohan/Documents/Rohan/Precip_project/github/output/persistent_forecasts/rolling_validation/meriden_allstations_persistent"
dir.create(persist_outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Read one station and rename variables
# ============================================================

read_station_data <- function(station_id, code) {
  
  station_path <- file.path(
    base_data_dir,
    paste0(station_id, "_merged_imputed_2004_2024.csv")
  )
  
  narr_path <- file.path(
    base_narr_dir,
    paste0(station_id, ".csv")
  )
  
  cols_needed <- c(
    "STATION",
    "DATE",
    "HourlyPrecipitation",
    "HourlyRelativeHumidity",
    "HourlyDewPointTemperature",
    "HourlySeaLevelPressure",
    "HourlyWindSpeed",
    "HourlyWindDirection"
  )
  
  dat <- fread(station_path, select = cols_needed)
  
  dat[, DATE := as.POSIXct(DATE, tz = "UTC")]
  setorder(dat, DATE)
  
  dat <- dat[
    DATE >= as.POSIXct("2004-01-01 00:00:00", tz = "UTC") &
      DATE <= as.POSIXct("2024-10-09 00:00:00", tz = "UTC")
  ]
  
  dat[, paste0("rain_", code) := numify(HourlyPrecipitation)]
  dat[, paste0("RH_", code)   := numify(HourlyRelativeHumidity)]
  dat[, paste0("dew_", code)  := numify(HourlyDewPointTemperature)]
  dat[, paste0("slp_", code)  := numify(HourlySeaLevelPressure)]
  
  dat[, paste0("S_h_", code)     := numify(HourlyWindSpeed)]
  dat[, paste0("theta_h_", code) := numify(HourlyWindDirection)]
  
  dat[, paste0("rad_", code) :=
        ((630 - get(paste0("theta_h_", code))) %% 360) * pi / 180]
  
  dat[, paste0("WindE_", code) :=
        get(paste0("S_h_", code)) * cos(get(paste0("rad_", code)))]
  
  dat[, paste0("WindN_", code) :=
        get(paste0("S_h_", code)) * sin(get(paste0("rad_", code)))]
  
  # NARR imputation
  narr_h <- narr_to_hourly(narr_path)
  
  dat <- dat[narr_h, on = "DATE"]
  
  dat[is.na(get(paste0("WindE_", code))) & !is.na(u),
      paste0("WindE_", code) := u]
  
  dat[is.na(get(paste0("WindN_", code))) & !is.na(v),
      paste0("WindN_", code) := v]
  
  dat[is.na(get(paste0("dew_", code))) & !is.na(dew),
      paste0("dew_", code) := dew]
  
  dat[is.na(get(paste0("RH_", code))) & !is.na(rh),
      paste0("RH_", code) := rh]
  
  dat[is.na(get(paste0("slp_", code))) & !is.na(slp),
      paste0("slp_", code) := slp]
  
  dat[, c("u", "v", "dew", "rh", "slp") := NULL]
  
  dat[, paste0("dew_WindE_", code) :=
        get(paste0("dew_", code)) * get(paste0("WindE_", code))]
  
  dat[, paste0("dew_WindN_", code) :=
        get(paste0("dew_", code)) * get(paste0("WindN_", code))]
  
  keep <- c(
    "DATE",
    paste0("rain_", code),
    paste0("RH_", code),
    paste0("dew_", code),
    paste0("slp_", code),
    paste0("WindE_", code),
    paste0("WindN_", code),
    paste0("dew_WindE_", code),
    paste0("dew_WindN_", code)
  )
  
  dat[, ..keep]
}

# ============================================================
# Read and merge all stations
# ============================================================

station_data_list <- list()

for (r in seq_len(nrow(station_map))) {
  code_r <- station_map$code[r]
  id_r   <- station_map$id[r]
  
  cat("Reading station:", code_r, id_r, "\n")
  
  station_data_list[[code_r]] <- read_station_data(
    station_id = id_r,
    code = code_r
  )
}

DT <- Reduce(
  function(x, y) merge(x, y, by = "DATE", all = FALSE),
  station_data_list
)

setorder(DT, DATE)

DT[, doy := yday(DATE)]

# ============================================================
# Feature construction
# ============================================================

L <- 6
L_sigma <- 2

all_codes <- station_map$code

for (code in all_codes) {
  
  rain_var <- paste0("rain_", code)
  
  for (h in 1:L) {
    DT[, paste0(rain_var, "_lag", h) := shift(get(rain_var), n = h, type = "lag")]
    DT[, paste0(rain_var, "_lag", h) := log1p(get(paste0(rain_var, "_lag", h)))]
  }
  
  if (code != target_code) {
    DT[, paste0(rain_var, "_log1p") := log1p(get(rain_var))]
  }
}

exo_base <- c("RH", "dew", "slp", "WindE", "WindN", "dew_WindE", "dew_WindN")

for (code in all_codes) {
  for (v in exo_base) {
    base_var <- paste0(v, "_", code)
    
    for (h in 1:L) {
      DT[, paste0(base_var, "_lag", h) := shift(get(base_var), n = h, type = "lag")]
    }
  }
}

target_rain_lags <- paste0("rain_", target_code, "_lag", 1:L)

target_exo_current <- paste0(exo_base, "_", target_code)

target_exo_lags <- unlist(
  lapply(1:L, function(h) {
    paste0(target_exo_current, "_lag", h)
  })
)

other_rain_current <- paste0("rain_", other_codes, "_log1p")

other_rain_lags <- unlist(
  lapply(other_codes, function(code) {
    paste0("rain_", code, "_lag", 1:L)
  })
)

other_exo_current <- unlist(
  lapply(other_codes, function(code) {
    paste0(exo_base, "_", code)
  })
)

other_exo_lags <- unlist(
  lapply(other_codes, function(code) {
    vars_code <- paste0(exo_base, "_", code)
    unlist(
      lapply(1:L, function(h) {
        paste0(vars_code, "_lag", h)
      })
    )
  })
)

vars <- c(
  "DATE",
  "doy",
  paste0("rain_", target_code),
  target_rain_lags,
  target_exo_current,
  target_exo_lags,
  other_rain_current,
  other_rain_lags,
  other_exo_current,
  other_exo_lags
)

DT <- DT[, ..vars]
DT <- na.omit(DT)

setnames(DT, paste0("rain_", target_code), "rain_M")

cat("Final modelling dataset dimensions:\n")
print(dim(DT))

# ============================================================
# Rolling validation settings
# ============================================================

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

# ============================================================
# Formula construction
# ============================================================

df_bbs <- 3

bbs_term <- function(v) {
  paste0("bbs(", v, ", df = ", df_bbs, ")")
}

bols_term <- function(v) {
  paste0("bols(", v, ", intercept = TRUE)")
}

target_rain_lags_model <- target_rain_lags

target_exo_current_model <- target_exo_current
target_exo_lags_model <- target_exo_lags

other_rain_current_model <- other_rain_current
other_rain_lags_model <- other_rain_lags
other_exo_current_model <- other_exo_current
other_exo_lags_model <- other_exo_lags

# -----------------------------
# mu formula
# -----------------------------
rhs_mu <- paste(
  c(
    bbs_term("doy"),
    
    vapply(target_exo_current_model, bbs_term, character(1)),
    vapply(other_exo_current_model, bbs_term, character(1)),
    vapply(other_rain_current_model, bbs_term, character(1)),
    
    bols_term(target_rain_lags_model[1]),
    vapply(target_rain_lags_model[-1], bbs_term, character(1)),
    
    unlist(
      lapply(other_codes, function(code) {
        lag_vars <- paste0("rain_", code, "_lag", 1:L)
        c(
          bols_term(lag_vars[1]),
          vapply(lag_vars[-1], bbs_term, character(1))
        )
      })
    ),
    
    vapply(target_exo_lags_model, bbs_term, character(1)),
    vapply(other_exo_lags_model, bbs_term, character(1))
  ),
  collapse = " + "
)

# -----------------------------
# sigma formula
# -----------------------------
target_rain_lags_sigma <- paste0("rain_", target_code, "_lag", 1:L_sigma)

other_rain_lags_sigma <- unlist(
  lapply(other_codes, function(code) {
    paste0("rain_", code, "_lag", 1:L_sigma)
  })
)

target_exo_lags_sigma <- unlist(
  lapply(1:L_sigma, function(h) {
    paste0(target_exo_current, "_lag", h)
  })
)

other_exo_lags_sigma <- unlist(
  lapply(other_codes, function(code) {
    vars_code <- paste0(exo_base, "_", code)
    unlist(
      lapply(1:L_sigma, function(h) {
        paste0(vars_code, "_lag", h)
      })
    )
  })
)

rhs_sigma <- paste(
  c(
    bbs_term("doy"),
    
    vapply(other_rain_current_model, bbs_term, character(1)),
    
    bols_term(target_rain_lags_sigma[1]),
    vapply(target_rain_lags_sigma[-1], bbs_term, character(1)),
    
    unlist(
      lapply(other_codes, function(code) {
        lag_vars <- paste0("rain_", code, "_lag", 1:L_sigma)
        c(
          bols_term(lag_vars[1]),
          vapply(lag_vars[-1], bbs_term, character(1))
        )
      })
    ),
    
    vapply(target_exo_lags_sigma, bbs_term, character(1)),
    vapply(other_exo_lags_sigma, bbs_term, character(1))
  ),
  collapse = " + "
)

# -----------------------------
# nu formula
# -----------------------------
other_lag1 <- paste0("rain_", other_codes, "_lag1")

rhs_nu <- paste(
  c(
    bbs_term("doy"),
    bbs_term(paste0("rain_", target_code, "_lag1")),
    vapply(other_lag1, bbs_term, character(1))
  ),
  collapse = " + "
)

form_list <- list(
  mu    = as.formula(paste("rain_M ~", rhs_mu)),
  sigma = as.formula(paste("rain_M ~", rhs_sigma)),
  nu    = as.formula(paste("rain_M ~", rhs_nu))
)

cat("\nmu formula length:", nchar(rhs_mu), "\n")
cat("sigma formula length:", nchar(rhs_sigma), "\n")
cat("nu formula length:", nchar(rhs_nu), "\n")

# ============================================================
# Parallel rolling validation
# ============================================================

t0_total <- Sys.time()

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "10"))
ncores <- max(1, ncores)

cl <- parallel::makeCluster(ncores)

parallel::clusterExport(cl, "userlib", envir = environment())

parallel::clusterEvalQ(cl, {
  .libPaths(c(userlib, .libPaths()))
  library(data.table)
  library(lubridate)
  library(gamboostLSS)
  library(mboost)
  library(gamlss.dist)
  library(zoo)
})
parallel::clusterExport(cl, "DT", envir = environment())

doSNOW::registerDoSNOW(cl)

cat("Using", ncores, "cores\n")
#forecast_origins = forecast_origins[1000:1005]
pb <- txtProgressBar(max = length(forecast_origins), style = 3)

progress <- function(n) {
  setTxtProgressBar(pb, n)
}

opts <- list(progress = progress)

results <- foreach(
  k = seq_along(forecast_origins),
  .options.snow = opts,
  .packages = character(0),
  .export = c(
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
    
    Y_te <- test_dt$rain_M
    
    naive_te <- expm1(test_dt[[paste0("rain_", target_code, "_lag1")]])
    
    persist_k <- data.table(
      origin_id = k,
      origin = origin,
      DATE = test_dt$DATE,
      actual = Y_te,
      persistence = naive_te
    )
    
    naive_tr_actual <- train_dt$rain_M
    naive_tr_pred   <- expm1(train_dt[[paste0("rain_", target_code, "_lag1")]])
    
    fit.gb <- gamboostLSS(
      formula  = form_list,
      data     = dtr,
      families = as.families("ZAGA"),
      control  = boost_control(mstop = 50, nu = 0.01)
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

# ============================================================
# Collect outputs
# ============================================================

all_metrics     <- lapply(results, `[[`, "metrics")
all_predictions <- lapply(results, `[[`, "predictions")
all_runtime     <- lapply(results, `[[`, "runtime")
all_errors      <- lapply(results, `[[`, "error")
all_persistence <- lapply(results, `[[`, "persistence")

metrics <- rbindlist(all_metrics, fill = TRUE)
predictions <- rbindlist(all_predictions, fill = TRUE)
runtime <- rbindlist(all_runtime, fill = TRUE)
persistence <- rbindlist(all_persistence, fill = TRUE)

if (nrow(metrics) == 0) {
  errors <- rbindlist(all_errors, fill = TRUE)
  fwrite(errors, file.path(outdir, "errors.csv"))
  stop("All all-station gamboostLSS fits failed. Check errors.csv.")
}

# ============================================================
# Save outputs
# ============================================================

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

# ============================================================
# Seasonal summary
# ============================================================

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

# ============================================================
# Overall summary
# ============================================================

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

# ============================================================
# Bucket-wise summary by observed rainfall
# ============================================================

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

# ============================================================
# Total runtime
# ============================================================

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