# ============================================================
# HPC FAILSAFE (GAMBOOSTLSS)
# ZAGA (gamboostLSS) + NARR imputation (3-hourly -> hourly)
# MODEL: Groton ~ Groton (lags<=6) + Windham (current + lags<=6)
# Save: model, metrics, predictions, plot PDF, runtime CSVs
# ============================================================

# -----------------------------
# Personal library + installs
# -----------------------------
userlib <- "/home/fbs24003/Rlibs/r452"
dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(userlib, .libPaths()))

req <- c(
  "data.table",
  "lubridate",
  "gamboostLSS",
  "mboost",
  "gamlss.dist",
  "ggplot2",
  "zoo"
)

for (p in req) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing missing package: ", p)
    install.packages(p, lib = userlib, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(gamboostLSS)
  library(mboost)
  library(gamlss.dist)
  library(ggplot2)
  library(zoo)
})

# -----------------------------
# Helpers
# -----------------------------
rmse <- function(a, p) sqrt(mean((a - p)^2, na.rm = TRUE))

rmse_peaks <- function(a, p, thr = NULL, q = 0.95) {
  if (is.null(thr)) thr <- as.numeric(quantile(a, probs = q, na.rm = TRUE))
  idx <- which(a >= thr)
  if (length(idx) < 5) return(list(threshold = thr, n = length(idx), rmse = NA_real_))
  list(threshold = thr, n = length(idx), rmse = rmse(a[idx], p[idx]))
}

numify <- function(x) {
  if (is.character(x)) x[x == ""] <- NA_character_
  as.numeric(x)
}

# -----------------------------
# NARR imputation helper (3-hourly -> hourly) for u/v/dew/rh/slp
# Returns data.table(DATE, u, v, dew, rh, slp)
# -----------------------------
narr_to_hourly <- function(narr_path) {
  
  narr <- fread(narr_path)
  
  # Detect time column
  if ("time" %in% names(narr)) {
    narr[, time := as.POSIXct(time, tz = "UTC")]
    setnames(narr, "time", "TIME_3H")
  } else if ("DATE" %in% names(narr)) {
    narr[, DATE := as.POSIXct(DATE, tz = "UTC")]
    setnames(narr, "DATE", "TIME_3H")
  } else {
    stop("NARR file has no 'time' or 'DATE' column.")
  }
  
  # Detect u/v columns
  if (all(c("uwnd","vwnd") %in% names(narr))) {
    narr[, `:=`(u_3h = numify(uwnd), v_3h = numify(vwnd))]
  } else if (all(c("u_3h","v_3h") %in% names(narr))) {
    narr[, `:=`(u_3h = numify(u_3h), v_3h = numify(v_3h))]
  } else {
    uv_cols <- tail(names(narr), 2)
    narr[, `:=`(u_3h = numify(get(uv_cols[1])), v_3h = numify(get(uv_cols[2])))]
  }
  
  # Detect dewpoint column
  if ("dpt" %in% names(narr)) {
    narr[, dpt_3h := numify(dpt)]
  } else if ("dpt_narr" %in% names(narr)) {
    narr[, dpt_3h := numify(dpt_narr)]
  } else {
    narr[, dpt_3h := NA_real_]
  }
  
  # Detect RH column (rhum)
  if ("rhum" %in% names(narr)) {
    narr[, rh_3h := numify(rhum)]
  } else if ("RH" %in% names(narr)) {
    narr[, rh_3h := numify(RH)]
  } else if ("rh_3h" %in% names(narr)) {
    narr[, rh_3h := numify(rh_3h)]
  } else {
    narr[, rh_3h := NA_real_]
  }
  
  # Detect sea-level pressure column (prmsl)
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
  
  # Kelvin -> Celsius if needed
  if (!all(is.na(narr$dpt_3h))) {
    if (median(narr$dpt_3h, na.rm = TRUE) > 100) narr[, dpt_3h := dpt_3h - 273.15]
  }
  
  # Pa -> hPa if needed
  if (!all(is.na(narr$slp_3h))) {
    if (median(narr$slp_3h, na.rm = TRUE) > 2000) narr[, slp_3h := slp_3h / 100]
  }
  
  # Hourly grid
  hour_grid <- data.table(TIME_3H = seq(min(narr$TIME_3H), max(narr$TIME_3H), by = "1 hour"))
  narr_h <- narr[hour_grid, on = "TIME_3H"]
  
  narr_h[, u   := zoo::na.approx(u_3h,   x = as.numeric(TIME_3H), na.rm = FALSE)]
  narr_h[, v   := zoo::na.approx(v_3h,   x = as.numeric(TIME_3H), na.rm = FALSE)]
  narr_h[, dew := zoo::na.approx(dpt_3h, x = as.numeric(TIME_3H), na.rm = FALSE)]
  narr_h[, rh  := zoo::na.approx(rh_3h,  x = as.numeric(TIME_3H), na.rm = FALSE)]
  narr_h[, slp := zoo::na.approx(slp_3h, x = as.numeric(TIME_3H), na.rm = FALSE)]
  
  narr_h[, .(DATE = TIME_3H, u, v, dew, rh, slp)]
}

# -----------------------------
# Runtime start
# -----------------------------
t0_total <- Sys.time()

# -----------------------------
# Paths
# -----------------------------
# Groton station file + NARR
path_G   <- "/home/fbs24003/rlibs/Precip_project/nalini_ct_airport_met/concatenated_files/merged_imputed_final/USW00014707_merged_imputed_2004_2024.csv"
narr_G   <- "/home/fbs24003/rlibs/Precip_project/narr_reanalysis_data/USW00014707.csv"

# Windham station file + NARR (USW00054767)
path_W   <- "/home/fbs24003/rlibs/Precip_project/nalini_ct_airport_met/concatenated_files/merged_imputed_final/USW00054767_merged_imputed_2004_2024.csv"
narr_W   <- "/home/fbs24003/rlibs/Precip_project/narr_reanalysis_data/USW00054767.csv"

outdir <- "/home/fbs24003/rlibs/Precip_project/nalini_ct_airport_met/gam_lss paper/out_zaga_gamboostlss_groton_windham"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Read + base feature build
# -----------------------------
cols_needed <- c(
  "STATION","DATE","LATITUDE","LONGITUDE","ELEVATION",
  "HourlyPrecipitation","HourlyRelativeHumidity",
  "HourlyDewPointTemperature","HourlySeaLevelPressure",
  "HourlyWindSpeed","HourlyWindDirection"
)

datG <- fread(path_G, select = cols_needed)
datW <- fread(path_W, select = cols_needed)

datG[, DATE := as.POSIXct(DATE, tz = "UTC")]
datW[, DATE := as.POSIXct(DATE, tz = "UTC")]

datG <- datG[DATE >= as.POSIXct("2004-01-01 00:00:00", tz = "UTC")]
datW <- datW[DATE >= as.POSIXct("2004-01-01 00:00:00", tz = "UTC")]

setorder(datG, DATE)
setorder(datW, DATE)

# Merge on common timeline
DT <- merge(
  datG, datW,
  by = "DATE",
  suffixes = c("_G", "_W")
)
setorder(DT, DATE)

# number of lags
L <- 6

# Time features
DT[, hour := hour(DATE)]
DT[, doy  := yday(DATE)]

# -----------------------------
# Groton core vars
# -----------------------------
DT[, rain_G := numify(HourlyPrecipitation_G)]
DT[, RH_G   := numify(HourlyRelativeHumidity_G)]
DT[, dew_G  := numify(HourlyDewPointTemperature_G)]
DT[, slp_G  := numify(HourlySeaLevelPressure_G)]

DT[, `:=`(
  alpha_tilde_G = 0,
  S_h_G     = numify(HourlyWindSpeed_G),
  theta_h_G = numify(HourlyWindDirection_G)
)]
DT[, rad_G := ((alpha_tilde_G + 630 - theta_h_G) %% 360) * pi/180]
DT[, `:=`(
  WindE_G = S_h_G * cos(rad_G),
  WindN_G = S_h_G * sin(rad_G)
)]

# -----------------------------
# Windham core vars
# -----------------------------
DT[, rain_W := numify(HourlyPrecipitation_W)]
DT[, RH_W   := numify(HourlyRelativeHumidity_W)]
DT[, dew_W  := numify(HourlyDewPointTemperature_W)]
DT[, slp_W  := numify(HourlySeaLevelPressure_W)]

DT[, `:=`(
  alpha_tilde_W = 0,
  S_h_W     = numify(HourlyWindSpeed_W),
  theta_h_W = numify(HourlyWindDirection_W)
)]
DT[, rad_W := ((alpha_tilde_W + 630 - theta_h_W) %% 360) * pi/180]
DT[, `:=`(
  WindE_W = S_h_W * cos(rad_W),
  WindN_W = S_h_W * sin(rad_W)
)]

# --------------------------------
# NARR imputation (Groton)
# --------------------------------
narr_h_G <- narr_to_hourly(narr_G)
DT <- DT[narr_h_G, on = "DATE"]

DT[is.na(WindE_G) & !is.na(u),   WindE_G := u]
DT[is.na(WindN_G) & !is.na(v),   WindN_G := v]
DT[is.na(dew_G)  & !is.na(dew),  dew_G  := dew]
DT[is.na(RH_G)   & !is.na(rh),   RH_G   := rh]
DT[is.na(slp_G)  & !is.na(slp),  slp_G  := slp]

DT[, c("u","v","dew","rh","slp") := NULL]

# --------------------------------
# NARR imputation (Windham)
# --------------------------------
narr_h_W <- narr_to_hourly(narr_W)
setnames(narr_h_W, c("u","v","dew","rh","slp"),
         c("u_W","v_W","dew_W_narr","rh_W_narr","slp_W_narr"))

DT <- DT[narr_h_W, on = "DATE"]

DT[is.na(WindE_W) & !is.na(u_W),        WindE_W := u_W]
DT[is.na(WindN_W) & !is.na(v_W),        WindN_W := v_W]
DT[is.na(dew_W)  & !is.na(dew_W_narr),  dew_W  := dew_W_narr]
DT[is.na(RH_W)   & !is.na(rh_W_narr),   RH_W   := rh_W_narr]
DT[is.na(slp_W)  & !is.na(slp_W_narr),  slp_W  := slp_W_narr]

DT[, c("u_W","v_W","dew_W_narr","rh_W_narr","slp_W_narr") := NULL]

# --------------------------------
# Interaction terms (time t)
# --------------------------------
DT[, `:=`(
  dew_WindE_G = dew_G * WindE_G,
  dew_WindN_G = dew_G * WindN_G,
  dew_WindE_W = dew_W * WindE_W,
  dew_WindN_W = dew_W * WindN_W
)]

# --------------------------------
# Lags (Groton + Windham)
# --------------------------------
for (i in 1:L) {
  # Groton lags
  DT[, paste0("rain_G_lag",    i) := shift(rain_G,    n = i, type = "lag")]
  DT[, paste0("RH_G_lag",      i) := shift(RH_G,      n = i, type = "lag")]
  DT[, paste0("dew_G_lag",     i) := shift(dew_G,     n = i, type = "lag")]
  DT[, paste0("slp_G_lag",     i) := shift(slp_G,     n = i, type = "lag")]
  DT[, paste0("WindE_G_lag",   i) := shift(WindE_G,   n = i, type = "lag")]
  DT[, paste0("WindN_G_lag",   i) := shift(WindN_G,   n = i, type = "lag")]
  DT[, paste0("dew_WindE_G_lag", i) := shift(dew_WindE_G, n = i, type = "lag")]
  DT[, paste0("dew_WindN_G_lag", i) := shift(dew_WindN_G, n = i, type = "lag")]
  
  # Windham lags
  DT[, paste0("rain_W_lag",    i) := shift(rain_W,    n = i, type = "lag")]
  DT[, paste0("RH_W_lag",      i) := shift(RH_W,      n = i, type = "lag")]
  DT[, paste0("dew_W_lag",     i) := shift(dew_W,     n = i, type = "lag")]
  DT[, paste0("slp_W_lag",     i) := shift(slp_W,     n = i, type = "lag")]
  DT[, paste0("WindE_W_lag",   i) := shift(WindE_W,   n = i, type = "lag")]
  DT[, paste0("WindN_W_lag",   i) := shift(WindN_W,   n = i, type = "lag")]
  DT[, paste0("dew_WindE_W_lag", i) := shift(dew_WindE_W, n = i, type = "lag")]
  DT[, paste0("dew_WindN_W_lag", i) := shift(dew_WindN_W, n = i, type = "lag")]
}

# log1p on rainfall lags (both stations) + OPTIONAL current Windham rain
for (i in 1:L) {
  DT[, paste0("rain_G_lag", i) := log1p(get(paste0("rain_G_lag", i)))]
  DT[, paste0("rain_W_lag", i) := log1p(get(paste0("rain_W_lag", i)))]
}
DT[, rain_W_log1p := log1p(rain_W)]

# --------------------------------
# Keep only modeling columns + NA trim
# --------------------------------
vars <- c(
  "DATE","hour","doy","rain_G",
  
  # Groton (current + lags)
  paste0("rain_G_lag", 1:L),
  "RH_G", paste0("RH_G_lag", 1:L),
  "dew_G", paste0("dew_G_lag", 1:L),
  "slp_G", paste0("slp_G_lag", 1:L),
  "WindE_G", paste0("WindE_G_lag", 1:L),
  "WindN_G", paste0("WindN_G_lag", 1:L),
  "dew_WindE_G", paste0("dew_WindE_G_lag", 1:L),
  "dew_WindN_G", paste0("dew_WindN_G_lag", 1:L),
  
  # Windham (current + lags)
  "rain_W_log1p",
  paste0("rain_W_lag", 1:L),
  "RH_W", paste0("RH_W_lag", 1:L),
  "dew_W", paste0("dew_W_lag", 1:L),
  "slp_W", paste0("slp_W_lag", 1:L),
  "WindE_W", paste0("WindE_W_lag", 1:L),
  "WindN_W", paste0("WindN_W_lag", 1:L),
  "dew_WindE_W", paste0("dew_WindE_W_lag", 1:L),
  "dew_WindN_W", paste0("dew_WindN_W_lag", 1:L)
)

DT <- DT[, ..vars]
DT <- na.omit(DT)
DT <- DT[DATE <= as.POSIXct("2024-10-09 00:00:00", tz = "UTC")]

# -----------------------------
# Train/Test split
# -----------------------------
split_date <- as.POSIXct("2024-10-08 13:00:00", tz="UTC")
i_train <- DT$DATE < split_date
i_test  <- DT$DATE >= split_date

pred <- setdiff(names(DT), c("DATE","rain_G"))

train_dt <- DT[i_train, c("DATE","rain_G", pred), with=FALSE]
test_dt  <- DT[i_test,  c("DATE","rain_G", pred), with=FALSE]

train_dt <- na.omit(train_dt)
test_dt  <- na.omit(test_dt)

dtr <- as.data.frame(train_dt[, c("rain_G", pred), with=FALSE])
dte <- as.data.frame(test_dt[,  c("rain_G", pred), with=FALSE])

DATE_tr <- train_dt$DATE; Y_tr <- train_dt$rain_G
DATE_te <- test_dt$DATE;  Y_te <- test_dt$rain_G

# -----------------------------
# GAMBOOSTLSS ZAGA fit (bbs for ALL terms, df=3)
# -----------------------------
gc()
df_bbs <- 3
bbs_term <- function(v) paste0("bbs(", v, ", df=", df_bbs,")")
bols_term_rain <- function(v) paste0("bols(", v,",intercept=TRUE)")

rhs_mu_gb <- paste(
  c(
    bbs_term("doy"),
    
    # --- current-time Groton met terms ---
    bbs_term("RH_G"),
    bbs_term("dew_G"),
    bbs_term("slp_G"),
    bbs_term("WindE_G"),
    bbs_term("WindN_G"),
    bbs_term("dew_WindE_G"),
    bbs_term("dew_WindN_G"),
    
    # --- current-time Windham rain (log1p) + met terms ---
    bbs_term("rain_W_log1p"),
    bbs_term("RH_W"),
    bbs_term("dew_W"),
    bbs_term("slp_W"),
    bbs_term("WindE_W"),
    bbs_term("WindN_W"),
    bbs_term("dew_WindE_W"),
    bbs_term("dew_WindN_W"),
    
    # --- Groton lag terms ---
    vapply(paste0("rain_G_lag", 1:1), bols_term_rain, character(1)),
    vapply(paste0("rain_G_lag", 2:2), bbs_term, character(1)),
    vapply(paste0("rain_G_lag", 3:L), bbs_term, character(1)),
    vapply(paste0("RH_G_lag",  1:L), bbs_term, character(1)),
    vapply(paste0("dew_G_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("slp_G_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("WindE_G_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("WindN_G_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("dew_WindE_G_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("dew_WindN_G_lag", 1:L), bbs_term, character(1)),
    
    # --- Windham lag terms ---
    vapply(paste0("rain_W_lag", 1:1), bols_term_rain, character(1)),
    vapply(paste0("rain_W_lag", 2:L), bbs_term, character(1)),
    vapply(paste0("RH_W_lag",  1:L), bbs_term, character(1)),
    vapply(paste0("dew_W_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("slp_W_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("WindE_W_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("WindN_W_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("dew_WindE_W_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("dew_WindN_W_lag", 1:L), bbs_term, character(1))
  ),
  collapse = " + "
)

L_sigma <- 2
rhs_sigma_gb <- paste(
  c(
    bbs_term("doy"),
    bbs_term("rain_W_log1p"),
    
    vapply(paste0("rain_G_lag", 1:1), bols_term_rain, character(1)),
    vapply(paste0("rain_G_lag", 2:L_sigma), bbs_term, character(1)),
    vapply(paste0("RH_G_lag",  1:L_sigma), bbs_term, character(1)),
    vapply(paste0("dew_G_lag", 1:L_sigma), bbs_term, character(1)),
    vapply(paste0("slp_G_lag", 1:L_sigma), bbs_term, character(1)),
    vapply(paste0("WindE_G_lag", 1:L_sigma), bbs_term, character(1)),
    vapply(paste0("WindN_G_lag", 1:L_sigma), bbs_term, character(1)),
    vapply(paste0("dew_WindE_G_lag", 1:L_sigma), bbs_term, character(1)),
    vapply(paste0("dew_WindN_G_lag", 1:L_sigma), bbs_term, character(1)),
    
    # --- Windham lagged terms ---
    vapply(paste0("rain_W_lag", 1:L_sigma), bbs_term, character(1)),
    vapply(paste0("RH_W_lag",  1:L), bbs_term, character(1)),
    vapply(paste0("dew_W_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("slp_W_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("WindE_W_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("WindN_W_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("dew_WindE_W_lag", 1:L), bbs_term, character(1)),
    vapply(paste0("dew_WindN_W_lag", 1:L), bbs_term, character(1))
  ),
  collapse = " + "
)

rhs_nu_gb <- paste(
  c(
    bbs_term("doy"),
    bbs_term("rain_G_lag1"),
    bbs_term("rain_W_lag1")
  ),
  collapse = " + "
)

form_list <- list(
  mu    = as.formula(paste("rain_G ~", rhs_mu_gb)),
  sigma = as.formula(paste("rain_G ~", rhs_sigma_gb)),
  nu    = as.formula(paste("rain_G ~", rhs_nu_gb))
)

t0_zaga <- Sys.time()

fit.gb <- gamboostLSS(
  formula  = form_list,
  data     = dtr,
  families = as.families("ZAGA"),
  control  = boost_control(mstop = 2000, nu = 0.01)
)

t1_zaga <- Sys.time()
time_zaga <- as.numeric(difftime(t1_zaga, t0_zaga, units = "secs"))

saveRDS(fit.gb, file = file.path(outdir, "fit_zaga_gamboostlss_groton_windham.rds"))

# -----------------------------
# Save summaries
# -----------------------------
sink(file.path(outdir, "fit_zaga_summary.txt")); print(summary(fit.gb)); sink()
sink(file.path(outdir, "fit_zaga_summary_mu.txt")); print(summary(fit.gb$mu)); sink()
sink(file.path(outdir, "fit_zaga_summary_sigma.txt")); print(summary(fit.gb$sigma)); sink()
sink(file.path(outdir, "fit_zaga_summary_nu.txt")); print(summary(fit.gb$nu)); sink()

# -----------------------------
# Predictions
# -----------------------------
mu_tr <- as.numeric(predict(fit.gb, newdata = dtr, parameter = "mu", type = "response"))
sigma_tr <- as.numeric(predict(fit.gb, newdata = dtr, parameter = "sigma", type = "response"))
nu_tr <- as.numeric(predict(fit.gb, newdata = dtr, parameter = "nu", type = "response"))
pred_tr <- (1 - nu_tr) * mu_tr

lo80_tr  <- gamlss.dist::qZAGA(0.10, mu = mu_tr, sigma = sigma_tr, nu = nu_tr)
hi80_tr  <- gamlss.dist::qZAGA(0.90, mu = mu_tr, sigma = sigma_tr, nu = nu_tr)
lo95_tr  <- gamlss.dist::qZAGA(0.025, mu = mu_tr, sigma = sigma_tr, nu = nu_tr)
hi95_tr  <- gamlss.dist::qZAGA(0.975, mu = mu_tr, sigma = sigma_tr, nu = nu_tr)

mu_te <- as.numeric(predict(fit.gb, newdata = dte, parameter = "mu", type = "response"))
sigma_te <- as.numeric(predict(fit.gb, newdata = dte, parameter = "sigma", type = "response"))
nu_te <- as.numeric(predict(fit.gb, newdata = dte, parameter = "nu", type = "response"))
pred_te <- (1 - nu_te) * mu_te

lo80_te  <- gamlss.dist::qZAGA(0.10, mu = mu_te, sigma = sigma_te, nu = nu_te)
hi80_te  <- gamlss.dist::qZAGA(0.90, mu = mu_te, sigma = sigma_te, nu = nu_te)
lo95_te  <- gamlss.dist::qZAGA(0.025, mu = mu_te, sigma = sigma_te, nu = nu_te)
hi95_te  <- gamlss.dist::qZAGA(0.975, mu = mu_te, sigma = sigma_te, nu = nu_te)

# -----------------------------
# Metrics (quantile-wise peaks)
# -----------------------------
rmse_tr <- rmse(Y_tr, pred_tr)
rmse_te <- rmse(Y_te, pred_te)

qs <- c(0.5,0.6,0.7,0.75,0.8,0.85,0.90,0.95,0.99)

peak_list <- lapply(qs, function(q){
  peaks_tr <- rmse_peaks(Y_tr, pred_tr, q = q)
  peaks_te <- rmse_peaks(Y_te, pred_te, q = q)
  
  data.table(
    split = c("train","test"),
    quantile = q,
    peak_threshold = c(peaks_tr$threshold, peaks_te$threshold),
    n_peaks = c(peaks_tr$n, peaks_te$n),
    rmse_peaks = c(peaks_tr$rmse, peaks_te$rmse)
  )
})

peak_metrics <- rbindlist(peak_list)

metrics <- data.table(
  split = c("train","test"),
  rmse = c(rmse_tr, rmse_te)
)

metrics <- merge(metrics, peak_metrics, by = "split", allow.cartesian = TRUE)
fwrite(metrics, file.path(outdir, "metrics.csv"))

# Save series predictions
series_predictions <- rbindlist(list(
  data.table(
    DATE = DATE_tr, actual = Y_tr,
    mu = mu_tr, sigma = sigma_tr, nu = nu_tr,
    zaga = pred_tr,
    lo80 = lo80_tr, hi80 = hi80_tr,
    lo95 = lo95_tr, hi95 = hi95_tr,
    part = "train"
  ),
  data.table(
    DATE = DATE_te, actual = Y_te,
    mu = mu_te, sigma = sigma_te, nu = nu_te,
    zaga = pred_te,
    lo80 = lo80_te, hi80 = hi80_te,
    lo95 = lo95_te, hi95 = hi95_te,
    part = "test"
  )
))
fwrite(series_predictions, file.path(outdir, "series_predictions.csv"))

# -----------------------------
# Plot
# -----------------------------
peaks_te_095 <- rmse_peaks(Y_te, pred_te, q = 0.95)

plot_df <- rbindlist(list(
  data.table(DATE = DATE_tr, actual = Y_tr, fitted = pred_tr,
             lo80 = lo80_tr, hi80 = hi80_tr, lo95 = lo95_tr, hi95 = hi95_tr,
             type = "In-sample fit"),
  data.table(DATE = DATE_te, actual = Y_te, fitted = pred_te,
             lo80 = lo80_te, hi80 = hi80_te, lo95 = lo95_te, hi95 = hi95_te,
             type = "Out-of-sample pred")
))

p <- ggplot(plot_df, aes(x = DATE)) +
  geom_ribbon(data = plot_df[type == "Out-of-sample pred"],
              aes(ymin = lo95, ymax = hi95),
              fill = "darkred", alpha = 0.15) +
  geom_ribbon(data = plot_df[type == "Out-of-sample pred"],
              aes(ymin = lo80, ymax = hi80),
              fill = "maroon", alpha = 0.30) +
  geom_line(aes(y = actual, color = "Actual"), linewidth = 0.45) +
  geom_line(data = plot_df[type=="In-sample fit"],
            aes(y = fitted, color = "In-sample fit"),
            linewidth = 0.35, linetype = "dashed") +
  geom_line(data = plot_df[type=="Out-of-sample pred"],
            aes(y = fitted, color = "Out-of-sample pred"),
            linewidth = 0.35, linetype = "dotted") +
  geom_vline(xintercept = split_date, linetype = "dashed") +
  scale_color_manual(values = c("Actual"="black","In-sample fit"="blue","Out-of-sample pred"="red")) +
  labs(
    title = "gamboostLSS ZAGA — Groton with Windham exogenous rain (current + lags)",
    subtitle = sprintf("Test RMSE = %.4f | Peak-RMSE(95%%) = %.4f", rmse_te, peaks_te_095$rmse),
    x = "Date", y = "Hourly precipitation", color = ""
  ) +
  theme_minimal(base_size = 12)

pdf(file.path(outdir, "actual_vs_zaga_groton_windham.pdf"), width = 12, height = 6)
print(p)
dev.off()

# -----------------------------
# Runtime end
# -----------------------------
t1_total <- Sys.time()
time_total <- as.numeric(difftime(t1_total, t0_total, units = "secs"))

runtime <- data.table(
  component = c("ZAGA_fit", "TOTAL_pipeline"),
  time_seconds = c(time_zaga, time_total)
)
fwrite(runtime, file.path(outdir, "runtime_seconds.csv"))

cat("\nDONE.\n")
cat("Outputs in: ", normalizePath(outdir), "\n")
cat("Runtime (sec): ZAGA=", time_zaga, " TOTAL=", time_total, "\n")