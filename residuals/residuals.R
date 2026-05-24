library(data.table)
library(ggplot2)
library(gamlss.dist)   # for pGA

save_zaga_residual_plots <- function(csv_path,
                                     file_name,
                                     out_dir = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/residuals",
                                     seed = 2,
                                     width = 12,
                                     height = 10,
                                     dpi = 300) {
  
  # -----------------------------
  # Create output directory if needed
  # -----------------------------
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # -----------------------------
  # Read predictions
  # -----------------------------
  preds <- fread(csv_path)
  preds[, DATE := as.POSIXct(DATE, tz = "UTC")]
  
  # -----------------------------
  # Keep only training data
  # -----------------------------
  train_preds <- preds[part == "train"]
  
  if (nrow(train_preds) == 0) {
    stop("No training rows found: part == 'train'")
  }
  
  # -----------------------------
  # Check required columns
  # -----------------------------
  needed_cols <- c("actual", "mu", "sigma", "nu", "zaga")
  missing_cols <- setdiff(needed_cols, names(train_preds))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # -----------------------------
  # Randomized quantile residuals function
  # -----------------------------
  rqres_zaga <- function(y, mu, sigma, nu, seed = 2) {
    n <- length(y)
    r <- numeric(n)
    
    # Positive observations
    pos <- y > 0
    if (any(pos)) {
      Fy_gamma <- pGA(y[pos], mu = mu[pos], sigma = sigma[pos])
      Fy <- nu[pos] + (1 - nu[pos]) * Fy_gamma
      Fy <- pmin(pmax(Fy, 1e-12), 1 - 1e-12)
      r[pos] <- qnorm(Fy)
    }
    
    # Zero observations
    zer <- y == 0
    if (any(zer)) {
      set.seed(seed)
      u <- runif(sum(zer), min = 0, max = nu[zer])
      u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
      r[zer] <- qnorm(u)
    }
    
    r
  }
  seed=sample(1:100000,1)
  # -----------------------------
  # Compute residuals
  # -----------------------------
  rq <- rqres_zaga(
    y     = train_preds$actual,
    mu    = train_preds$mu,
    sigma = train_preds$sigma,
    nu    = train_preds$nu,
    seed  = seed
  )
  
  # -----------------------------
  # Save plot
  # -----------------------------
  out_file <- file.path(out_dir, paste0(file_name, ".png"))
  
  png(out_file, width = width, height = height, units = "in", res = dpi)
  par(mfrow = c(1, 1))
  # 1. Histogram
  #hist(rq,
  #     freq = FALSE,
  #     xlab = "Normalized Q Residual",
  #     main = "Histogram",
  #     xlim = c(-8, 8),cex.lab=1.5,cex.main=2)
  #curve(dnorm(x),
  #      col = "darkblue",
  #      lwd = 2,
  #      add = TRUE)
  #lines(density(rq), col = "red")
  #legend(x = 2, y = 0.4, bty = "n",
  #       legend = c("Normal Distribution", "Density of QRes"),
  #      fill = c("darkblue", "red"),cex = 0.8)
  
  # 2. QQ plot
  car::qqPlot(rq,distribution='norm', main = "Normal Q-Q Plot",ylab='Quantile Residuals')
  #qqnorm(rq, main = "Normal Q-Q Plot",cex.lab=1.5,cex.main=2)
  #qqline(rq, col = "red", lwd = 2)
  #grid(col='darkgray')
  
  # 3. Residual vs fitted
  #plot((train_preds$zaga), rq,
  #     pch = 19,
  #     main = "Residual vs Fitted Plot",
  #     xlab = "Predicted",
  #     ylab = "Normalized Q Residual",cex.lab=1.5,cex.main=2)
  #abline(h = 0, col = "red", lwd = 2)
  
  # 4. Residual vs index
  #plot(rq, pch = 19,
  #     main = "Residual vs Index",
  #     xlab = "Index",
  #     ylab = "Normalized Q Residual",cex.lab=1.5,cex.main=2)
  #abline(h = 0, col = "red", lwd = 2,)
  
  dev.off()
  
  message("Residual plot saved to: ", out_file)
  
  invisible(list(
    residuals = rq,
    train_data = train_preds,
    output_file = out_file
  ))
}

#Model2 Case1
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_12hours_test/out_zaga_gamboostlss_all_preds/series_predictions.csv",
  file_name = "zaga_residuals_m2_case1"
)

#Model2 Case2
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_120hours_test/out_zaga_gamboostlss_all_preds/series_predictions.csv",
  file_name = "zaga_residuals_m2_case2"
)


#Model2 Case3
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/out_zaga_gamboostlss_all_preds/series_predictions.csv",
  file_name = "zaga_residuals_m2_case3"
)

#Model1 Case1
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_12hours_test/out_zaga_gamboostlss_all_preds/series_predictions.csv",
  file_name = "zaga_residuals_m1_case1"
)

#Model1 Case2
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_120hours_test/out_zaga_gamboostlss_all_preds/series_predictions.csv",
  file_name = "zaga_residuals_m1_case2"
)

#Model1 Case3
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/out_zaga_gamboostlss_all_preds_rainonly/series_predictions.csv",
  file_name = "zaga_residuals_m1_case3"
)

#Model3 Case1  Danbury -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_12hours_test/out_zaga_gamboostlss_meriden_danbury/series_predictions.csv",
  file_name = "zaga_residuals_m3_danbury_case1"
)
#Model3 Case2  Danbury -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_120hours_test/out_zaga_gamboostlss_meriden_danbury/series_predictions.csv",
  file_name = "zaga_residuals_m3_danbury_case2"
)
#Model3 Case3  Danbury -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/out_zaga_gamboostlss_meriden_danbury/series_predictions.csv",
  file_name = "zaga_residuals_m3_danbury_case3"
)

#Model3 Case1  Bradley -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_12hours_test/out_zaga_gamboostlss_meriden_bradley/series_predictions.csv",
  file_name = "zaga_residuals_m3_bradley_case1"
)
#Model3 Case2  bradley -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_120hours_test/out_zaga_gamboostlss_meriden_bradley/series_predictions.csv",
  file_name = "zaga_residuals_m3_bradley_case2"
)
#Model3 Case3  bradley -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/out_zaga_gamboostlss_meriden_bradley/series_predictions.csv",
  file_name = "zaga_residuals_m3_bradley_case3"
)

#Model3 Case1  newhaven -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_12hours_test/out_zaga_gamboostlss_meriden_newhaven/series_predictions.csv",
  file_name = "zaga_residuals_m3_newhaven_case1"
)
#Model3 Case2  newhaven -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_120hours_test/out_zaga_gamboostlss_meriden_newhaven/series_predictions.csv",
  file_name = "zaga_residuals_m3_newhaven_case2"
)
#Model3 Case3  newhaven -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/out_zaga_gamboostlss_meriden_newhaven/series_predictions.csv",
  file_name = "zaga_residuals_m3_newhaven_case3"
)

#Model3 Case1  windham -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_12hours_test/out_zaga_gamboostlss_meriden_windham/series_predictions.csv",
  file_name = "zaga_residuals_m3_windham_case1"
)
#Model3 Case2  windham -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_120hours_test/out_zaga_gamboostlss_meriden_windham/series_predictions.csv",
  file_name = "zaga_residuals_m3_windham_case2"
)
#Model3 Case3  windham -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/out_zaga_gamboostlss_meriden_windham/series_predictions.csv",
  file_name = "zaga_residuals_m3_windham_case3"
)

#Model3 Case1  groton -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_12hours_test/out_zaga_gamboostlss_meriden_groton/series_predictions.csv",
  file_name = "zaga_residuals_m3_groton_case1"
)
#Model3 Case2  groton -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_120hours_test/out_zaga_gamboostlss_meriden_groton/series_predictions.csv",
  file_name = "zaga_residuals_m3_groton_case2"
)
#Model3 Case3  groton -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/out_zaga_gamboostlss_meriden_groton/series_predictions.csv",
  file_name = "zaga_residuals_m3_groton_case3"
)

#Model3 Case1  sikorsky -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_12hours_test/out_zaga_gamboostlss_meriden_sikorsky/series_predictions.csv",
  file_name = "zaga_residuals_m3_sikorsky_case1"
)
#Model3 Case2  sikorsky -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_120hours_test/out_zaga_gamboostlss_meriden_sikorsky/series_predictions.csv",
  file_name = "zaga_residuals_m3_sikorsky_case2"
)
#Model3 Case3  sikorsky -
save_zaga_residual_plots(
  csv_path = "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/out_zaga_gamboostlss_meriden_sikorsky/series_predictions.csv",
  file_name = "zaga_residuals_m3_sikorsky_case3"
)


'/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/out_zaga_gamboostlss_meriden_allstations/series_predictions.csv'


#Model4 Case1 -
save_zaga_residual_plots(
  csv_path = '/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_12hours_test/out_zaga_gamboostlss_meriden_allstations/series_predictions.csv',
  file_name = "zaga_residuals_m4_case1"
)
#Model4 Case2 -
save_zaga_residual_plots(
  csv_path = '/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/predictions_120hours_test/out_zaga_gamboostlss_meriden_allstations/series_predictions.csv',
  file_name = "zaga_residuals_m4_case2"
)
#Model4 Case3 -
save_zaga_residual_plots(
  csv_path = '/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/gam_lss paper/out_zaga_gamboostlss_meriden_allstations/series_predictions.csv',
  file_name = "zaga_residuals_m4_case3"
)
