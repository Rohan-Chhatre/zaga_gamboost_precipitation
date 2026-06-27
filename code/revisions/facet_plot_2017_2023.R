library(data.table)

path_M <- "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/concatenated_files/merged_imputed_final/USW00054788_merged_imputed_2004_2024.csv"

DT <- fread(path_M)

DT[, DATE := as.POSIXct(DATE, tz = "UTC")]
DT[, Year := format(DATE, "%Y")]

# Keep only 2017-2023
DT <- DT[
  DATE >= as.POSIXct("2017-01-01", tz = "UTC") &
    DATE <  as.POSIXct("2024-01-01", tz = "UTC")
]

yrs <- sort(unique(DT$Year))

# 7 panels
par(
  mfrow = c(4, 2),
  mar = c(3, 3, 3, 1),
  oma = c(2, 2, 3, 1)
)

for (yr in yrs) {
  
  x <- DT[Year == yr, HourlyPrecipitation]
  
  ts.plot(
    x,
    ylab = "",
    xlab = "",
    main = yr,
    lwd = 1,
    ylim=c(0,30)
  )
}

mtext(
  "Meriden Hourly Precipitation (2017–2023)",
  outer = TRUE,
  cex = 1.3
)

#--------------------------------------------------------------------------------------------
library(data.table)

path_M <- "/Users/rohan/Documents/Rohan/Precip_project/nalini_ct_airport_met/concatenated_files/merged_imputed_final/USW00054788_merged_imputed_2004_2024.csv"

DT <- fread(path_M)

DT[, DATE := as.POSIXct(DATE, tz = "UTC")]

DT <- DT[
  DATE >= as.POSIXct("2004-01-01", tz = "UTC") &
    DATE <  as.POSIXct("2024-01-01", tz = "UTC")
]

DT[, rain := as.numeric(HourlyPrecipitation)]
DT[, year := format(DATE, "%Y")]

# Rainfall categories
DT[, rain_cat := fifelse(
  rain == 0, "Dry: 0 mm",
  fifelse(rain <= 2, "Light: (0, 2] mm",
          fifelse(rain <= 10, "Moderate: (2, 10] mm",
                  "Heavy: >10 mm"))
)]

tab <- DT[
  ,
  .N,
  by = .(year, rain_cat)
]

tab[, prop := N / sum(N), by = year]

# Make wide matrix for barplot
wide <- dcast(tab, rain_cat ~ year, value.var = "prop", fill = 0)

cat_order <- c(
  "Dry: 0 mm",
  "Light: (0, 2] mm",
  "Moderate: (2, 10] mm",
  "Heavy: >10 mm"
)

wide <- wide[match(cat_order, rain_cat)]

mat <- as.matrix(wide[, -1])
rownames(mat) <- wide$rain_cat
