library(ggplot2)
library(maps)

# Station location data
station_locations <- data.frame(
  Station = c("USW00014707", "USW00014740", "USW00014758", "USW00054734", 
              "USW00054767", "USW00054788", "USW00094702"),
  Name = c("Groton", "Bradley", "New Haven", "Danbury", 
           "Windham", "Meriden", "Sikorsky"),
  lat = c(41.3275, 41.9375, 41.26389, 41.37139, 41.74194, 41.50972, 41.1583),
  long = c(-72.04944, -72.6819, -72.88722, -73.48278, -72.18361, -72.82778, -73.1289)
)

# Get US map data
us_map <- map_data("state")

# Extract Connecticut
ct_map <- subset(us_map, region == "connecticut")

# Plot
ggplot() +
  geom_polygon(data = ct_map,
               aes(x = long, y = lat, group = group),
               fill = "grey90",
               color = "black") +
  geom_point(data = station_locations,
             aes(x = long, y = lat),
             color = "red",
             size = 3) +
  geom_text(data = station_locations,
            aes(x = long, y = lat, label = Name),
            vjust = -1,
            size = 4) +
  coord_fixed(1.3) +
  labs(
    title = "CT Weather Stations",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()

