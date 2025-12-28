load(file="Objects.RData")
load(file="PRISMelevation.RData")


library(maps)
library(fields)

pdf(file="~/Downloads/GammaBarHat.pdf")
map('state',c("wyoming","montana","idaho"),mar = c(4.1, 4.1, 12, 0.1))
quilt.plot(c(wy.stn.info$lon,id.stn.info$lon,mt.stn.info$lon),
           c(wy.stn.info$lat,id.stn.info$lat,mt.stn.info$lat),
           1 - c(unlist(wy.gamma),unlist(id.gamma),unlist(mt.gamma)),
           add=TRUE,nx=100,ny=100)
dev.off()

pdf(file="~/Downloads/ChiHat.pdf")
map('state',c("wyoming","montana","idaho"),mar = c(4.1, 4.1, 12, 0.1))
quilt.plot(c(wy.stn.info$lon,id.stn.info$lon,mt.stn.info$lon),
           c(wy.stn.info$lat,id.stn.info$lat,mt.stn.info$lat),
           c(unlist(wy.chi),unlist(id.chi),unlist(mt.chi)),
           add=TRUE,nx=100,ny=100)
dev.off()

pdf(file="~/Downloads/ChiHatMinusGammaBarHat.pdf")
map('state',c("wyoming","montana","idaho"),mar = c(4.1, 4.1, 12, 0.1))
quilt.plot(c(wy.stn.info$lon,id.stn.info$lon,mt.stn.info$lon),
           c(wy.stn.info$lat,id.stn.info$lat,mt.stn.info$lat),
           c(unlist(wy.chi),unlist(id.chi),unlist(mt.chi)) - 
             (1 - c(unlist(wy.gamma),unlist(id.gamma),unlist(mt.gamma))),
           add=TRUE,nx=100,ny=100)
dev.off()



inds.x <- which(PRISMelevation$x >= -117 & PRISMelevation$x <= -104,arr.ind=TRUE)
inds.y <- which(PRISMelevation$y >= 41 & PRISMelevation$y <= 49,arr.ind=TRUE)
lon_prism <- PRISMelevation$x[inds.x]
lat_prism <- PRISMelevation$y[inds.y]
elev_prism <- PRISMelevation$z[inds.x,inds.y]
get_utm_zone <- function(lon) {
  # Handle special cases if needed, but for the US, this is often sufficient
  as.integer(floor((lon + 180) / 6) + 1)
}

latlon_prism <- expand.grid(lat_prism,lon_prism)

get_utm_zone <- function(lon) {
  # Handle special cases if needed, but for the US, this is often sufficient
  as.integer(floor((lon + 180) / 6) + 1)
}

latlon_prism <- expand.grid(lat_prism,lon_prism)
#longitude of SNOTEL station locations
lon_stns <- c(wy.stn.info$lon,mt.stn.info$lon,id.stn.info$lon)

#latitude of SNOTEL station locations
lat_stns <- c(wy.stn.info$lat,mt.stn.info$lat,id.stn.info$lat)

library(sf)
# 1. Create dataframe
data <- data.frame(#first 261 are from SNOTEL stations
  location_id = 1:length(c(lat_stns,latlon_prism[,1])),
  latitude = c(lat_stns,latlon_prism[,1]),
  longitude = c(lon_stns,latlon_prism[,2])
)

# 2. Convert the data frame to an sf object
# The initial CRS is WGS 84 (EPSG: 4326)
sf_points <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)

# 3. Reproject to a suitable CRS for the continental US 
# USA Contiguous Albers Equal Area Conic (EPSG: 5070)
# This CRS has units in meters, which is ideal for distance calculations.
sf_projected <- st_transform(sf_points, crs = 5070)

# 4. Calculate the pairwise distance matrix
# st_distance() returns a matrix of class 'units' with meters as the unit.
distance_matrix_meters <- st_distance(sf_projected[1:261,])

# 5. Convert the distance matrix from meters to kilometers
# The 'units' package handles this conversion automatically.
distance_matrix_km <- distance_matrix_meters / 1000

# Convert the sf object to a matrix of coordinates
xy_matrix <- st_coordinates(sf_projected)/1000
xy_matrix_stns <- xy_matrix[1:261,]
xy_matrix_prism <- xy_matrix[-(1:261),]

distance_matrix_km_numeric <- as.matrix(distance_matrix_km)

xvar <- xy_matrix[1:261,1]
yvar <- xy_matrix[1:261,2]

zvar <- 1 - c(unlist(wy.gamma),unlist(mt.gamma),unlist(id.gamma))
zvar2 <- c(unlist(wy.chi),unlist(mt.chi),unlist(id.chi))


library(sp)
library(spdep)
library(geoR)
d1 <- 0
d2 <- 300
# Conduct the Moran's I test for gamma bar
spatial_data <- SpatialPointsDataFrame(
  coords = xy_matrix_stns,
  data = data.frame(zvar = zvar)
)
neighbors_list <- dnearneigh(spatial_data, d1 = d1, d2 = d2)
distances_list <- nbdists(neighbors_list, spatial_data@coords)
inverse_distances_list <- lapply(distances_list, function(x) 1/x)
weights_matrix <- nb2listw(
  neighbors_list,
  glist = inverse_distances_list,
  style = "W"
)
moran_result <- moran.test(spatial_data$zvar, weights_matrix)
print(moran_result)

# Conduct the Moran's I test for chi
spatial_data <- SpatialPointsDataFrame(
  coords = xy_matrix_stns,
  data = data.frame(zvar = zvar2)
)
neighbors_list <- dnearneigh(spatial_data, d1 = d1, d2 = d2)
distances_list <- nbdists(neighbors_list, spatial_data@coords)
inverse_distances_list <- lapply(distances_list, function(x) 1/x)
weights_matrix <- nb2listw(
  neighbors_list,
  glist = inverse_distances_list,
  style = "W"
)
moran_result <- moran.test(spatial_data$zvar, weights_matrix)
print(moran_result)






create_custom_pairs_plot <- function(data) {
  # --- Input Validation ---
  # Ensure the input is a data frame or matrix
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input must be a data frame or matrix.")
  }
  
  # Ensure there are at least 2 variables for a meaningful pairs plot
  if (ncol(data) < 2) {
    stop("Input data must have at least 2 columns (variables).")
  }
  
  # --- Custom Panel Functions ---
  
  # 1. Panel for the diagonal: Histograms with density plots
  panel.hist <- function(x, ...) {
    # Set the user coordinate system to match the data range
    usr <- par("usr")
    on.exit(par(usr = usr))
    
    # Calculate histogram and density
    h <- hist(x, breaks = "FD", plot = FALSE)
    d <- density(x, na.rm = TRUE)
    
    # Dynamically set the y-axis limit to accommodate the density plot and label
    max_y <- max(h$density, d$y) * 1.25
    par(usr = c(usr[1:2], 0, max_y))
    
    # Plot histogram bars and add a density curve, removing both axes
    hist(x, breaks = "FD", col = "lightblue", border = "black", prob = TRUE, add = TRUE, xaxt = "n", yaxt = "n")
    lines(d, col = "red", lwd = 2)
    
    # Add the variable name at the top of the panel
    var_name <- colnames(data)[par("mfg")[2]]
    text(x = mean(usr[1:2]), y = max(h$density, d$y) * 1.05, labels = var_name, font = 2, cex = 1.2, pos = 3)
  }
  
  # 2. Panel for the upper triangle: Correlation heatplot with text
  panel.cor <- function(x, y, ...) {
    # Get the user coordinate system for the current panel
    usr <- par("usr")
    
    # Calculate the correlation coefficient
    r <- cor(x, y, use = "complete.obs")
    
    # Create a color palette from blue to red, centered on white at 0
    pal <- colorRampPalette(c("dodgerblue4", "white", "firebrick4"))(100)
    
    # Determine the color based on the correlation value
    # Scale r from [-1, 1] to [1, 100] for the palette index
    color_index <- round((r + 1) * 50)
    color_index <- max(1, min(100, color_index)) # Clamp to valid range
    rect(usr[1], usr[3], usr[2], usr[4], col = pal[color_index])
    
    # Add the correlation value as text in the center
    text(mean(usr[1:2]), mean(usr[3:4]), labels = format(r, digits = 2), col = "black", cex = 1.5)
  }
  
  # 3. Panel for the lower triangle: Scatterplots with LOESS smooth
  panel.scatter <- function(x, y, ...) {
    # Add points, removing both axes
    points(x, y, pch = 20, xaxt = "n", yaxt = "n", ...)
    # Add a LOESS smooth curve to the scatterplot
    lines(lowess(x, y), col = "red", lwd = 2)
  }
  
  # --- Main Plotting Call ---
  
  # Set up the plotting window and layout. We use smaller outer margins
  # since we're no longer placing labels there.
  par(mar = c(1, 1, 1, 1), oma = c(2, 2, 2, 2))
  
  # Call the pairs function with our custom panels
  pairs(data,
        lower.panel = panel.scatter,
        upper.panel = panel.cor,
        diag.panel = panel.hist,
        labels = NULL, # Remove the labels from the diagonal
        main = "", # Ensure no main title is drawn
        cex = 0.8,
        gap = 0.5 # Gap between panels
  )
}

dat_mat <- cbind(zvar,zvar2,xy_matrix[1:261,1],xy_matrix[1:261,2],elevs,avgprecip,avgtemp)
colnames(dat_mat) <- c("Gamma Bar Hat","Chi Hat","Easting","Northing","Elevation","Total Precip","Avg Temp")

pdf(file="~/Downloads/ScatterplotCorMatrix.pdf",w=1.3*8.5,h=1.2*9)
create_custom_pairs_plot(dat_mat)
dev.off()





pdf(file="SpatialPlot4c.pdf")
par(mfrow=c(1,1),mar = c(5.1, 4.1, 4.1, 2.1))
map('state',c("wyoming","montana","idaho"),mar = c(4.1, 4.1, 12, 0.1))
quilt.plot(lon_stns,lat_stns,elevs,xlab="Longitude",ylab="Latitude",nx=128,ny=128,add=TRUE)
map('state',c("wyoming","montana","idaho"),add = TRUE)
dev.off()

pdf(file="SpatialPlot4d.pdf")
par(mfrow=c(1,1),mar = c(5.1, 4.1, 4.1, 2.1))
map('state',c("wyoming","montana","idaho"),mar = c(4.1, 4.1, 12, 0.1))
quilt.plot(lon_stns,lat_stns,avgprecip,xlab="Longitude",ylab="Latitude",nx=128,ny=128,add=TRUE)
map('state',c("wyoming","montana","idaho"),add = TRUE)
dev.off()

pdf(file="SpatialPlot4h.pdf")
par(mfrow=c(1,1),mar = c(5.1, 4.1, 4.1, 2.1))
map('state',c("wyoming","montana","idaho"),mar = c(4.1, 4.1, 12, 0.1))
quilt.plot(lon_stns,lat_stns,avgtemp,xlab="Longitude",ylab="Latitude",nx=128,ny=128,add=TRUE)
map('state',c("wyoming","montana","idaho"),add = TRUE)
dev.off()




pdf(file="~/Downloads/Variograms.pdf",w=8.5,h=4.5)
par(mfrow=c(1,2))
geodata <- as.geodata(data.frame(x=xvar[1:261], y=yvar[1:261], z=zvar, elev=elevs, avgprecip=avgprecip,avgtemp=avgtemp),
                      coords.col=1:2, data.col=3, covar.col=4:6)
vario <- variog(geodata, trend="cte", max.dist=350) # Using "1st" for a linear trend
ml_fit_7 <- likfit(geodata,
                   trend = "cte", # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
plot(vario)
lines(ml_fit_7)

geodata <- as.geodata(data.frame(x=xvar[1:261], y=yvar[1:261], z=zvar2, elev=elevs, avgprecip=avgprecip,avgtemp=avgtemp),
                      coords.col=1:2, data.col=3, covar.col=4:6)
vario <- variog(geodata, trend="cte", max.dist=350) # Using "1st" for a linear trend
ml_fit_7 <- likfit(geodata,
                   trend = "cte", # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
plot(vario)
lines(ml_fit_7)
dev.off()




library(geoR)
geodata <- as.geodata(data.frame(x=xvar[1:261], y=yvar[1:261], z=zvar, elev=scale(elevs), avgprecip=scale(avgprecip), avgtemp=scale(avgtemp)),
                      coords.col=1:2, data.col=3, covar.col=4:6)

ml_fit_0 <- likfit(geodata,
                   trend = "cte", # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_0)

ml_fit_1 <- likfit(geodata,
                   trend = ~ elev, # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_1)

ml_fit_2 <- likfit(geodata,
                   trend = ~ avgtemp, # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_2)

ml_fit_3 <- likfit(geodata,
                   trend = ~ avgprecip, # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_3)

ml_fit_4 <- likfit(geodata,
                   trend = ~ avgprecip + avgtemp , # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_4)

ml_fit_5 <- likfit(geodata,
                   trend = ~ elev + avgtemp , # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_5)

ml_fit_6 <- likfit(geodata,
                   trend = ~ elev + avgprecip, # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_6)

ml_fit_7 <- likfit(geodata,
                   trend = ~ elev + avgprecip + avgtemp , # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_7)



geodata <- as.geodata(data.frame(x=xvar[1:261], y=yvar[1:261], z=zvar2, elev=scale(elevs), avgprecip=scale(avgprecip), avgtemp=scale(avgtemp)),
                      coords.col=1:2, data.col=3, covar.col=4:6)

ml_fit_0 <- likfit(geodata,
                   trend = "cte", # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_0)

ml_fit_1 <- likfit(geodata,
                   trend = ~ elev, # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_1)

ml_fit_2 <- likfit(geodata,
                   trend = ~ avgtemp, # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_2)

ml_fit_3 <- likfit(geodata,
                   trend = ~ avgprecip, # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_3)

ml_fit_4 <- likfit(geodata,
                   trend = ~ avgprecip + avgtemp , # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_4)

ml_fit_5 <- likfit(geodata,
                   trend = ~ elev + avgtemp , # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_5)

ml_fit_6 <- likfit(geodata,
                   trend = ~ elev + avgprecip, # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_6)

ml_fit_7 <- likfit(geodata,
                   trend = ~ elev + avgprecip + avgtemp , # or trend=~elev_stns for your specific case
                   cov.model = "exponential",
                   ini.cov.pars = c(.15, 100)) # Initial guesses for (partial sill, range)
AIC(ml_fit_7)
