# =====================================================
# 01_mdwkde.R
# Computes a normalized Magnitude–Depth Weighted KDE
# =====================================================

# ---- Libraries ----
library(sf)
library(spatstat.geom)
library(spatstat.explore)
library(raster)

# ---- 1. Load data ----
# Earthquake CSV must have: longitude, latitude, magnitude, depth
eq <- read.csv("C:/Users/archi/Desktop/mdwkde/Earthquake Data (Mw 5.0 above, 2000–2024) - cleaned.csv")

# Philippines Level-2 shapefile
phl <- st_read("C:/Users/archi/Desktop/mdwkde/gadm41_PHL_shp/gadm41_PHL_2.shp", quiet = TRUE)

# ---- 2. Convert earthquakes to sf and project ----
eq_sf <- st_as_sf(
  eq,
  coords = c("longitude", "latitude"),
  crs = 4326
)

eq_utm  <- st_transform(eq_sf, 32651)
phl_utm <- st_transform(phl, 32651)

# ---- 3. Create spatstat window (Philippines boundary) ----
phl_union <- st_union(phl_utm)
phl_win   <- spatstat.geom::as.owin(phl_union)

# ---- 4. Extract coordinates ----
coords <- st_coordinates(eq_utm)

# ---- 5. Keep only earthquakes inside PH boundary ----
inside <- spatstat.geom::inside.owin(
  x = coords[,1],
  y = coords[,2],
  w = phl_win
)

coords_in <- coords[inside, ]
mag_in    <- eq$magnitude[inside]
depth_in  <- eq$depth[inside]

# ---- 6. Create ppp object ----
ppp <- spatstat.geom::ppp(
  x = coords_in[,1],
  y = coords_in[,2],
  window = phl_win
)

# ---- 7. Magnitude–Depth weights ----
depth_safe <- pmax(depth_in, 5)        # km
weights    <- (mag_in^2) / depth_safe  # MD weight

# ---- 8. MDWKDE ----
sigma_km <- 50000  # 50 km bandwidth

mdwkde <- density(
  ppp,
  weights = weights,
  sigma = sigma_km
)


# ---- 9. Normalize to 0–1 ----
v <- mdwkde$v

mdwkde_norm <- (mdwkde - min(v, na.rm = TRUE)) /
  (max(v, na.rm = TRUE) - min(v, na.rm = TRUE))

# ---- 10. Convert to raster and save ----
md_raster <- raster(mdwkde_norm)

crs(md_raster) <- "+proj=utm +zone=51 +datum=WGS84 +units=m +no_defs"

writeRaster(
  md_raster,
  filename = "mdwkde.tif",
  overwrite = TRUE
)

# ---- 11. Quick sanity check plot ----
plot(
  mdwkde_norm,
  main = "Normalized Magnitude–Depth Weighted KDE",
  col = hcl.colors(100, "Spectral", rev = TRUE)
)

cat("MDWKDE raster written to mdwkde.tif\n")
cat("Number of points used:", ppp$n, "\n")
