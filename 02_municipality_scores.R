# =====================================================
# 02_municipality_scores.R
# Simplify polygons + export web-ready GeoJSON
# =====================================================

library(sf)
library(raster)
library(exactextractr)
library(rmapshaper)

# ---- 1. Load MDWKDE raster ----
md <- raster("mdwkde.tif")

# ---- 2. Load PH Level-2 boundaries ----
phl <- st_read("gadm41_PHL_shp/gadm41_PHL_2.shp", quiet = TRUE)

# ---- 3. Match CRS ----
phl <- st_transform(phl, crs(md))

# ---- 4. Extract mean hazard per municipality ----
phl$hazard <- exact_extract(md, phl, "mean")

# ---- 5. Reproject for web ----
phl_web <- st_transform(phl, 4326)

# ---- 6. SIMPLIFY GEOMETRY (THIS IS THE KEY) ----
# keep = 0.05 is usually safe for municipality maps
phl_simple <- rmapshaper::ms_simplify(
  phl_web,
  keep = 0.05,
  keep_shapes = TRUE
)

# ---- 7. Write GeoJSON ----
st_write(
  phl_simple,
  "municipal_hazard.geojson",
  delete_dsn = TRUE
)

cat("Simplified municipal_hazard.geojson written\n")