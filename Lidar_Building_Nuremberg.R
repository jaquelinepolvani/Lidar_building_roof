# =============================================================================
# Building Extraction and Roof Geometry from Urban ALS
# Author: Jaqueline Lopes Polvani
# =============================================================================

# -----------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------

required_pkgs <- c("lidR", "terra", "sf", "dplyr", "purrr", "tidyr", "ggplot2",
                   "ggnewscale", "maptiles", "tidyterra", "patchwork")

missing_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(lidR)
library(terra)
library(sf)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(maptiles)
library(tidyterra)
library(patchwork)

# ---- Paths ------------------------------------------

las_path  <- "C:/Users/Jaque/Documents/Lidar class/Building/649_5480.laz"
dtm_path  <- "C:/Users/Jaque/Documents/Lidar class/Building/649_5480_DGM1.tif"

out_dir <- file.path(dirname(las_path), "outputs")
fig_dir <- file.path(dirname(las_path), "figures")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Write figures through temporary files to avoid file-locking issues.

open_png <- function(final_path, width = 900, height = 900) {
  while (!is.null(dev.list())) dev.off()
  tmp_path <- file.path(tempdir(), basename(final_path))
  assign(".png_tmp_path", tmp_path, envir = .GlobalEnv)
  assign(".png_final_path", final_path, envir = .GlobalEnv)
  png(tmp_path, width = width, height = height)
}
close_png <- function() {
  dev.off()
  file.copy(.png_tmp_path, .png_final_path, overwrite = TRUE)
}
safe_ggsave <- function(final_path, plot, ...) {
  tmp_path <- file.path(tempdir(), basename(final_path))
  ggsave(tmp_path, plot, ...)
  file.copy(tmp_path, final_path, overwrite = TRUE)
}

target_crs <- 25832

# Analysis parameters

min_height          <- 2.5
chm_res             <- 0.5
min_building_area   <- 20
footprint_p75_min   <- 3.0

roof_min_height     <- min_height
plane_dist_thresh   <- 0.15
min_plane_pts       <- 25
ransac_iterations   <- 300
max_facets_per_bldg <- 12
max_wall_slope_deg  <- 75
min_facet_area_m2   <- 5
min_facet_pts       <- 25

merge_slope_tol_deg  <- 8
merge_aspect_tol_deg <- 20
merge_low_slope_deg  <- 5
merge_max_gap_m      <- 2.0

gable_slope_tol_deg      <- 15
gable_aspect_tol_deg     <- 30
hip_slope_spread_tol_deg <- 20
facet_adjacency_tol_m    <- 1.5


# -----------------------------------------------------------------------
# 1. Load input data
# -----------------------------------------------------------------------

las <- readLAS(las_path)
if (is.empty(las)) stop("LAS file could not be read or is empty: ", las_path)
if (is.na(st_crs(las))) {
  st_crs(las) <- target_crs
}

dtm <- rast(dtm_path)
if (is.na(crs(dtm)) || crs(dtm) == "") {
  crs(dtm) <- paste0("EPSG:", target_crs)
}


# -----------------------------------------------------------------------
# 2. Normalize point heights
# -----------------------------------------------------------------------

las_norm <- normalize_height(las, dtm)
las_norm <- filter_poi(las_norm, Z >= -1, Z <= 100)


# -----------------------------------------------------------------------
# 3. Select building-classified points
# -----------------------------------------------------------------------

building_pts <- filter_poi(las_norm, Classification == 6, Z >= min_height)
if (npoints(building_pts) == 0) {
  stop("No points with Classification == 6 above min_height. ",
       "Check table(las_norm$Classification).")
}

chm <- rasterize_canopy(las_norm, res = chm_res, algorithm = dsmtin())

open_png(file.path(fig_dir, "00_raw_surface_model.png"))
plot(chm, main = "Raw ALS surface height model (unfiltered)",
     col = hcl.colors(50, "Inferno"), colNA = "steelblue3",
     plg = list(title = "Height (m)", title.cex = 0.9))
close_png()

# Subcircle radius preserves narrow building edges during rasterization.
building_mask <- rasterize_canopy(building_pts, res = chm_res,
                                  algorithm = p2r(subcircle = 0.4))
building_mask <- !is.na(building_mask)
names(building_mask) <- "building"

open_png(file.path(fig_dir, "01_building_points_mask.png"))
plot(building_mask, main = "Building-classified points, rasterized (before closing)",
     col = c("grey90", "goldenrod2"), legend = FALSE)
close_png()


# -----------------------------------------------------------------------
# 4. Extract building footprints
# -----------------------------------------------------------------------

building_mask_binary <- classify(building_mask,
                                 rcl = matrix(c(0, NA, 1, 1), ncol = 2, byrow = TRUE))

close_w <- matrix(1, 3, 3)
building_mask_binary <- focal(building_mask_binary, w = close_w, fun = "max", na.rm = TRUE)
building_mask_binary <- focal(building_mask_binary, w = close_w, fun = "min", na.rm = TRUE)

building_patches <- patches(building_mask_binary, directions = 8, allowGaps = FALSE)

buildings_poly <- as.polygons(building_patches, dissolve = TRUE) |> st_as_sf()
names(buildings_poly)[1] <- "building_id"
buildings_poly$building_id <- as.integer(buildings_poly$building_id)

buildings_poly <- buildings_poly |>
  st_make_valid() |>
  st_collection_extract("POLYGON") |>
  st_cast("MULTIPOLYGON") |>
  mutate(area_m2 = as.numeric(st_area(geometry))) |>
  filter(area_m2 >= min_building_area, !st_is_empty(geometry))

# Remove low-height objects that may be included in the building-classified mask.
z_p75 <- terra::extract(chm, terra::vect(buildings_poly),
                        fun = function(v) quantile(v, 0.75, na.rm = TRUE))
buildings_poly$z_p75 <- z_p75[, 2]
buildings_poly <- buildings_poly |> filter(!is.na(z_p75), z_p75 >= footprint_p75_min)

buildings_poly <- st_simplify(buildings_poly,
                              dTolerance = 0.3,
                              preserveTopology = TRUE) |>
  st_make_valid() |>
  st_collection_extract("POLYGON") |>
  st_cast("MULTIPOLYGON") |>
  filter(!st_is_empty(geometry))

# Recalculate area after geometry simplification.
buildings_poly <- buildings_poly |>
  mutate(area_m2 = as.numeric(st_area(geometry))) |>
  filter(area_m2 >= min_building_area, !st_is_empty(geometry)) |>
  mutate(building_id = row_number())


open_png(file.path(fig_dir, "02_building_footprints.png"))
plot(chm, main = "Building footprints over surface height model",
     col = hcl.colors(50, "Inferno"), colNA = "steelblue3",
     plg = list(title = "Height (m)", title.cex = 0.9))
plot(st_geometry(buildings_poly), border = "red", lwd = 1.5, add = TRUE)
legend("bottomright", legend = "Detected building footprint",
       col = "red", lwd = 1.5, bty = "o", bg = "white", cex = 0.8)
close_png()


# -----------------------------------------------------------------------
# 5. Assign points to building footprints
# -----------------------------------------------------------------------

las_high <- merge_spatial(building_pts, buildings_poly["building_id"], attribute = "building_id")
las_bldg <- filter_poi(las_high, !is.na(building_id))


# -----------------------------------------------------------------------
# 6. Segment and classify roof facets
# -----------------------------------------------------------------------

# RANSAC identifies the largest planar subset within a roof point cloud.
ransac_plane_inliers <- function(xyz, dist_thresh, iterations = 300) {
  xyz <- as.matrix(xyz)
  n <- nrow(xyz)
  if (n < 3) return(NULL)
  best_inliers <- integer(0)
  for (i in seq_len(iterations)) {
    idx <- sample.int(n, 3)
    p1 <- xyz[idx[1], ]; p2 <- xyz[idx[2], ]; p3 <- xyz[idx[3], ]
    v1 <- p2 - p1; v2 <- p3 - p1
    normal <- c(v1[2]*v2[3] - v1[3]*v2[2],
                v1[3]*v2[1] - v1[1]*v2[3],
                v1[1]*v2[2] - v1[2]*v2[1])
    norm_len <- sqrt(sum(normal^2))
    if (norm_len < 1e-8) next
    normal <- normal / norm_len
    d <- -sum(normal * p1)
    dist <- abs(xyz %*% normal + d)
    inliers <- which(dist <= dist_thresh)
    if (length(inliers) > length(best_inliers)) best_inliers <- inliers
  }
  if (length(best_inliers) < 3) return(NULL)
  best_inliers
}

segment_roof_facets <- function(xyz, dist_thresh, min_pts, max_facets = 12,
                                iterations = 300) {
  n <- nrow(xyz)
  facet_id <- rep(NA_integer_, n)
  remaining_idx <- seq_len(n)
  current_id <- 0L
  for (f in seq_len(max_facets)) {
    if (length(remaining_idx) < min_pts) break
    inliers_rel <- ransac_plane_inliers(xyz[remaining_idx, c("X","Y","Z")],
                                        dist_thresh, iterations = iterations)
    if (is.null(inliers_rel) || length(inliers_rel) < min_pts) break
    current_id <- current_id + 1L
    facet_id[remaining_idx[inliers_rel]] <- current_id
    remaining_idx <- remaining_idx[-inliers_rel]
  }
  facet_id
}

# PCA estimates the plane normal; slope and aspect are derived from that normal.
fit_plane_geometry <- function(xyz) {
  if (nrow(xyz) < 3) return(NULL)
  centered <- scale(as.matrix(xyz[, c("X","Y","Z")]), center = TRUE, scale = FALSE)
  pca <- prcomp(centered)
  normal <- pca$rotation[, 3]
  if (normal[3] < 0) normal <- -normal
  slope_deg  <- acos(normal[3]) * 180 / pi
  aspect_rad <- atan2(normal[1], normal[2])
  aspect_deg <- (aspect_rad * 180 / pi + 360) %% 360
  hull <- chull(xyz$X, xyz$Y)
  hull_pts <- xyz[c(hull, hull[1]), c("X","Y")]
  footprint_area <- abs(sum(hull_pts$X[-nrow(hull_pts)] * hull_pts$Y[-1] -
                              hull_pts$X[-1] * hull_pts$Y[-nrow(hull_pts)])) / 2
  true_area <- footprint_area / cos(slope_deg * pi / 180)
  list(slope_deg  = round(slope_deg, 1),
       aspect_deg = round(aspect_deg, 1),
       n_points   = nrow(xyz),
       area_m2    = round(true_area, 1),
       z_mean     = round(mean(xyz$Z), 2))
}

circular_angle_diff <- function(a, b) {
  d <- abs(a - b) %% 360
  pmin(d, 360 - d)
}

min_gap_between_pointsets <- function(pi, pj, n_sub = 400) {
  pi <- as.matrix(pi); pj <- as.matrix(pj)
  if (nrow(pi) == 0 || nrow(pj) == 0) return(Inf)
  if (nrow(pi) > n_sub) pi <- pi[sample(nrow(pi), n_sub), , drop = FALSE]
  if (nrow(pj) > n_sub) pj <- pj[sample(nrow(pj), n_sub), , drop = FALSE]
  d <- sqrt(outer(pi[,1], pj[,1], "-")^2 + outer(pi[,2], pj[,2], "-")^2)
  min(d)
}

merge_similar_facets <- function(xyz, facet_id, slope_tol, aspect_tol,
                                 low_slope, max_gap = 2.0) {
  ids <- setdiff(unique(facet_id), NA)
  if (length(ids) < 2) return(facet_id)
  geoms <- lapply(ids, function(id) {
    fit_plane_geometry(xyz[!is.na(facet_id) & facet_id == id, ])
  })
  names(geoms) <- as.character(ids)
  for (i in seq_along(ids)) {
    for (j in seq_along(ids)) {
      if (j <= i) next
      gi <- geoms[[as.character(ids[i])]]
      gj <- geoms[[as.character(ids[j])]]
      if (is.null(gi) || is.null(gj)) next
      both_flat <- gi$slope_deg < low_slope && gj$slope_deg < low_slope
      slope_diff  <- abs(gi$slope_deg - gj$slope_deg)
      aspect_diff <- circular_angle_diff(gi$aspect_deg, gj$aspect_deg)
      geom_ok <- both_flat || (slope_diff <= slope_tol && aspect_diff <= aspect_tol)
      if (!geom_ok) next
      # Require spatial proximity to avoid merging separated parallel planes.
      pi <- xyz[!is.na(facet_id) & facet_id == ids[i], c("X","Y")]
      pj <- xyz[!is.na(facet_id) & facet_id == ids[j], c("X","Y")]
      gap <- min_gap_between_pointsets(pi, pj)
      if (gap > max_gap) next
      facet_id[facet_id == ids[j]] <- ids[i]
    }
  }
  facet_id
}

compute_facet_adjacency <- function(xyz, facet_id, plane_ids, tol = 1.5) {
  n <- length(plane_ids)
  adj <- matrix(FALSE, n, n, dimnames = list(as.character(plane_ids),
                                             as.character(plane_ids)))
  if (n < 2) return(adj)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (j <= i) next
      pi <- xyz[!is.na(facet_id) & facet_id == plane_ids[i], c("X","Y")]
      pj <- xyz[!is.na(facet_id) & facet_id == plane_ids[j], c("X","Y")]
      gap <- min_gap_between_pointsets(pi, pj)
      if (gap < tol) adj[i, j] <- adj[j, i] <- TRUE
    }
  }
  adj
}

# Roof types are assigned from the number, orientation, slope, and adjacency of facets.
classify_roof_type <- function(facets, adj = NULL) {
  n <- nrow(facets)
  if (n == 0) return("unclassified")
  
  if (n == 1) {
    return(if (facets$slope_deg[1] < 7) "flat" else "mono-pitch (shed)")
  }
  
  if (n == 2) {
    slope_diff  <- abs(facets$slope_deg[1] - facets$slope_deg[2])
    aspect_diff <- circular_angle_diff(facets$aspect_deg[1], facets$aspect_deg[2])
    slope_ok  <- slope_diff <= gable_slope_tol_deg
    aspect_ok <- abs(aspect_diff - 180) <= gable_aspect_tol_deg
    adj_ok    <- is.null(adj) || isTRUE(adj[1, 2])
    if (slope_ok && aspect_ok && adj_ok) return("gable (opposite slopes)")
    return("complex/uncertain")
  }
  
  if (n %in% 3:4) {
    slope_spread <- max(facets$slope_deg) - min(facets$slope_deg)
    adj_ok <- TRUE
    if (!is.null(adj)) adj_ok <- all(rowSums(adj) > 0)
    if (slope_spread <= hip_slope_spread_tol_deg && adj_ok) {
      return("hip (similar slopes)")
    }
    return("complex/uncertain")
  }
  
  "complex/uncertain"
}

building_ids <- sort(unique(las_bldg@data$building_id))
roof_results <- list()

# Detect roof planes independently for each segmented building.
for (bid in building_ids) {
  
  b_pts <- filter_poi(las_bldg, building_id == bid)
  if (npoints(b_pts) < min_plane_pts) next
  
  roof_pts <- filter_poi(b_pts, Z >= roof_min_height)
  if (npoints(roof_pts) < min_plane_pts) next
  
  roof_xyz <- roof_pts@data[, c("X", "Y", "Z")]
  
  facet_id <- segment_roof_facets(
    roof_xyz,
    dist_thresh = plane_dist_thresh,
    min_pts     = min_facet_pts,
    max_facets  = max_facets_per_bldg,
    iterations  = ransac_iterations
  )
  
  facet_id <- merge_similar_facets(
    roof_xyz, facet_id,
    slope_tol  = merge_slope_tol_deg,
    aspect_tol = merge_aspect_tol_deg,
    low_slope  = merge_low_slope_deg,
    max_gap    = merge_max_gap_m
  )
  
  plane_ids <- sort(unique(facet_id[!is.na(facet_id)]))
  if (length(plane_ids) == 0) next
  
  facets <- map_dfr(plane_ids, function(pid) {
    sub <- roof_xyz[!is.na(facet_id) & facet_id == pid, ]
    if (nrow(sub) < min_facet_pts) return(NULL)
    geo <- fit_plane_geometry(sub)
    if (is.null(geo)) return(NULL)
    as.data.frame(c(building_id = bid, plane_id = pid, geo))
  })
  
  facets <- facets |> filter(slope_deg <= max_wall_slope_deg,
                             area_m2 >= min_facet_area_m2)
  if (nrow(facets) == 0) next
  
  kept_ids <- facets$plane_id
  adj <- compute_facet_adjacency(roof_xyz, facet_id, kept_ids,
                                 tol = facet_adjacency_tol_m)
  
  roof_type <- classify_roof_type(facets, adj)
  facets$roof_type <- roof_type
  
  roof_results[[as.character(bid)]] <- facets
}

roof_planes <- bind_rows(roof_results)


if (nrow(roof_planes) == 0) {
  stop(
    "No roof facets were detected for any building.\n",
    "Try loosening plane_dist_thresh (e.g. 0.25) or lowering min_facet_pts ",
    "(e.g. 10-15), then re-run section 6."
  )
}


# -----------------------------------------------------------------------
# 7. Join roof attributes to building footprints
# -----------------------------------------------------------------------

roof_summary <- roof_planes |>
  group_by(building_id) |>
  summarise(
    n_facets        = n(),
    mean_slope_deg  = round(mean(slope_deg), 1),
    roof_area_m2    = round(sum(area_m2), 1),
    roof_type       = first(roof_type),
    .groups = "drop"
  )

buildings_out <- buildings_poly |> left_join(roof_summary, by = "building_id")


# -----------------------------------------------------------------------
# 8. Export
# -----------------------------------------------------------------------

tmp_gpkg <- file.path(tempdir(), "buildings.gpkg")
st_write(buildings_out, tmp_gpkg, delete_dsn = TRUE, quiet = TRUE)
file.copy(tmp_gpkg, file.path(out_dir, "buildings.gpkg"), overwrite = TRUE)
write.csv(roof_planes, file.path(out_dir, "roof_planes.csv"), row.names = FALSE)


# -----------------------------------------------------------------------
# 9. Generate overview figures
# -----------------------------------------------------------------------

chm_df <- as.data.frame(aggregate(chm, fact = 4, fun = "mean"), xy = TRUE, na.rm = FALSE)
names(chm_df)[3] <- "height"

buildings_plot <- buildings_out |>
  mutate(roof_type_plot = if_else(is.na(roof_type), "no facet detected", roof_type))

roof_type_colors <- c(
  "flat"                    = "#66C2A5",
  "mono-pitch (shed)"       = "#FC8D62",
  "gable (opposite slopes)" = "#8DA0CB",
  "hip (similar slopes)"    = "#A6D854",
  "complex/uncertain"       = "#E78AC3",
  "no facet detected"       = "grey60"
)

p_rooftypes_map <- ggplot() +
  geom_raster(data = chm_df, aes(x, y, fill = height)) +
  scale_fill_gradient(low = "grey15", high = "grey85", na.value = "steelblue3",
                      name = "Height (m)") +
  ggnewscale::new_scale_fill() +
  geom_sf(data = buildings_plot, aes(fill = roof_type_plot),
          colour = "black", linewidth = 0.15) +
  scale_fill_manual(values = roof_type_colors, name = "Roof type") +
  labs(title = "Buildings coloured by roof type",
       subtitle = "Grey = mask detected an object here, but no valid roof facet was found") +
  theme_minimal() +
  theme(axis.title = element_blank())

safe_ggsave(file.path(fig_dir, "03_roof_types.png"), p_rooftypes_map,
            width = 9, height = 9, dpi = 150)

p_slope <- ggplot(roof_planes, aes(x = slope_deg)) +
  geom_histogram(binwidth = 3, fill = "steelblue", colour = "white") +
  labs(title = "Distribution of detected roof facet slopes",
       x = "Slope (degrees from horizontal)", y = "Number of facets") +
  theme_minimal()
safe_ggsave(file.path(fig_dir, "04_slope_histogram.png"), p_slope, width = 7, height = 5)


# -----------------------------------------------------------------------
# 10. Generate building-scale examples
# -----------------------------------------------------------------------

plot_building_example <- function(bid, label, title_text, buffer_m = 40) {
  b_geom     <- buildings_out[buildings_out$building_id == bid, ]
  b_centroid <- st_centroid(b_geom)
  
  b_bbox <- st_bbox(b_geom)
  b_diag <- sqrt((b_bbox["xmax"] - b_bbox["xmin"])^2 + (b_bbox["ymax"] - b_bbox["ymin"])^2)
  buffer_m <- max(buffer_m, b_diag * 0.9)
  
  zoom_bbox <- st_bbox(st_buffer(b_centroid, buffer_m))
  xlim <- c(zoom_bbox["xmin"], zoom_bbox["xmax"])
  ylim <- c(zoom_bbox["ymin"], zoom_bbox["ymax"])
  
  chm_zoom       <- crop(chm, ext(zoom_bbox))
  buildings_zoom <- st_crop(buildings_out, zoom_bbox)
  
  p_chm <- ggplot() +
    geom_spatraster(data = chm_zoom) +
    scale_fill_gradient(low = "grey15", high = "grey85", na.value = "steelblue3",
                        name = "Height (m)") +
    geom_sf(data = buildings_zoom, fill = NA, colour = "red", linewidth = 1) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = "Detected footprint (from ALS)") +
    theme_minimal() +
    theme(axis.title = element_blank(), axis.text = element_blank())
  
  p_sat <- tryCatch({
    zoom_sf  <- st_as_sfc(zoom_bbox, crs = st_crs(buildings_out))
    sat_tile <- maptiles::get_tiles(zoom_sf, provider = "Esri.WorldImagery",
                                    crop = TRUE, zoom = 19)
    ggplot() +
      tidyterra::geom_spatraster_rgb(data = sat_tile) +
      geom_sf(data = buildings_zoom, fill = NA, colour = "red", linewidth = 1) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      labs(title = "Satellite imagery (same extent)") +
      theme_minimal() +
      theme(axis.title = element_blank(), axis.text = element_blank())
  }, error = function(e) {
    message("Could not fetch satellite basemap: ", conditionMessage(e)); NULL
  })
  
  combined <- if (!is.null(p_sat)) (p_chm | p_sat) else p_chm
  combined <- combined + patchwork::plot_annotation(title = title_text)
  
  safe_ggsave(file.path(fig_dir, paste0("05_", label, "_building", bid, ".png")),
              combined, width = if (!is.null(p_sat)) 11 else 6, height = 6, dpi = 150)
  invisible(NULL)
}

representative_candidates <- roof_summary |>
  left_join(st_drop_geometry(buildings_out)[, c("building_id", "area_m2")],
            by = "building_id") |>
  filter(roof_type %in% c("gable (opposite slopes)", "hip (similar slopes)"),
         area_m2 >= 80, area_m2 <= 400)

if (nrow(representative_candidates) > 0) {
  representative_bid <- representative_candidates$building_id[
    which.max(representative_candidates$area_m2)]
  plot_building_example(
    representative_bid, "representative",
    title_text = "Representative example (80-400 sqm footprint, cleanly classified)"
  )
}

complex_bid <- roof_summary$building_id[which.max(roof_summary$n_facets)]
plot_building_example(complex_bid, "complex_case",
                      title_text = "Complex/atypical example (most facets detected)",
                      buffer_m = 60)

roof_type_gallery <- roof_summary |>
  left_join(st_drop_geometry(buildings_out)[, c("building_id", "area_m2")],
            by = "building_id") |>
  filter(roof_type != "unclassified") |>
  group_by(roof_type) |>
  slice_max(area_m2, n = 1, with_ties = FALSE) |>
  ungroup()


for (i in seq_len(nrow(roof_type_gallery))) {
  gallery_bid   <- roof_type_gallery$building_id[i]
  gallery_type  <- roof_type_gallery$roof_type[i]
  gallery_label <- paste0("gallery_", gsub("[^a-z0-9]+", "_", tolower(gallery_type)))
  plot_building_example(
    gallery_bid, gallery_label,
    title_text = paste0("Roof type: ", gallery_type, " (building ", gallery_bid, ")")
  )
}

# -----------------------------------------------------------------------
# 11. Zoom comparison: satellite vs roof-type classification
# -----------------------------------------------------------------------
# Two panels at the same extent: satellite only, then satellite with
# footprints coloured by roof type.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Centroid of the densest cluster of classified buildings.
pick_busy_area <- function(buildings, buffer_m = 150) {
  classified <- buildings[!is.na(buildings$roof_type), ]
  if (nrow(classified) < 2) return(NULL)
  coords <- st_coordinates(st_centroid(classified))
  counts <- vapply(seq_len(nrow(coords)), function(i) {
    d <- sqrt((coords[, 1] - coords[i, 1])^2 + (coords[, 2] - coords[i, 2])^2)
    sum(d <= buffer_m)
  }, integer(1))
  best <- which.max(counts)
  list(center = as.numeric(coords[best, ]), n_buildings = counts[best])
}

plot_zoom_comparison <- function(center_xy, buffer_m = 150, label = "area",
                                 title_text = NULL) {
  
  center_sf <- st_sfc(st_point(center_xy), crs = st_crs(buildings_out))
  zoom_bbox <- st_bbox(st_buffer(center_sf, buffer_m))
  
  xlim <- c(zoom_bbox["xmin"], zoom_bbox["xmax"])
  ylim <- c(zoom_bbox["ymin"], zoom_bbox["ymax"])
  
  zoom_sf  <- st_as_sfc(zoom_bbox, crs = st_crs(buildings_out))
  sat_tile <- tryCatch(
    maptiles::get_tiles(zoom_sf, provider = "Esri.WorldImagery",
                        crop = TRUE, zoom = 19),
    error = function(e) {
      message("Satellite basemap unavailable: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(sat_tile)) return(invisible(NULL))
  
  buildings_zoom <- st_crop(buildings_out, zoom_bbox) |>
    mutate(roof_type_plot = if_else(is.na(roof_type),
                                    "no facet detected", roof_type))
  
  p_sat <- ggplot() +
    tidyterra::geom_spatraster_rgb(data = sat_tile) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = "Satellite imagery") +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text  = element_blank(),
          plot.title = element_text(size = 11))
  
  # Partial alpha keeps roof features visible under the fills.
  p_class <- ggplot() +
    tidyterra::geom_spatraster_rgb(data = sat_tile, alpha = 0.45) +
    geom_sf(data = buildings_zoom, aes(fill = roof_type_plot),
            colour = "black", linewidth = 0.35, alpha = 0.75) +
    scale_fill_manual(values = roof_type_colors, name = "Roof type") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = "Same extent with roof-type classification") +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text  = element_blank(),
          plot.title = element_text(size = 11),
          legend.position = "right")
  
  combined <- (p_sat | p_class) +
    patchwork::plot_annotation(
      title   = title_text %||% paste0("Zoom comparison: ", label),
      caption = paste0("Extent: \u00b1", buffer_m, " m around ",
                       round(center_xy[1]), ", ", round(center_xy[2]),
                       "  |  classified buildings in view: ",
                       sum(!is.na(buildings_zoom$roof_type)))
    )
  
  safe_ggsave(file.path(fig_dir, paste0("08_zoom_", label, ".png")),
              combined, width = 14, height = 7.5, dpi = 150)
  invisible(NULL)
}

# Default comparisons: tight residential cluster and wider neighbourhood.
busy <- pick_busy_area(buildings_out, buffer_m = 150)
if (!is.null(busy)) {
  plot_zoom_comparison(
    center_xy  = busy$center,
    buffer_m   = 150,
    label      = "busy_residential",
    title_text = paste0("Dense residential cluster (",
                        busy$n_buildings, " classified buildings within 150 m)")
  )
  plot_zoom_comparison(
    center_xy  = busy$center,
    buffer_m   = 350,
    label      = "neighbourhood_wide",
    title_text = "Wider neighbourhood view (same centre, larger extent)"
  )
}
