# =============================================================================
# 07-08. cluster-direction and compositional reassembly
# =============================================================================
# Run the numbered sections sequentially from the project root used in 01-06.
# Workflow: local spatial clusters of richness change -> cluster-direction
# cross-classification -> temporal compositional reassembly within each guild.
# No plotting, graphics devices, or colour settings are included. GeoTIFFs are
# analytical outputs, not rendered figures. There is no workflow wrapper or
# separate configuration file. All settings are defined in Section 0.

# ---- 0. Paths and settings ---------------------------------------------------
project_dir <- getwd()
guild_file <- file.path(
  project_dir, "outputs", "thermal_guild_classification", "species_thermal_guilds.csv"
)
distribution_dir <- file.path(
  project_dir, "outputs", "distributional_change_overlap"
)
richness_dir <- file.path(distribution_dir, "richness_maps")
change_dir <- file.path(distribution_dir, "richness_change_maps")
stack_dir <- file.path(distribution_dir, "species_stacks")
basin_file <- file.path(project_dir, "data", "raw", "river_basins.shp")

# These attributes are taken from the ORIGINAL basin vector, not basin_zones.tif.
# Preserve one output row per polygon feature; do not dissolve by these labels.
basin_code_field <- "river_code"
basin_group_field <- "river_bas0"
input_encoding <- "EUC-KR"  
output_encoding <- "EUC-KR"

# Exact filenames produced by script 05-06; placeholders are filled in order.
richness_file_template <- "richness_%s_%s.tif"  # guild, scenario
change_file_template <- "richness_change_%s.tif"  # guild; six named layers
stack_file_template <- "species_binary_%s.tif"  # scenario; fcode-named layers

output_dir <- file.path(project_dir, "outputs", "cluster_direction_reassembly")
lisa_dir <- file.path(output_dir, "spatial_clusters")
cluster_direction_dir <- file.path(output_dir, "cluster_direction")
reassembly_dir <- file.path(output_dir, "compositional_reassembly")
table_dir <- file.path(output_dir, "tables")
# Optional diagnostic CSVs retain cell-level Ii, Z.Ii, p-values and classification.
save_lisa_diagnostics <- TRUE
diagnostic_dir <- file.path(output_dir, "lisa_diagnostics")
overwrite_output <- FALSE  # Set TRUE deliberately before replacing old results.

# Original local-cluster settings. Keep the matrix-based distance calculation:
# geographic coordinates are treated as planar coordinates in degrees, NOT km.
# Changing knn_longlat to TRUE changes the analysis rather than only formatting.
n_neighbors <- 8L
knn_longlat <- FALSE
significance_level <- 0.05
p_adjust_method <- "BH"  # Applied separately within each guild-scenario raster.

# Original spatial-analysis scope: three guilds x three SSPs x two future periods.
current_scenario <- "current_normal"
future_scenarios <- data.frame(
  scenario = c(paste0("nearf_ssp", c(126, 370, 585)),
               paste0("farf_ssp", c(126, 370, 585))),
  period = rep(c("near-future", "far-future"), each = 3L),
  period_code = rep(c("nearf", "farf"), each = 3L),
  SSP = rep(c("SSP126", "SSP370", "SSP585"), 2L),
  stringsAsFactors = FALSE
)
guild_lookup <- data.frame(
  thermal.guild = c("cold-affinity", "intermediate-affinity", "warm-affinity"),
  guild = c("cold", "intermediate", "warm"), stringsAsFactors = FALSE
)
# Original reassembly scope: current versus far-future SSP5-8.5 only.
# This does NOT request beta diversity for all 18 spatial-analysis combinations.
reassembly_scenario <- "farf_ssp585"

# Install required packages once in your local R console, not in this script:
# install.packages(c("terra", "sf", "spdep"))
# Optional: install.packages("dbscan") enables spdep's planar k-d-tree search.
required_packages <- c("terra", "sf", "spdep")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install the required R package before continuing: ", package)
  }
}
library(terra)
library(sf)
library(spdep)

# ---- 1. Validate inputs and prepare exact file paths --------------------------
if (!is.numeric(n_neighbors) || length(n_neighbors) != 1L ||
    is.na(n_neighbors) || !is.finite(n_neighbors) || n_neighbors < 1L ||
    n_neighbors != as.integer(n_neighbors)) {
  stop("n_neighbors must be a positive integer.")
}
if (length(significance_level) != 1L || !is.finite(significance_level) ||
    significance_level <= 0 || significance_level >= 1 ||
    length(p_adjust_method) != 1L || is.na(p_adjust_method) ||
    !p_adjust_method %in% stats::p.adjust.methods) {
  stop("Check significance_level and p_adjust_method in Section 0.")
}
if (!is.logical(knn_longlat) || length(knn_longlat) != 1L || is.na(knn_longlat)) {
  stop("knn_longlat must be TRUE or FALSE.")
}
if (!reassembly_scenario %in% future_scenarios$scenario) {
  stop("reassembly_scenario must be one of the six configured future scenarios.")
}
for (input_file in c(guild_file, basin_file)) {
  if (!file.exists(input_file)) stop("Missing input file: ", input_file)
}

fish.info <- read.csv(
  guild_file, fileEncoding = input_encoding, stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!all(c("fcode", "thermal.guild") %in% names(fish.info)) ||
    anyDuplicated(names(fish.info)) || !nrow(fish.info)) {
  stop("Guild metadata require unique headers and columns fcode and thermal.guild.")
}
fish.info$fcode <- as.character(fish.info$fcode)
fish.info$thermal.guild <- as.character(fish.info$thermal.guild)
if (anyNA(fish.info[, c("fcode", "thermal.guild")]) ||
    any(!nzchar(trimws(fish.info$fcode))) || anyDuplicated(fish.info$fcode) ||
    !setequal(unique(fish.info$thermal.guild), guild_lookup$thermal.guild)) {
  stop("Species codes must be complete and unique; all three guilds are required.")
}

basins <- terra::vect(basin_file)
basin_fields <- c(basin_code_field, basin_group_field)
if (anyDuplicated(basin_fields) || !all(basin_fields %in% names(basins)) ||
    !nrow(basins) || !nzchar(terra::crs(basins)) ||
    !all(terra::geomtype(basins) == "polygons")) {
  stop("Basin input requires polygons, a defined CRS, and the two distinct ",
       "attributes selected in Section 0. Missing river_code is not replaced ",
       "automatically by BASIN_ID.")
}
basins <- basins[, basin_fields]
names(basins) <- c("river_code", "river_bas0")
basin_attributes <- data.frame(
  polygon_id = seq_len(nrow(basins)), as.data.frame(basins), check.names = FALSE
)

# Do not crop/dissolve this vector or impose additional region exclusions here.
# The analytical domain is the non-NA domain of the supplied change raster.
all_scenarios <- c(current_scenario, future_scenarios$scenario)
richness_files <- matrix(
  NA_character_, nrow = nrow(guild_lookup), ncol = length(all_scenarios),
  dimnames = list(guild_lookup$guild, all_scenarios)
)
for (guild in guild_lookup$guild) {
  richness_files[guild, ] <- file.path(
    richness_dir, sprintf(richness_file_template, guild, all_scenarios)
  )
}
change_files <- setNames(
  file.path(change_dir, sprintf(change_file_template, guild_lookup$guild)),
  guild_lookup$guild
)
stack_files <- setNames(
  file.path(stack_dir, sprintf(stack_file_template,
                               c(current_scenario, reassembly_scenario))),
  c(current_scenario, reassembly_scenario)
)
required_files <- unique(c(as.vector(richness_files), change_files, stack_files))
if (any(!file.exists(required_files))) {
  stop("Missing input rasters: ", paste(required_files[!file.exists(required_files)],
                                        collapse = ", "))
}

spatial_combinations <- expand.grid(
  scenario = future_scenarios$scenario, guild = guild_lookup$guild,
  stringsAsFactors = FALSE
)
spatial_combinations <- cbind(
  spatial_combinations,
  future_scenarios[match(spatial_combinations$scenario, future_scenarios$scenario),
                   c("period", "period_code", "SSP"), drop = FALSE]
)
rownames(spatial_combinations) <- NULL
spatial_combinations$lisa_file <- file.path(
  lisa_dir, paste0("lisa_", spatial_combinations$guild, "_",
                   spatial_combinations$scenario, ".tif")
)
spatial_combinations$cluster_direction_file <- file.path(
  cluster_direction_dir, paste0("cluster_direction_", spatial_combinations$guild,
                                "_", spatial_combinations$scenario, ".tif")
)
spatial_combinations$diagnostic_file <- if (save_lisa_diagnostics) {
  file.path(diagnostic_dir, paste0("lisa_diagnostics_", spatial_combinations$guild,
                                   "_", spatial_combinations$scenario, ".csv"))
} else {
  rep(NA_character_, nrow(spatial_combinations))
}

# Reject accidental replacement before creating any analytical output.
if (!overwrite_output && dir.exists(output_dir) &&
    length(list.files(output_dir, recursive = TRUE, all.files = TRUE,
                      no.. = TRUE)) > 0L) {
  stop("output_dir already contains files. Choose an empty directory or set ",
       "overwrite_output <- TRUE deliberately. No files have been deleted.")
}
for (directory in c(output_dir, lisa_dir, cluster_direction_dir, reassembly_dir,
                    table_dir, if (save_lisa_diagnostics) diagnostic_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

# ---- 2. Save class definitions and input metadata -----------------------------
lisa_lookup <- data.frame(
  value = 1:5, class = c("LL", "LH", "NS", "HL", "HH"),
  description = c("Low-low cluster", "Low-high spatial outlier", "Not significant",
                  "High-low spatial outlier", "High-high cluster"),
  stringsAsFactors = FALSE
)
class_lut <- data.frame(
  value = 1:7, class = c("H.P", "H.0", "H.N", "C.P", "C.0", "C.N", "other"),
  spatial_cluster = c(rep("hotspot", 3L), rep("coldspot", 3L), "other"),
  richness_change = c("positive", "no change", "negative", "positive",
                      "no change", "negative", "not classified"),
  stringsAsFactors = FALSE
)
# HH/LL describe relative spatial clusters of change, NOT the sign of change.
# Codes 1-6 combine HH/LL with the actual sign; LH/HL/NS receive code 7.
write.csv(lisa_lookup, file.path(table_dir, "lisa_class_lookup.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(class_lut, file.path(table_dir, "cluster_direction_class_lookup.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(spatial_combinations, file.path(table_dir, "spatial_output_manifest.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(fish.info[, c("fcode", "thermal.guild"), drop = FALSE],
          file.path(table_dir, "analysed_species.csv"), row.names = FALSE,
          fileEncoding = output_encoding)
write.csv(basin_attributes, file.path(table_dir, "basin_feature_lookup.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 3. Identify local spatial clusters of richness change --------------------
lisa_area_rows <- lisa_count_rows <- spatial_summary_rows <-
  vector("list", nrow(spatial_combinations))

for (i in seq_len(nrow(spatial_combinations))) {
  guild <- spatial_combinations$guild[i]
  scenario <- spatial_combinations$scenario[i]
  scenario_info <- spatial_combinations[i, c("guild", "scenario", "period", "SSP"),
                                        drop = FALSE]
  message("Local spatial clusters: ", guild, " / ", scenario)
  
  current_richness <- terra::rast(richness_files[guild, current_scenario])
  future_richness <- terra::rast(richness_files[guild, scenario])
  if (terra::nlyr(current_richness) != 1L || terra::nlyr(future_richness) != 1L ||
      !nzchar(terra::crs(current_richness))) {
    stop("Richness inputs must be single-layer rasters with a defined CRS.")
  }
  terra::compareGeom(current_richness, future_richness, stopOnError = TRUE)
  change <- future_richness - current_richness
  names(change) <- scenario
  
  # Confirm that the change stack from 05-06 describes exactly the same cells.
  change_stack <- terra::rast(change_files[[guild]])
  if (anyDuplicated(names(change_stack)) || !scenario %in% names(change_stack)) {
    stop("Missing or duplicated scenario layer in: ", change_files[[guild]])
  }
  stored_change <- change_stack[[scenario]]
  terra::compareGeom(change, stored_change, stopOnError = TRUE)
  vchange <- terra::values(change, mat = FALSE)
  vstored <- terra::values(stored_change, mat = FALSE)
  if (!identical(is.na(vchange), is.na(vstored)) ||
      !isTRUE(all.equal(as.numeric(vchange), as.numeric(vstored),
                        check.attributes = FALSE))) {
    stop("Stored richness change differs from future minus current: ", guild,
         " / ", scenario)
  }
  cell_id <- which(!is.na(vchange))
  if (length(cell_id) <= n_neighbors || any(!is.finite(vchange[cell_id]))) {
    stop("Need more valid, finite cells than n_neighbors: ", guild, " / ", scenario)
  }
  if (!is.finite(stats::sd(vchange[cell_id])) || stats::sd(vchange[cell_id]) == 0) {
    stop("Local Moran analysis is undefined for constant richness change: ",
         guild, " / ", scenario, ". No automatic NS classification is made.")
  }
  
  # Preserve source point order and pass a coordinate MATRIX, not an sf object, to knearneigh. 
  # No projection or river-network barrier is introduced.
  points <- terra::as.points(change, na.rm = TRUE)
  points_sf <- sf::st_as_sf(points)
  coords <- sf::st_coordinates(points_sf)
  x <- as.numeric(terra::values(points)[, 1])
  point_cells <- terra::cellFromXY(change, coords)
  if (!identical(as.numeric(point_cells), as.numeric(cell_id)) ||
      !isTRUE(all.equal(x, as.numeric(vchange[cell_id])))) {
    stop("Raster cell order and point order do not agree.")
  }
  if (knn_longlat && !isTRUE(sf::st_is_longlat(points_sf))) {
    stop("knn_longlat = TRUE requires geographic coordinates.")
  }
  neighbors <- spdep::knn2nb(
    spdep::knearneigh(coords, k = n_neighbors, longlat = knn_longlat,
                      use_kd_tree = TRUE),
    sym = FALSE
  )
  weights <- spdep::nb2listw(neighbors, style = "W", zero.policy = TRUE)
  local_moran <- spdep::localmoran(
    x, weights, zero.policy = TRUE, conditional = TRUE,
    alternative = "two.sided", mlvar = TRUE, adjust.x = FALSE
  )
  
  # Analytical Z-based p-values, not permutation p-values (source procedure).
  z_x <- as.numeric(scale(x))
  spatial_lag <- spdep::lag.listw(weights, x, zero.policy = TRUE)
  z_lag <- as.numeric(scale(spatial_lag))
  z_i <- local_moran[, "Z.Ii"]
  if (any(!is.finite(z_x)) || any(!is.finite(z_lag)) || any(!is.finite(z_i))) {
    stop("Undefined local statistic or standardized lag: ", guild, " / ", scenario,
         ". Inspect the data; undefined tests are not relabelled NS.")
  }
  p_raw <- 2 * stats::pnorm(-abs(z_i))
  p_adjusted <- stats::p.adjust(p_raw, method = p_adjust_method)
  significant <- p_adjusted < significance_level
  
  # Preserve source zero-boundary rules and assignment order exactly.
  quadrant <- rep("NS", length(x))
  quadrant[z_x >= 0 & z_lag >= 0 & significant] <- "HH"
  quadrant[z_x <= 0 & z_lag <= 0 & significant] <- "LL"
  quadrant[z_x >= 0 & z_lag <= 0 & significant] <- "HL"
  quadrant[z_x <= 0 & z_lag >= 0 & significant] <- "LH"
  points_sf$quad_code <- match(quadrant, lisa_lookup$class)
  lisa <- terra::rasterize(terra::vect(points_sf), change, field = "quad_code")
  names(lisa) <- "LISA_quad"
  terra::writeRaster(lisa, spatial_combinations$lisa_file[i],
                     overwrite = overwrite_output, datatype = "INT1U")
  
  if (save_lisa_diagnostics) {
    lisa_diagnostics <- data.frame(
      cell = cell_id, x = coords[, 1], y = coords[, 2], change = x,
      spatial_lag = spatial_lag, z_change = z_x, z_spatial_lag = z_lag,
      Ii = local_moran[, "Ii"], E.Ii = local_moran[, "E.Ii"],
      Var.Ii = local_moran[, "Var.Ii"], Z.Ii = z_i,
      p_raw = p_raw, p_adjusted = p_adjusted, significant = significant,
      LISA.spot = points_sf$quad_code, LISA.class = quadrant,
      guild = guild, scenario = scenario, row.names = NULL
    )
    write.csv(lisa_diagnostics, spatial_combinations$diagnostic_file[i],
              row.names = FALSE, fileEncoding = output_encoding)
    rm(lisa_diagnostics)
  }
  
  # Original LISA area method: default polygon extraction, WITHOUT exact weights.
  basins_local <- basins
  if (!terra::same.crs(basins_local, lisa)) {
    basins_local <- terra::project(basins_local, terra::crs(lisa))
  }
  cell_area <- terra::cellSize(lisa, unit = "km")
  area_layers <- NULL
  for (j in seq_len(nrow(lisa_lookup))) {
    class_area <- (lisa == lisa_lookup$value[j]) * cell_area
    names(class_area) <- paste0("area_", lisa_lookup$class[j], "_km2")
    area_layers <- if (is.null(area_layers)) class_area else c(area_layers, class_area)
  }
  area_by_basin <- terra::extract(
    area_layers, basins_local, fun = sum, na.rm = TRUE,
    exact = FALSE, touches = FALSE, small = TRUE
  )
  basin_order <- match(basin_attributes$polygon_id, area_by_basin$ID)
  if (anyNA(basin_order) || anyDuplicated(area_by_basin$ID)) {
    stop("Cannot match extracted LISA areas to basin features.")
  }
  lisa_area_rows[[i]] <- cbind(
    basin_attributes, area_by_basin[basin_order, -1, drop = FALSE],
    scenario_info[rep(1L, nrow(basin_attributes)), , drop = FALSE]
  )
  lisa_count_rows[[i]] <- data.frame(
    guild = guild, scenario = scenario, value = lisa_lookup$value,
    class = lisa_lookup$class, n_cells = tabulate(points_sf$quad_code, nbins = 5L)
  )
  spatial_summary_rows[[i]] <- cbind(scenario_info, data.frame(
    n_cells = length(x), n_zero_change = sum(x == 0),
    mean_change = mean(x), sd_change = stats::sd(x),
    n_significant = sum(significant), n_neighbors = n_neighbors,
    knn_longlat = knn_longlat, p_adjust_method = p_adjust_method,
    significance_level = significance_level
  ))
  rm(points, points_sf, neighbors, weights, local_moran, coords, change_stack,
     current_richness, future_richness, stored_change, area_layers, class_area)
}

lisa_area_all <- do.call(rbind, lisa_area_rows)
lisa_area_columns <- paste0("area_", lisa_lookup$class, "_km2")
# Retain the original complete-case filtering for this table only.
n_incomplete_lisa_area_rows <- sum(!complete.cases(lisa_area_all[, lisa_area_columns]))
lisa_area_all <- lisa_area_all[
  complete.cases(lisa_area_all[, lisa_area_columns]), , drop = FALSE
]
lisa_area_all$area_Other_km2 <- lisa_area_all$area_LH_km2 +
  lisa_area_all$area_NS_km2 + lisa_area_all$area_HL_km2
rownames(lisa_area_all) <- NULL
lisa_counts <- do.call(rbind, lisa_count_rows)
spatial_summary <- do.call(rbind, spatial_summary_rows)
write.csv(lisa_area_all, file.path(table_dir, "lisa_area_by_basin.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(lisa_counts, file.path(table_dir, "lisa_class_counts.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(spatial_summary, file.path(table_dir, "spatial_analysis_summary.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 4. Cross-classify spatial clusters and richness-change direction ---------
cluster_area_rows <- cluster_count_rows <- vector("list", nrow(spatial_combinations))
for (i in seq_len(nrow(spatial_combinations))) {
  guild <- spatial_combinations$guild[i]
  scenario <- spatial_combinations$scenario[i]
  scenario_info <- spatial_combinations[i, c("guild", "scenario", "period", "SSP"),
                                        drop = FALSE]
  message("Cluster-direction classification: ", guild, " / ", scenario)
  
  L <- terra::rast(spatial_combinations$lisa_file[i])
  C_all <- terra::rast(change_files[[guild]])
  if (terra::nlyr(L) != 1L || anyDuplicated(names(C_all)) ||
      !scenario %in% names(C_all)) {
    stop("Check the LISA raster and the named richness-change layer.")
  }
  C <- C_all[[scenario]]
  terra::compareGeom(L, C, stopOnError = TRUE)
  v_lisa <- terra::values(L, mat = FALSE)
  v_change <- terra::values(C, mat = FALSE)
  if (!identical(is.na(v_lisa), is.na(v_change)) ||
      any(!v_lisa[!is.na(v_lisa)] %in% lisa_lookup$value) ||
      any(!is.finite(v_change[!is.na(v_change)]))) {
    stop("Invalid LISA codes or inconsistent valid-cell masks: ", guild, " / ", scenario)
  }
  
  cluster_direction <- terra::rast(L, vals = NA_integer_)
  names(cluster_direction) <- "SCC_class"
  valid <- !is.na(L) & !is.na(C)
  cluster_direction[valid & L == 5 & C > 0] <- 1L  # H.P
  cluster_direction[valid & L == 5 & C == 0] <- 2L # H.0
  cluster_direction[valid & L == 5 & C < 0] <- 3L  # H.N
  cluster_direction[valid & L == 1 & C > 0] <- 4L  # C.P
  cluster_direction[valid & L == 1 & C == 0] <- 5L # C.0
  cluster_direction[valid & L == 1 & C < 0] <- 6L  # C.N
  cluster_direction[valid & (L %in% c(2, 3, 4))] <- 7L # other
  # Zero means exactly zero richness change; no epsilon/tolerance is introduced.
  v_class <- terra::values(cluster_direction, mat = FALSE)
  if (!identical(is.na(v_class), is.na(v_change)) ||
      any(!v_class[!is.na(v_class)] %in% class_lut$value)) {
    stop("Incomplete cluster-direction classification: ", guild, " / ", scenario)
  }
  terra::writeRaster(cluster_direction,
                     spatial_combinations$cluster_direction_file[i],
                     overwrite = overwrite_output, datatype = "INT1U")
  
  basins_local <- basins
  if (!terra::same.crs(basins_local, cluster_direction)) {
    basins_local <- terra::project(basins_local, terra::crs(cluster_direction))
  }
  cell_area <- terra::cellSize(cluster_direction, unit = "km")
  area_layers <- NULL
  for (j in seq_len(nrow(class_lut))) {
    class_area <- terra::ifel(
      cluster_direction == class_lut$value[j], cell_area, 0
    )
    names(class_area) <- paste0("area_", class_lut$class[j], "_km2")
    area_layers <- if (is.null(area_layers)) class_area else c(area_layers, class_area)
  }
  # Unlike Section 3, the source uses exact cell coverage for these seven classes.
  area_by_basin <- terra::extract(
    area_layers, basins_local, fun = sum, na.rm = TRUE, exact = TRUE
  )
  basin_order <- match(basin_attributes$polygon_id, area_by_basin$ID)
  if (anyNA(basin_order) || anyDuplicated(area_by_basin$ID)) {
    stop("Cannot match extracted cross-classification areas to basin features.")
  }
  # Preserve every input polygon; do not replace missing areas with zero.
  cluster_area_rows[[i]] <- cbind(
    basin_attributes, area_by_basin[basin_order, -1, drop = FALSE],
    scenario_info[rep(1L, nrow(basin_attributes)), , drop = FALSE]
  )
  cluster_count_rows[[i]] <- data.frame(
    guild = guild, scenario = scenario, value = class_lut$value,
    class = class_lut$class,
    n_cells = tabulate(as.integer(v_class[!is.na(v_class)]), nbins = 7L)
  )
  rm(L, C, C_all, cluster_direction, area_layers, class_area)
}
cluster_area_all <- do.call(rbind, cluster_area_rows)
rownames(cluster_area_all) <- NULL
cluster_counts <- do.call(rbind, cluster_count_rows)
write.csv(cluster_area_all, file.path(table_dir, "cluster_direction_area_by_basin.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(cluster_counts, file.path(table_dir, "cluster_direction_class_counts.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 5. Quantify temporal compositional reassembly within each guild ----------
# These inputs must be the BASIN-CONSTRAINED, study-region species stacks from script 05-06, 
# not individual unconstrained predictions from script 02.
fmap_current <- terra::rast(stack_files[[current_scenario]])
fmap_future <- terra::rast(stack_files[[reassembly_scenario]])
terra::compareGeom(fmap_current[[1]], fmap_future[[1]], stopOnError = TRUE)
for (layer_names in list(names(fmap_current), names(fmap_future))) {
  if (anyDuplicated(layer_names) || length(layer_names) != nrow(fish.info) ||
      !setequal(layer_names, fish.info$fcode)) {
    stop("Species stack layer names must uniquely match all fcode values. ",
         "Do not relabel layers by position to hide a mismatch.")
  }
}
# Select by fcode, allowing a different storage order but never assuming identity.
reassembly_files <- setNames(character(nrow(guild_lookup)), guild_lookup$guild)
reassembly_summary_rows <- transition_count_rows <- vector("list", nrow(guild_lookup))
transition_levels <- c("both_empty", "complete_loss", "complete_gain",
                       "complete_replacement", "partial_or_no_change")

for (g in seq_len(nrow(guild_lookup))) {
  guild <- guild_lookup$guild[g]
  codes <- fish.info$fcode[fish.info$thermal.guild == guild_lookup$thermal.guild[g]]
  message("Compositional reassembly: ", guild, " / ", reassembly_scenario)
  combination_row <- which(spatial_combinations$guild == guild &
                             spatial_combinations$scenario == reassembly_scenario)
  if (!length(codes) || length(combination_row) != 1L) {
    stop("Missing guild species or a unique spatial-analysis combination.")
  }
  change_stack <- terra::rast(change_files[[guild]])
  change_map <- change_stack[[reassembly_scenario]]
  lisa_map <- terra::rast(spatial_combinations$lisa_file[combination_row])
  cc_map <- terra::rast(spatial_combinations$cluster_direction_file[combination_row])
  terra::compareGeom(fmap_current[[1]], change_map, stopOnError = TRUE)
  terra::compareGeom(fmap_current[[1]], lisa_map, stopOnError = TRUE)
  terra::compareGeom(fmap_current[[1]], cc_map, stopOnError = TRUE)
  vchange <- terra::values(change_map, mat = FALSE)
  v_lisa <- terra::values(lisa_map, mat = FALSE)
  v_class <- terra::values(cc_map, mat = FALSE)
  
  # Preserve the source domain, including unchanged cells and both-empty cells.
  # Never filter only HH/LL, positive changes, or rows with defined beta components.
  cell_id <- which(!is.na(vchange))
  if (!length(cell_id) || anyNA(v_lisa[cell_id]) || anyNA(v_class[cell_id])) {
    stop("The change domain is empty or missing spatial classifications.")
  }
  current_values <- terra::values(fmap_current[[codes]], mat = TRUE)
  future_values <- terra::values(fmap_future[[codes]], mat = TRUE)
  current_values <- current_values[cell_id, , drop = FALSE]
  future_values <- future_values[cell_id, , drop = FALSE]
  if (!identical(dim(current_values), dim(future_values)) ||
      nrow(current_values) != length(cell_id)) {
    stop("Species matrices must have equal dimensions and aligned raster cells.")
  }
  n_na_current <- sum(is.na(current_values))
  n_na_future <- sum(is.na(future_values))
  # Source convention: species NAs inside the change domain represent assemblage absence (e.g., outside species-specific basin constraints). 
  # This is NOT a generic missing-data imputation rule; input preparation must justify it.
  current_values[is.na(current_values)] <- 0
  future_values[is.na(future_values)] <- 0
  if (!all(current_values %in% c(0, 1)) || !all(future_values %in% c(0, 1))) {
    stop("Species matrices contain values other than 0, 1, or the replaced NAs.")
  }
  
  # Vectorized evaluation of the SAME source presence/absence formulas.
  shared <- rowSums(current_values == 1 & future_values == 1)
  Loss <- rowSums(current_values == 1 & future_values == 0)
  Gain <- rowSums(current_values == 0 & future_values == 1)
  rich_current <- rowSums(current_values)
  rich_future <- rowSums(future_values)
  sorensen_denominator <- 2 * shared + Loss + Gain
  turnover_denominator <- shared + pmin(Loss, Gain)
  
  beta_total <- rep(NA_real_, length(cell_id))
  turnover <- rep(NA_real_, length(cell_id))
  ok_total <- sorensen_denominator > 0
  ok_turnover <- turnover_denominator > 0
  beta_total[ok_total] <- (Loss[ok_total] + Gain[ok_total]) /
    sorensen_denominator[ok_total]
  turnover[ok_turnover] <- pmin(Loss[ok_turnover], Gain[ok_turnover]) /
    turnover_denominator[ok_turnover]
  nestedness <- beta_total - turnover
  # Both empty: all three metrics are NA.
  # Complete loss/gain: beta_total = 1, turnover and nestedness remain NA.
  # Both nonempty with no shared species: turnover = 1, nestedness = 0.
  # Undefined components are deliberately NOT changed to zero.
  delta_calculated <- rich_future - rich_current
  if (!isTRUE(all.equal(as.numeric(vchange[cell_id]), as.numeric(delta_calculated))) ||
      !all(delta_calculated == Gain - Loss) ||
      !all(rich_current == shared + Loss) ||
      !all(rich_future == shared + Gain)) {
    stop("Richness and compositional counts disagree for guild ", guild,
         ". Check scenario identity, species names, NA meaning, and input masks.")
  }
  
  transition <- rep("partial_or_no_change", length(cell_id))
  transition[rich_current == 0 & rich_future == 0] <- "both_empty"
  transition[rich_current > 0 & rich_future == 0] <- "complete_loss"
  transition[rich_current == 0 & rich_future > 0] <- "complete_gain"
  transition[shared == 0 & rich_current > 0 & rich_future > 0] <- "complete_replacement"
  xy <- terra::xyFromCell(fmap_current, cell_id)
  
  # Original analytical column names are retained for downstream compatibility;
  # guild and scenario identifiers are appended to make standalone CSVs explicit.
  reassembly <- data.frame(
    cell = cell_id, x = xy[, 1], y = xy[, 2],
    CC.class = v_class[cell_id], LISA.spot = v_lisa[cell_id],
    change = vchange[cell_id], rich.current = rich_current, rich.future = rich_future,
    turnover = turnover, nestedness = nestedness, beta_total = beta_total,
    shared = shared, Loss = Loss, Gain = Gain, transition = transition,
    guild = guild, current_scenario = current_scenario,
    future_scenario = reassembly_scenario, stringsAsFactors = FALSE
  )
  reassembly_files[[guild]] <- file.path(
    reassembly_dir, paste0("compositional_reassembly_", guild, "_",
                           reassembly_scenario, ".csv")
  )
  write.csv(reassembly, reassembly_files[[guild]], row.names = FALSE,
            fileEncoding = output_encoding, na = "NA")
  reassembly_summary_rows[[g]] <- data.frame(
    guild = guild, current_scenario = current_scenario,
    future_scenario = reassembly_scenario, n_species = length(codes),
    n_cells = length(cell_id), n_current_species_NAs_as_absence = n_na_current,
    n_future_species_NAs_as_absence = n_na_future,
    n_beta_total_NA = sum(is.na(beta_total)), n_turnover_NA = sum(is.na(turnover)),
    n_nestedness_NA = sum(is.na(nestedness)), richness_consistency_passed = TRUE,
    output_file = reassembly_files[[guild]]
  )
  transition_count_rows[[g]] <- data.frame(
    guild = guild, future_scenario = reassembly_scenario,
    transition = transition_levels,
    n_cells = as.integer(table(factor(transition, levels = transition_levels)))
  )
  rm(current_values, future_values, change_stack, change_map, lisa_map, cc_map, reassembly)
}

# ---- 6. Export summaries, provenance, and run information ---------------------
reassembly_summary <- do.call(rbind, reassembly_summary_rows)
transition_counts <- do.call(rbind, transition_count_rows)
write.csv(reassembly_summary, file.path(table_dir, "reassembly_summary.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(transition_counts, file.path(table_dir, "reassembly_transition_counts.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

input_manifest <- data.frame(
  role = c("guild_metadata", "basin_polygons",
           rep("richness_raster", length(richness_files)),
           rep("richness_change_stack", length(change_files)),
           rep("species_binary_stack", length(stack_files))),
  file = c(guild_file, basin_file, as.vector(richness_files),
           unname(change_files), unname(stack_files)), stringsAsFactors = FALSE
)
write.csv(input_manifest, file.path(table_dir, "input_manifest.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

analysis_settings <- list(
  n_neighbors = n_neighbors, knn_longlat = knn_longlat,
  weights_style = "W", symmetric_neighbors = FALSE,
  localmoran_conditional = TRUE, localmoran_mlvar = TRUE,
  p_value_method = "two-sided analytical normal approximation from Z.Ii",
  p_adjust_method = p_adjust_method, p_adjust_scope = "one guild-scenario raster",
  significance_level = significance_level,
  quadrant_rule = "scale(change) and scale(lag(change)); source boundary assignment order",
  lisa_area_exact = FALSE, cluster_direction_area_exact = TRUE,
  n_incomplete_lisa_area_rows = n_incomplete_lisa_area_rows,
  current_scenario = current_scenario, reassembly_scenario = reassembly_scenario,
  domain_rule = "non-NA richness-change cells; retain both-empty and other cells",
  species_NA_rule = "replace with 0 inside the change domain, as in source",
  undefined_beta_rule = "retain NA components; do not impute zero",
  basin_code_field = basin_code_field, basin_group_field = basin_group_field,
  save_lisa_diagnostics = save_lisa_diagnostics
)

# Save plain R tables/settings/file paths, not external-pointer raster objects.
save(
  guild_lookup, future_scenarios, spatial_combinations, lisa_lookup, class_lut,
  basin_attributes, lisa_area_all, cluster_area_all, lisa_counts, cluster_counts,
  spatial_summary, reassembly_summary, transition_counts, reassembly_files,
  input_manifest, analysis_settings,
  file = file.path(output_dir, "cluster_direction_reassembly.RData")
)
writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "session_info.txt"))
version_packages <- c(required_packages, "dbscan")
package_versions <- data.frame(
  package = c("R", version_packages),
  version = c(as.character(getRversion()), vapply(version_packages, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p))
    else "not installed (optional)"
  }, character(1)))
)
write.csv(package_versions, file.path(output_dir, "package_versions.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
expected_outputs <- c(spatial_combinations$lisa_file,
                      spatial_combinations$cluster_direction_file, reassembly_files)
if (any(!file.exists(expected_outputs))) {
  stop("Some expected analytical outputs were not created: ",
       paste(expected_outputs[!file.exists(expected_outputs)], collapse = ", "))
}
message("Cluster-direction and compositional reassembly analyses complete: ", output_dir)
