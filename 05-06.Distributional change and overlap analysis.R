# =============================================================================
# 05-06. Distributional change and overlap analysis
# =============================================================================
# Run the numbered sections sequentially from the project root used in 01-04.
# Source workflow: basin-constrained distributions -> richness and change ->
# spatial overlap -> suitable-habitat area and change -> far-future Gaussian LMM.
# No figure drawing or graphics-device calls are included. GeoTIFFs are retained
# as analytical data, including the species stacks needed by later analyses.
# Settings are defined here; there is no configuration file or workflow wrapper.

# ---- 0. Paths and settings ---------------------------------------------------
project_dir <- getwd()
raw_dir <- file.path(project_dir, "data", "raw")
occurrence_file <- file.path(
  project_dir, "outputs", "data_preprocessing", "combined_occurrences.RData"
)
guild_file <- file.path(
  project_dir, "outputs", "thermal_guild_classification", "species_thermal_guilds.csv"
)
sdm_map_dir <- file.path(
  project_dir, "outputs", "species_distribution_modelling", "prediction_maps"
)
basin_file <- file.path(raw_dir, "river_basins.shp")

# Input fields in the supplied river-basin vector. 
# These represent two distinct spatial units: river_bas0 for dispersal constraints; BASIN_ID for area and LMM.
basin_id_field <- "BASIN_ID"
basin_name_field <- "BASIN_NM"
basin_group_field <- "river_bas0"

# The two placeholders are the species fcode and the scenario, respectively.
# For archived original SDM maps, use "2.binary_%s_MaxEnt_%s.tif" instead.
binary_file_template <- "binary_%s_MaxEnt_%s.tif"
input_encoding <- "EUC-KR"  
output_encoding <- "EUC-KR"

output_dir <- file.path(project_dir, "outputs", "distributional_change_overlap")
constrained_dir <- file.path(output_dir, "constrained_maps")
stack_dir <- file.path(output_dir, "species_stacks")
richness_dir <- file.path(output_dir, "richness_maps")
change_dir <- file.path(output_dir, "richness_change_maps")
reference_dir <- file.path(output_dir, "spatial_reference")
table_dir <- file.path(output_dir, "tables")
lmm_dir <- file.path(output_dir, "lmm")
overwrite_output <- FALSE  # Set TRUE deliberately before replacing old results.

# Retain the source's spatial exclusions. 
# The Unicode escape is the original Korean attribute value for Jeju; change this setting if labels were translated.
excluded_basin_groups <- "\uc81c\uc8fc"
analysis_longitude_max <- 130.5

# Use explicit scenario names, never alphabetical file order, for subtraction.
scenarios <- data.frame(
  scenario = c("current_normal", paste0("nearf_ssp", c(126, 370, 585)),
               paste0("farf_ssp", c(126, 370, 585))),
  period = c("current", rep("near-future", 3), rep("far-future", 3)),
  period_code = c("current", rep("nearf", 3), rep("farf", 3)),
  sce = c("normal", rep(c("ssp126", "ssp370", "ssp585"), 2)),
  SSP = c("current", rep(c("SSP126", "SSP370", "SSP585"), 2)),
  stringsAsFactors = FALSE
)
current_scenario <- "current_normal"
future_scenarios <- scenarios$scenario[scenarios$scenario != current_scenario]
guild_lookup <- data.frame(
  thermal.guild = c("cold-affinity", "intermediate-affinity", "warm-affinity"),
  guild = c("cold", "intermediate", "warm"), stringsAsFactors = FALSE
)
guild_levels <- guild_lookup$thermal.guild
ssp_levels <- c("SSP126", "SSP370", "SSP585")
run_far_future_lmm <- TRUE  # FALSE runs only the spatial analyses (Sections 1-8).
singularity_tolerance <- 1e-4

# Install missing packages once in your local R console, not inside the workflow:
# install.packages(c("terra", "lme4", "lmerTest", "emmeans", "pbkrtest"))
required_packages <- "terra"
if (run_far_future_lmm) {
  required_packages <- c(required_packages, "lme4", "lmerTest", "emmeans", "pbkrtest")
}
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install the required R package before continuing: ", package)
  }
}
library(terra)

# ---- 1. Read and check metadata, occurrences, and basin polygons --------------
required_files <- c(guild_file, occurrence_file, basin_file)
if (any(!file.exists(required_files))) {
  stop("Missing input files: ", paste(required_files[!file.exists(required_files)],
                                      collapse = ", "))
}
if (!dir.exists(sdm_map_dir)) stop("SDM map directory not found: ", sdm_map_dir)
if (!overwrite_output && dir.exists(output_dir) &&
    length(list.files(output_dir, recursive = TRUE)) > 0L) {
  stop("output_dir already contains files. Use a new directory or explicitly ",
       "set overwrite_output <- TRUE. No existing files have been deleted.")
}
for (directory in c(output_dir, constrained_dir, stack_dir, richness_dir,
                    change_dir, reference_dir, table_dir, lmm_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

fish.info <- read.csv(guild_file, fileEncoding = input_encoding,
                      stringsAsFactors = FALSE, check.names = FALSE)
# Accept original metadata headers, while preserving the species-name values.
names(fish.info)[names(fish.info) == "\uad6d\uba85"] <- "common_name"
names(fish.info)[names(fish.info) == "\ud559\uba85"] <- "scientific_name"
required_fields <- c("fcode", "common_name", "scientific_name", "thermal.guild",
                     "PC1.centroid", "PC2.centroid")
if (!all(required_fields %in% names(fish.info)) || !nrow(fish.info) ||
    anyDuplicated(names(fish.info))) {
  stop("Guild metadata require unique headers and these fields: ",
       paste(required_fields, collapse = ", "))
}
if (anyNA(fish.info[, required_fields])) stop("Missing values in guild metadata.")
for (field in c("fcode", "common_name", "scientific_name", "thermal.guild")) {
  fish.info[[field]] <- as.character(fish.info[[field]])
  if (any(!nzchar(trimws(fish.info[[field]])))) {
    stop("Blank metadata values in: ", field)
  }
}
if (anyDuplicated(fish.info$fcode) || anyDuplicated(fish.info$common_name) ||
    any(!grepl("^[A-Za-z0-9-]+$", fish.info$fcode))) {
  stop("Species codes and common-name keys must be unique; use filename-safe codes.")
}
if (!setequal(unique(fish.info$thermal.guild), guild_levels)) {
  stop("The input must represent each of the three expected thermal guilds.")
}
if (!is.numeric(fish.info$PC1.centroid) || !is.numeric(fish.info$PC2.centroid) ||
    any(!is.finite(as.matrix(fish.info[, c("PC1.centroid", "PC2.centroid")])))) {
  stop("PC1.centroid and PC2.centroid must be finite numeric values.")
}

# Use combined occurrences BEFORE environmental filtering/thinning, as in source.
occurrence_data <- new.env(parent = emptyenv())
load(occurrence_file, envir = occurrence_data)
if (!exists("DB.fish", envir = occurrence_data, inherits = FALSE)) {
  stop("The occurrence RData must contain the data frame DB.fish.")
}
DB.fish <- occurrence_data$DB.fish
rm(occurrence_data)
if (!is.data.frame(DB.fish) ||
    !all(c("species", "lon", "lat") %in% names(DB.fish))) {
  stop("DB.fish requires species, lon, and lat.")
}
DB.fish <- DB.fish[DB.fish$species %in% fish.info$common_name, , drop = FALSE]
if (!all(fish.info$common_name %in% DB.fish$species)) {
  stop("Some guild-table species have no matching combined occurrence records.")
}
if (!is.numeric(DB.fish$lon) || !is.numeric(DB.fish$lat) ||
    any(!is.finite(DB.fish$lon)) || any(!is.finite(DB.fish$lat)) ||
    any(abs(DB.fish$lon) > 180) || any(abs(DB.fish$lat) > 90)) {
  stop("Retained occurrences require valid WGS84 longitude/latitude coordinates.")
}

rb.shp <- terra::vect(basin_file)
source_basin_fields <- c(basin_id_field, basin_name_field, basin_group_field)
if (anyDuplicated(source_basin_fields) ||
    !all(source_basin_fields %in% names(rb.shp)) ||
    !nrow(rb.shp) || !nzchar(terra::crs(rb.shp)) ||
    !all(terra::geomtype(rb.shp) == "polygons")) {
  stop("The basin vector requires polygon geometry, a defined CRS, and the ",
       "three distinct fields named in Section 0.")
}
# Keep a compact, consistent attribute schema without translating label values.
rb.shp <- rb.shp[, source_basin_fields]
names(rb.shp) <- c("BASIN_ID", "BASIN_NM", "river_bas0")
basin_attributes <- as.data.frame(rb.shp)
if (anyNA(basin_attributes) || anyDuplicated(basin_attributes$BASIN_ID)) {
  stop("Basin attributes must be complete, with one feature per BASIN_ID. ",
       "Multipart geometry within a feature is supported.")
}
for (field in names(basin_attributes)) {
  if (any(!nzchar(trimws(as.character(basin_attributes[[field]]))))) {
    stop("Blank basin attribute values in: ", field)
  }
}

# The source's longitude cutoff assumes geographic coordinates. 
# Transform only vectors when necessary; binary rasters are never resampled or reprojected.
rb.shp <- terra::project(rb.shp, "EPSG:4326")
rb.shp.used <- rb.shp[!rb.shp$river_bas0 %in% excluded_basin_groups, ]
if (!nrow(rb.shp.used) || !is.finite(analysis_longitude_max)) {
  stop("No analysis basins remain, or the longitude cutoff is invalid.")
}
analysis_extent <- terra::ext(rb.shp.used)
if (analysis_extent[1] >= analysis_longitude_max) {
  stop("The longitude cutoff lies outside the analysis region.")
}
analysis_extent[2] <- analysis_longitude_max
rb.shp.used2 <- terra::crop(rb.shp.used, analysis_extent)
if (!nrow(rb.shp.used2)) stop("No basin polygons remain after spatial cropping.")

# Exact filenames prevent a short species code from matching another species.
input_map_files <- constrained_files <- matrix(
  NA_character_, nrow = nrow(fish.info), ncol = nrow(scenarios),
  dimnames = list(fish.info$fcode, scenarios$scenario)
)
for (j in seq_len(nrow(scenarios))) {
  scenario <- scenarios$scenario[j]
  input_map_files[, j] <- file.path(
    sdm_map_dir, scenario, sprintf(binary_file_template, fish.info$fcode, scenario)
  )
  constrained_files[, j] <- file.path(
    constrained_dir, scenario,
    paste0("constrained_", fish.info$fcode, "_MaxEnt_", scenario, ".tif")
  )
}
if (any(!file.exists(input_map_files))) {
  stop("Missing species/scenario binary maps: ",
       paste(input_map_files[!file.exists(input_map_files)], collapse = ", "))
}
input_template <- terra::rast(input_map_files[1, current_scenario])
if (terra::nlyr(input_template) != 1L || !nzchar(terra::crs(input_template))) {
  stop("Each input binary map must be single-layer with a defined CRS.")
}
write.csv(fish.info, file.path(table_dir, "analysed_species.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(scenarios, file.path(table_dir, "scenario_lookup.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 2. Identify the occupied basin groups for each species -------------------
species_basin_groups <- setNames(vector("list", nrow(fish.info)), fish.info$fcode)
basin_assignment_rows <- basin_assignment_summary <- vector("list", nrow(fish.info))
for (i in seq_len(nrow(fish.info))) {
  species_records <- DB.fish[DB.fish$species == fish.info$common_name[i], , drop = FALSE]
  species_xy <- unique(species_records[, c("lon", "lat"), drop = FALSE])
  species_points <- terra::vect(species_xy, geom = c("lon", "lat"), crs = "EPSG:4326")
  basin_hits <- terra::extract(rb.shp, species_points)
  occupied_groups <- sort(unique(as.character(
    basin_hits$river_bas0[!is.na(basin_hits$river_bas0)]
  )))
  if (!length(occupied_groups)) {
    stop("No occurrence intersects a basin for ", fish.info$fcode[i],
         ". Check coordinates and the basin vector; no nearest-basin assignment is made.")
  }
  species_basin_groups[[i]] <- occupied_groups
  basin_assignment_rows[[i]] <- data.frame(
    fcode = fish.info$fcode[i], common_name = fish.info$common_name[i],
    scientific_name = fish.info$scientific_name[i], river_bas0 = occupied_groups
  )
  basin_assignment_summary[[i]] <- data.frame(
    fcode = fish.info$fcode[i], n_records = nrow(species_records),
    n_unique_locations = nrow(species_xy), n_basin_groups = length(occupied_groups),
    n_selected_polygons = sum(rb.shp$river_bas0 %in% occupied_groups)
  )
}
species_basin_assignments <- do.call(rbind, basin_assignment_rows)
basin_assignment_summary <- do.call(rbind, basin_assignment_summary)
write.csv(species_basin_assignments, file.path(table_dir, "species_basin_groups.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(basin_assignment_summary,
          file.path(table_dir, "species_basin_assignment_summary.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 3. Apply the fixed, occurrence-based dispersal constraint ----------------
# A species retains all polygons sharing an occupied river_bas0 label, not only individual BASIN_ID polygons containing a record. 
# Apply the same mask in all periods. Masking preserves the input extent and makes other cells NA.
constrained_manifest <- vector("list", nrow(fish.info) * nrow(scenarios))
manifest_index <- 0L
for (i in seq_len(nrow(fish.info))) {
  message("Constraining distributions: ", fish.info$fcode[i])
  species_basins <- rb.shp[rb.shp$river_bas0 %in% species_basin_groups[[i]], ]
  species_basins <- terra::project(species_basins, terra::crs(input_template))
  for (j in seq_len(nrow(scenarios))) {
    scenario <- scenarios$scenario[j]
    binary_map <- terra::rast(input_map_files[i, j])
    if (terra::nlyr(binary_map) != 1L) stop("Input map has multiple layers: ", input_map_files[i, j])
    terra::compareGeom(input_template, binary_map, stopOnError = TRUE)
    invalid_values <- terra::global(
      !is.na(binary_map) & binary_map != 0 & binary_map != 1, "sum", na.rm = TRUE
    )[1, 1]
    if (!is.finite(invalid_values) || invalid_values > 0) {
      stop("Input map contains values other than 0, 1, or NA: ", input_map_files[i, j])
    }
    dir.create(dirname(constrained_files[i, j]), recursive = TRUE, showWarnings = FALSE)
    constrained_map <- terra::mask(binary_map, species_basins, touches = TRUE)
    names(constrained_map) <- fish.info$fcode[i]
    terra::writeRaster(constrained_map, constrained_files[i, j], overwrite = overwrite_output)
    manifest_index <- manifest_index + 1L
    constrained_manifest[[manifest_index]] <- data.frame(
      fcode = fish.info$fcode[i], scenario = scenario,
      input_file = input_map_files[i, j], constrained_file = constrained_files[i, j]
    )
  }
}
constrained_manifest <- do.call(rbind, constrained_manifest)
write.csv(constrained_manifest, file.path(table_dir, "constrained_map_manifest.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 4. Create analysis-region species stacks and guild richness --------------
stack_files <- setNames(
  file.path(stack_dir, paste0("species_binary_", scenarios$scenario, ".tif")),
  scenarios$scenario
)
richness_files <- matrix(
  NA_character_, nrow = nrow(guild_lookup), ncol = nrow(scenarios),
  dimnames = list(guild_lookup$guild, scenarios$scenario)
)
richness_summary <- vector("list", nrow(guild_lookup) * nrow(scenarios))
summary_index <- 0L
analysis_basins <- terra::project(rb.shp.used2, terra::crs(input_template))
for (j in seq_len(nrow(scenarios))) {
  scenario <- scenarios$scenario[j]
  message("Stacking species and calculating richness: ", scenario)
  species_stack <- terra::rast(constrained_files[, j])
  names(species_stack) <- fish.info$fcode
  species_stack <- terra::crop(species_stack, analysis_basins, mask = TRUE, touches = TRUE)
  if (j == 1L) {
    analysis_template <- species_stack[[1]]
  } else {
    terra::compareGeom(analysis_template, species_stack[[1]], stopOnError = TRUE)
  }
  terra::writeRaster(species_stack, stack_files[[scenario]], overwrite = overwrite_output)
  for (g in seq_len(nrow(guild_lookup))) {
    guild <- guild_lookup$guild[g]
    codes <- fish.info$fcode[fish.info$thermal.guild == guild_lookup$thermal.guild[g]]
    # Do not impose a new habitat mask after summing or change the interpretation of species NAs.
    richness <- terra::app(species_stack[[codes]], fun = sum, na.rm = TRUE)
    names(richness) <- scenario
    richness_files[g, j] <- file.path(richness_dir, paste0("richness_", guild, "_", scenario, ".tif"))
    terra::writeRaster(richness, richness_files[g, j], overwrite = overwrite_output)
    summary_index <- summary_index + 1L
    richness_summary[[summary_index]] <- data.frame(
      guild = guild, scenario = scenario, n_species = length(codes),
      min_richness = terra::global(richness, "min", na.rm = TRUE)[1, 1],
      max_richness = terra::global(richness, "max", na.rm = TRUE)[1, 1]
    )
  }
}
richness_summary <- do.call(rbind, richness_summary)
write.csv(richness_summary, file.path(table_dir, "richness_ranges.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 5. Calculate guild richness change: future minus current -----------------
change_files <- setNames(character(nrow(guild_lookup)), guild_lookup$guild)
for (g in seq_len(nrow(guild_lookup))) {
  guild <- guild_lookup$guild[g]
  current_richness <- terra::rast(richness_files[g, current_scenario])
  change_stack <- NULL
  for (scenario in future_scenarios) {
    future_richness <- terra::rast(richness_files[g, scenario])
    terra::compareGeom(current_richness, future_richness, stopOnError = TRUE)
    richness_change <- future_richness - current_richness
    names(richness_change) <- scenario
    change_stack <- if (is.null(change_stack)) richness_change else c(change_stack, richness_change)
  }
  change_files[[guild]] <- file.path(change_dir, paste0("richness_change_", guild, ".tif"))
  terra::writeRaster(change_stack, change_files[[guild]], overwrite = overwrite_output)
}

# ---- 6. Spatial overlap between guilds and within guilds -----------------------
# One reusable numerical helper replaces the source's duplicated pairwise wrappers. 
# All species/scenario iteration remains in the top-level script.
spatial_metrics <- function(r1, r2, area_unit = "km") {
  if (terra::nlyr(r1) != 1L || terra::nlyr(r2) != 1L) {
    stop("Spatial overlap requires two single-layer rasters.")
  }
  terra::compareGeom(r1, r2, stopOnError = TRUE)
  area <- terra::cellSize(r1, unit = area_unit, mask = FALSE)
  v1 <- terra::values(r1, mat = FALSE)
  v2 <- terra::values(r2, mat = FALSE)
  va <- terra::values(area, mat = FALSE)
  keep <- !(is.na(v1) | is.na(v2) | is.na(va))
  v1 <- v1[keep]; v2 <- v2[keep]; va <- va[keep]
  if (any(!is.finite(v1)) || any(!is.finite(v2)) || any(!is.finite(va)) ||
      any(v1 < 0) || any(v2 < 0) || any(va <= 0)) {
    stop("Overlap inputs must contain finite, non-negative richness and positive cell areas.")
  }
  v1w <- v1 * va
  v2w <- v2 * va
  s1 <- sum(v1w); s2 <- sum(v2w)
  out <- data.frame(
    D = NA_real_, I = NA_real_, BC = NA_real_, D_abs = NA_real_, I_abs = NA_real_,
    n_cells = length(v1w), sum_r1 = s1, sum_r2 = s2, area_unit = area_unit,
    status = "ok", stringsAsFactors = FALSE
  )
  # Preserve undefined overlap when either total is zero, but always return a
  # complete row rather than the source's incomplete list in this edge case.
  if (!length(v1w) || s1 == 0 || s2 == 0) {
    out$status <- if (!length(v1w)) "no_shared_cells" else "zero_total_richness"
    return(out)
  }
  p <- v1w / s1
  q <- v2w / s2
  out$D <- 1 - 0.5 * sum(abs(p - q))
  out$I <- sum(sqrt(p * q))
  out$BC <- 2 * sum(pmin(v1w, v2w)) / sum(v1w + v2w)
  total_ratio <- min(s1, s2) / max(s1, s2)
  out$D_abs <- out$D * total_ratio
  out$I_abs <- out$I * total_ratio
  out
}

# 6.1. Between-guild overlap within each climate condition.
# Source convention: divide each guild map by ITS OWN SPATIAL MAXIMUM.
guild_pairs <- t(utils::combn(seq_len(nrow(guild_lookup)), 2))
inter_rows <- vector("list", nrow(scenarios) * nrow(guild_pairs))
row_index <- 0L
for (j in seq_len(nrow(scenarios))) {
  scenario <- scenarios$scenario[j]
  guild_maps <- vector("list", nrow(guild_lookup))
  map_maxima <- numeric(nrow(guild_lookup))
  for (g in seq_len(nrow(guild_lookup))) {
    guild_maps[[g]] <- terra::rast(richness_files[g, j])
    map_maxima[g] <- terra::global(guild_maps[[g]], "max", na.rm = TRUE)[1, 1]
    # Zero maps are not divided by zero; the helper reports undefined overlap.
    if (is.finite(map_maxima[g]) && map_maxima[g] > 0) {
      guild_maps[[g]] <- guild_maps[[g]] / map_maxima[g]
    }
  }
  for (k in seq_len(nrow(guild_pairs))) {
    a <- guild_pairs[k, 1]; b <- guild_pairs[k, 2]
    row_index <- row_index + 1L
    inter_rows[[row_index]] <- cbind(
      data.frame(r1 = guild_lookup$guild[a], r2 = guild_lookup$guild[b]),
      spatial_metrics(guild_maps[[a]], guild_maps[[b]]),
      data.frame(period = scenarios$period_code[j], sce = scenarios$sce[j],
                 scenario = scenario, normalization = "per_map_spatial_maximum",
                 original_max_r1 = map_maxima[a], original_max_r2 = map_maxima[b])
    )
  }
}
DB.comp.inter <- do.call(rbind, inter_rows)

# 6.2. Within-guild current-future overlap using UNNORMALIZED richness maps.
intra_rows <- vector("list", nrow(guild_lookup) * length(future_scenarios))
row_index <- 0L
for (g in seq_len(nrow(guild_lookup))) {
  current_richness <- terra::rast(richness_files[g, current_scenario])
  for (scenario in future_scenarios) {
    future_richness <- terra::rast(richness_files[g, scenario])
    row_index <- row_index + 1L
    intra_rows[[row_index]] <- cbind(
      data.frame(r1 = current_scenario, r2 = scenario),
      spatial_metrics(current_richness, future_richness),
      data.frame(type = guild_lookup$guild[g], normalization = "none")
    )
  }
}
DB.comp.intra <- do.call(rbind, intra_rows)
write.csv(DB.comp.inter, file.path(table_dir, "overlap_between_guilds.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(DB.comp.intra, file.path(table_dir, "overlap_within_guilds.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 7. Calculate species-specific suitable-habitat area by basin -------------
map.current <- terra::rast(stack_files[[current_scenario]])
if (!identical(names(map.current), fish.info$fcode)) {
  stop("Species-stack layer names/order do not match the guild table.")
}
rb.area <- terra::project(rb.shp.used2, terra::crs(map.current))
# A one-to-one numeric index supports either numeric or character BASIN_ID
# values without changing the spatial grouping. Original IDs remain in tables.
basin.lookup <- data.frame(
  basin_index = seq_len(nrow(rb.area)), BASIN_ID = as.character(rb.area$BASIN_ID),
  BASIN_NM = as.character(rb.area$BASIN_NM), river_bas0 = as.character(rb.area$river_bas0),
  stringsAsFactors = FALSE
)
if (anyDuplicated(basin.lookup$BASIN_ID)) stop("BASIN_ID is not unique after cropping.")
rb.area$basin_index <- basin.lookup$basin_index
# Retain rasterize + zonal, not fractional polygon-cell area extraction.
basin.rast <- terra::rasterize(rb.area, map.current[[1]], field = "basin_index",
                               touches = FALSE, background = NA)
names(basin.rast) <- "basin_index"
cell.area <- terra::cellSize(map.current[[1]], unit = "km")
names(cell.area) <- "cell_area_km2"
terra::writeRaster(basin.rast, file.path(reference_dir, "basin_zones.tif"),
                   overwrite = overwrite_output)
terra::writeRaster(cell.area, file.path(reference_dir, "cell_area_km2.tif"),
                   overwrite = overwrite_output)
write.csv(basin.lookup, file.path(reference_dir, "basin_lookup.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

basin_area_rows <- vector("list", nrow(scenarios))
for (j in seq_len(nrow(scenarios))) {
  scenario <- scenarios$scenario[j]
  species_stack <- terra::rast(stack_files[[scenario]])
  terra::compareGeom(map.current[[1]], species_stack[[1]], stopOnError = TRUE)
  if (!identical(names(species_stack), fish.info$fcode)) {
    stop("Species-stack layer names/order differ in: ", scenario)
  }
  habitat_area <- species_stack * cell.area
  names(habitat_area) <- fish.info$fcode
  z <- as.data.frame(terra::zonal(habitat_area, basin.rast, fun = "sum", na.rm = TRUE))
  names(z)[1] <- "basin_index"
  if (!all(fish.info$fcode %in% names(z))) stop("Species columns are missing from zonal output.")
  basin_index <- match(basin.lookup$basin_index, z$basin_index)
  rows_this_scenario <- vector("list", nrow(fish.info))
  for (i in seq_len(nrow(fish.info))) {
    area_values <- rep(0, nrow(basin.lookup))
    matched <- !is.na(basin_index)
    area_values[matched] <- z[[fish.info$fcode[i]]][basin_index[matched]]
    # Source convention: missing/all-NA basin sums mean zero suitable area.
    area_values[is.na(area_values)] <- 0
    if (any(!is.finite(area_values)) || any(area_values < 0)) {
      stop("Non-finite or negative suitable-habitat area for: ", fish.info$fcode[i])
    }
    rows_this_scenario[[i]] <- data.frame(
      fcode = fish.info$fcode[i], PC1.centroid = fish.info$PC1.centroid[i],
      PC2.centroid = fish.info$PC2.centroid[i], type = fish.info$thermal.guild[i],
      BASIN_ID = basin.lookup$BASIN_ID, BASIN_NM = basin.lookup$BASIN_NM,
      river_bas0 = basin.lookup$river_bas0, period = scenarios$period[j],
      SSP = scenarios$SSP[j], scenario = scenario, area.km2 = area_values,
      stringsAsFactors = FALSE
    )
  }
  basin_area_rows[[j]] <- do.call(rbind, rows_this_scenario)
}
DB.basin.area <- do.call(rbind, basin_area_rows)
rownames(DB.basin.area) <- NULL
write.csv(DB.basin.area, file.path(table_dir, "species_habitat_area_by_basin.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 8. Calculate basin-level and national suitable-habitat change ------------
DB.current.basin <- DB.basin.area[DB.basin.area$scenario == current_scenario, , drop = FALSE]
DB.basin.diff <- DB.basin.area[DB.basin.area$scenario != current_scenario, , drop = FALSE]
# Use a delimiter that cannot be confused with underscores in basin identifiers.
current_keys <- paste(DB.current.basin$fcode, DB.current.basin$BASIN_ID, sep = "\r")
future_keys <- paste(DB.basin.diff$fcode, DB.basin.diff$BASIN_ID, sep = "\r")
if (anyDuplicated(current_keys)) stop("Duplicate species-basin keys in current area table.")
current_index <- match(future_keys, current_keys)
if (anyNA(current_index)) stop("Some future species-basin rows lack a current reference.")
names(DB.basin.diff)[names(DB.basin.diff) == "area.km2"] <- "area.future.km2"
DB.basin.diff$area.current.km2 <- DB.current.basin$area.km2[current_index]
DB.basin.diff$area.diff.km2 <- DB.basin.diff$area.future.km2 - DB.basin.diff$area.current.km2
DB.basin.diff$Change <- NA_real_
positive_current <- DB.basin.diff$area.current.km2 > 0
DB.basin.diff$Change[positive_current] <-
  DB.basin.diff$area.diff.km2[positive_current] /
  DB.basin.diff$area.current.km2[positive_current] * 100

# Current area = 0 leaves percentage change undefined (NA), including gains.
# Preserve the executable source logic, not its conflicting comment for "new".
DB.basin.diff$change.status <- "persisting"
DB.basin.diff$change.status[
  DB.basin.diff$area.current.km2 == 0 & DB.basin.diff$area.future.km2 == 0
] <- "absent"
DB.basin.diff$change.status[
  DB.basin.diff$area.current.km2 == 0 & DB.basin.diff$area.future.km2 > 0
] <- "new"
DB.basin.diff$change.status[
  DB.basin.diff$area.current.km2 > 0 & DB.basin.diff$area.future.km2 == 0
] <- "lost"
write.csv(DB.basin.diff, file.path(table_dir, "species_habitat_change_by_basin.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# Sum areas BEFORE calculating national percentages; do not average basin changes.
DB.total.diff <- stats::aggregate(
  cbind(area.current.km2, area.future.km2) ~
    fcode + PC1.centroid + PC2.centroid + type + period + SSP + scenario,
  data = DB.basin.diff, FUN = sum, na.rm = TRUE
)
DB.total.diff$area.diff.km2 <- DB.total.diff$area.future.km2 - DB.total.diff$area.current.km2
DB.total.diff$Change <- NA_real_
positive_total <- DB.total.diff$area.current.km2 > 0
DB.total.diff$Change[positive_total] <-
  DB.total.diff$area.diff.km2[positive_total] /
  DB.total.diff$area.current.km2[positive_total] * 100
write.csv(DB.total.diff, file.path(table_dir, "species_habitat_change_national.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# Save spatial tables and paths, not external-pointer SpatRaster objects.
analysis_settings <- list(
  excluded_basin_groups = excluded_basin_groups,
  analysis_longitude_max = analysis_longitude_max,
  constraint_field = basin_group_field, area_group_field = basin_id_field,
  between_guild_normalization = "per_map_spatial_maximum",
  within_guild_normalization = "none", run_far_future_lmm = run_far_future_lmm
)
save(fish.info, scenarios, guild_lookup, species_basin_groups,
     species_basin_assignments, basin_assignment_summary, constrained_manifest,
     stack_files, richness_files, change_files, richness_summary, basin.lookup,
     DB.comp.inter, DB.comp.intra, DB.basin.area, DB.basin.diff, DB.total.diff,
     analysis_settings, file = file.path(output_dir, "distributional_change_overlap.RData"))

# ---- 9. Far-future suitable-habitat change: Gaussian mixed model --------------
# This section reads the saved area-change table, so it can also be executed
# separately after the Section 0 settings/packages have been evaluated.
if (run_far_future_lmm) {
  DB <- read.csv(file.path(table_dir, "species_habitat_change_by_basin.csv"),
                 fileEncoding = output_encoding, stringsAsFactors = FALSE,
                 check.names = FALSE, colClasses = c(BASIN_ID = "character"))
  DB.glmm <- DB[DB$period == "far-future" & DB$area.current.km2 > 0, , drop = FALSE]
  if (!nrow(DB.glmm)) stop("No far-future rows have positive current habitat area.")
  if (anyNA(DB.glmm[, c("fcode", "BASIN_ID", "type", "SSP",
                        "area.current.km2", "area.future.km2")]) ||
      any(!is.finite(DB.glmm$area.current.km2)) ||
      any(!is.finite(DB.glmm$area.future.km2)) || any(DB.glmm$area.future.km2 < 0)) {
    stop("The LMM input contains missing, non-finite, or invalid area values.")
  }
  positive_area <- c(DB.glmm$area.current.km2[DB.glmm$area.current.km2 > 0],
                     DB.glmm$area.future.km2[DB.glmm$area.future.km2 > 0])
  epsilon <- min(positive_area) / 2
  # Apply the source correction ONLY to zero future area, not to both periods.
  DB.glmm$area.future.adj <- ifelse(DB.glmm$area.future.km2 == 0, epsilon,
                                    DB.glmm$area.future.km2)
  DB.glmm$logRR <- log(DB.glmm$area.future.adj / DB.glmm$area.current.km2)
  DB.glmm$type <- factor(DB.glmm$type, levels = guild_levels)
  DB.glmm$SSP <- factor(DB.glmm$SSP, levels = ssp_levels)
  DB.glmm$fcode <- factor(DB.glmm$fcode)
  DB.glmm$BASIN_ID <- factor(DB.glmm$BASIN_ID)
  if (anyNA(DB.glmm$type) || anyNA(DB.glmm$SSP) ||
      any(!is.finite(DB.glmm$logRR)) || any(table(DB.glmm$type, DB.glmm$SSP) == 0)) {
    stop("The LMM needs finite logRR values and all guild-by-SSP combinations.")
  }
  lmm_counts <- as.data.frame(table(DB.glmm$type, DB.glmm$SSP), stringsAsFactors = FALSE)
  names(lmm_counts) <- c("type", "SSP", "n_observations")
  write.csv(DB.glmm, file.path(lmm_dir, "far_future_lmm_input.csv"),
            row.names = FALSE, fileEncoding = output_encoding)
  write.csv(lmm_counts, file.path(lmm_dir, "far_future_lmm_sample_counts.csv"),
            row.names = FALSE, fileEncoding = output_encoding)
  write.csv(data.frame(epsilon_km2 = epsilon,
                       n_zero_future_rows = sum(DB.glmm$area.future.km2 == 0)),
            file.path(lmm_dir, "log_response_ratio_correction.csv"),
            row.names = FALSE, fileEncoding = output_encoding)
  
  # Equivalent sum contrasts to the source, specified locally to avoid changing global R options. 
  # Thermal PCA centroid columns are NOT fixed-effect terms.
  model.full <- lmerTest::lmer(
    logRR ~ type * SSP + (1 | fcode) + (1 | BASIN_ID) + (1 | fcode:BASIN_ID),
    data = DB.glmm, REML = TRUE,
    contrasts = list(type = "contr.sum", SSP = "contr.sum"),
    control = lme4::lmerControl(optimizer = "bobyqa"), na.action = na.fail
  )
  lmm_singular <- lme4::isSingular(model.full, tol = singularity_tolerance)
  lmm_variance <- as.data.frame(lme4::VarCorr(model.full))
  fixed_effects <- as.data.frame(coef(summary(model.full)))
  fixed_effects <- cbind(term = rownames(fixed_effects), fixed_effects)
  rownames(fixed_effects) <- NULL
  writeLines(capture.output(summary(model.full)), file.path(lmm_dir, "lmm_summary.txt"))
  write.csv(fixed_effects, file.path(lmm_dir, "lmm_fixed_effects.csv"),
            row.names = FALSE, fileEncoding = output_encoding)
  write.csv(lmm_variance, file.path(lmm_dir, "lmm_variance_components.csv"),
            row.names = FALSE, fileEncoding = output_encoding)
  writeLines(capture.output(list(
    singular = lmm_singular, tolerance = singularity_tolerance,
    convergence = model.full@optinfo$conv
  )), file.path(lmm_dir, "lmm_diagnostics.txt"))
  save(model.full, DB.glmm, epsilon, lmm_singular,
       file = file.path(lmm_dir, "far_future_lmm_model.RData"))
  if (lmm_singular) {
    warning("The source LMM is singular. Its random-effects structure has NOT been simplified.")
  }
  
  # H1: thermal-guild effect; H2: guild-by-SSP interaction, as in the source.
  Type3test <- stats::anova(model.full, type = 3, ddf = "Kenward-Roger")
  type3_table <- data.frame(term = rownames(Type3test), as.data.frame(Type3test),
                            row.names = NULL, check.names = FALSE)
  write.csv(type3_table, file.path(lmm_dir, "lmm_type_III_tests.csv"),
            row.names = FALSE, fileEncoding = output_encoding)
  
  # Keep the source's emmeans defaults. These do not guarantee the same df
  # method as the explicitly Kenward-Roger Type III test, especially at large N.
  emm.type <- emmeans::emmeans(model.full, ~ type | SSP)
  emm.pairs <- emmeans::contrast(emm.type, method = "pairwise", adjust = "tukey")
  emm_summary <- summary(emm.type)
  pair_summary <- summary(emm.pairs)
  write.csv(as.data.frame(emm_summary), file.path(lmm_dir, "lmm_guild_means_by_SSP.csv"),
            row.names = FALSE, fileEncoding = output_encoding)
  write.csv(as.data.frame(pair_summary), file.path(lmm_dir, "lmm_guild_pairs_by_SSP.csv"),
            row.names = FALSE, fileEncoding = output_encoding)
  writeLines(capture.output(list(
    marginal_means = emm_summary, pairwise_tests = pair_summary,
    emmeans_options = emmeans::emm_options()
  )), file.path(lmm_dir, "lmm_emmeans_details.txt"))
  save(model.full, DB.glmm, epsilon, Type3test, emm.type, emm.pairs,
       lmm_singular, lmm_variance,
       file = file.path(lmm_dir, "far_future_lmm_results.RData"))
}

# ---- 10. Record the software environment ------------------------------------
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
package_versions <- data.frame(
  package = required_packages,
  version = vapply(required_packages, function(x) as.character(utils::packageVersion(x)),
                   character(1))
)
write.csv(package_versions, file.path(output_dir, "package_versions.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
message("Distributional change and overlap analyses complete. Outputs: ", output_dir)
