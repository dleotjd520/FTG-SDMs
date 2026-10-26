# =============================================================================
# Data collection and preprocessing for freshwater fish distribution analyses
# =============================================================================
# Execute the numbered sections sequentially in an interactive R session.
# Settings are defined below; there is no separate configuration file or main().
# Two small sampling helpers are defined next to the steps that use them.

# ---- 0. Paths, settings, and packages -----------------------------------------
project_dir <- getwd()  # Start R in the directory containing this script.
raw_dir <- file.path(project_dir, "data", "raw")
output_dir <- file.path(project_dir, "outputs", "data_preprocessing")
gbif_dir <- file.path(output_dir, "gbif_download")
background_dir <- file.path(output_dir, "background_pool")
species_dir <- file.path(output_dir, "species_datasets")
model_input_dir <- file.path(output_dir, "model_inputs")

monitoring_records_file <- file.path(raw_dir, "biomonitoring_records.csv")
monitoring_sites_file <- file.path(raw_dir, "biomonitoring_sites.csv")
species_information_file <- file.path(raw_dir, "fish_species_information.csv")
climate_dir <- file.path(raw_dir, "CHELSA_1981_2010")
elevation_file <- file.path(raw_dir, "wc2.1_30s_elev.tif")
boundary_file <- file.path(raw_dir, "korea_boundary_buffer.shp")
coastline_file <- file.path(raw_dir, "ne_10m_coastline.shp")

# USER ACTION REQUIRED FOR A NEW GBIF DOWNLOAD ONLY:
# Enter YOUR OWN GBIF username, password, and registered email
# in your local R console BEFORE running this section.
# Sys.setenv(
#   GBIF_USER  = "YOUR_GBIF_USERNAME",
#   GBIF_PWD   = "YOUR_GBIF_PASSWORD",
#   GBIF_EMAIL = "YOUR_GBIF_ACCOUNT_EMAIL"
# )
input_encoding <- "EUC-KR"  # Change to "UTF-8" for re-encoded input tables.
output_encoding <- "EUC-KR"
request_new_gbif_download <- FALSE
gbif_download_key <- "0004225-250811113504898" # Historical download key 
generate_background_pool <- TRUE
background_pool_file <- file.path(background_dir, "fish_background.csv")

initial_gbif_year <- 1950L
analysis_start_year <- 1981L
coordinate_precision_limit <- 0.01
coordinate_uncertainty_limit_m <- 10000
excluded_uncertainty_values <- c(301, 3036, 999, 9999)
coordinate_cleaning_buffer_m <- 2000
allowed_record_types <- c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "OBSERVATION", "OCCURRENCE")
background_pool_target <- 20000L
background_cores <- 5L
background_latitude_limits <- c(-60, 80)
environmental_bins <- 25L
pca_variance_cutoff <- 0.95
minimum_distance_m <- 5000
maximum_background_per_species <- 10000L
minimum_presence_records <- 10L
thinning_batch_size <- 100L
coastal_cluster_radius_m <- 1000
coastal_distance_cutoff_m <- 5000
correlation_threshold <- 0.9
vif_threshold <- 5
partition_repeats <- 10L
test_fraction <- 0.2
background_seed <- 910520L  # Explicit seed added; archive the resulting pool.
thinning_seed <- 910520L
partition_seed <- 910520L
variable_selection_seed <- 910520L  # Explicit seed added for usdm sampling.

required_packages <- c("rgbif", "dplyr", "CoordinateCleaner", "terra", "sf",
                       "dbscan", "geodist", "usdm")
if (generate_background_pool) required_packages <- c(required_packages, "megaSDM")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install the required R package before continuing: ", package)
  }
}
library(dplyr)
sf::sf_use_s2(TRUE)  # Explicitly use spherical distances in the sf steps.

for (directory in c(output_dir, gbif_dir, background_dir,
                    species_dir, model_input_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
required_files <- c(monitoring_records_file, monitoring_sites_file,
                    species_information_file, elevation_file,
                    boundary_file, coastline_file)
if (any(!file.exists(required_files))) {
  stop("Missing input files: ",
       paste(required_files[!file.exists(required_files)], collapse = "; "))
}
if (!dir.exists(climate_dir)) stop("Climate directory does not exist: ", climate_dir)

# ---- 1. Import national monitoring records and the species metadata -----------
BM.sp.long <- read.csv(monitoring_records_file, fileEncoding = input_encoding,
                       stringsAsFactors = FALSE, check.names = FALSE)
BM.site <- read.csv(monitoring_sites_file, fileEncoding = input_encoding,
                    stringsAsFactors = FALSE, check.names = FALSE)
DB.sp.info <- read.csv(species_information_file, fileEncoding = input_encoding,
                       stringsAsFactors = FALSE, check.names = FALSE)

legacy_headers <- c(
  "Stream.code", "LNG", "LAT", "\uad6d\uba85", "\ud559\uba85",
  "Abundance.\uc5f0\ud3c9\uade0.", "Abundance(\uc5f0\ud3c9\uade0)",
  "SDM\ub300\uc0c1\uc5ec\ubd80", "Fresh", "Brack", "Saltwater"
)
english_headers <- c(
  "site_id", "lon", "lat", "common_name", "scientific_name",
  "abundance", "abundance", "freshwater", "brackish", "saltwater"
)
input_tables <- list(BM.sp.long, BM.site, DB.sp.info)
for (j in seq_along(input_tables)) {
  header_match <- match(names(input_tables[[j]]), legacy_headers)
  rename <- !is.na(header_match)
  names(input_tables[[j]])[rename] <- english_headers[header_match[rename]]
  if (anyDuplicated(names(input_tables[[j]]))) {
    stop("Duplicate column names after header conversion in input table ", j)
  }
}
BM.sp.long <- input_tables[[1]]
BM.site <- input_tables[[2]]
DB.sp.info <- input_tables[[3]]
rm(input_tables)

required_headers <- list(
  c("site_id", "common_name", "scientific_name", "abundance"),
  c("site_id", "lon", "lat"),
  c("fcode", "common_name", "scientific_name")
)
for (j in seq_along(required_headers)) {
  table_names <- names(list(BM.sp.long, BM.site, DB.sp.info)[[j]])
  missing_headers <- setdiff(required_headers[[j]], table_names)
  if (length(missing_headers)) {
    stop("Missing headers in input table ", j, ": ",
         paste(missing_headers, collapse = ", "))
  }
}
if (!is.numeric(BM.sp.long$abundance) ||
    !is.numeric(BM.site$lon) || !is.numeric(BM.site$lat)) {
  stop("Abundance, longitude, and latitude must be numeric.")
}
if (anyNA(BM.site$site_id) || anyDuplicated(BM.site$site_id)) {
  stop("The site table must have one non-missing row per site_id.")
}
if (anyNA(DB.sp.info[, required_headers[[3]]]) ||
    anyDuplicated(DB.sp.info$common_name) || anyDuplicated(DB.sp.info$fcode)
    ) {
  stop("Species metadata must have unique names/codes.")
}

# Keep only fields needed here; unrelated traits from the metadata are excluded.
metadata_fields <- c(required_headers[[3]],
                     intersect(c("freshwater", "brackish", "saltwater"),
                               names(DB.sp.info)))
DB.sp.info <- DB.sp.info[, metadata_fields, drop = FALSE]
BM.sp.long2 <- BM.sp.long[!is.na(BM.sp.long$abundance) &
                            BM.sp.long$abundance > 0, , drop = FALSE]
if (!nrow(BM.sp.long2) ||
    anyNA(BM.sp.long2[, c("site_id", "common_name", "scientific_name")]) ||
    any(!nzchar(trimws(BM.sp.long2$common_name))) ||
    any(!nzchar(trimws(BM.sp.long2$scientific_name)))) {
  stop("Positive monitoring records require non-missing site and taxon names.")
}
sp.list <- unique(BM.sp.long2[, c("common_name", "scientific_name")])
if (anyDuplicated(sp.list$common_name)) {
  stop("A monitoring common name maps to multiple scientific names; reconcile first.")
}
DB.BM <- merge(BM.sp.long2[, c("site_id", "common_name")],
               BM.site[, c("site_id", "lon", "lat")], by = "site_id", all.x = TRUE)

# ---- 2. Retrieve a GBIF occurrence download -----------------------------------
f (request_new_gbif_download) {
  query_names <- sort(unique(BM.sp.long$scientific_name))
  query_names <- query_names[!is.na(query_names) & nzchar(trimws(query_names))]
  taxon_lookup <- data.frame(query_name = query_names, usage_key = NA_real_,
                             matched_name = NA_character_,
                             match_type = NA_character_)
  for (i in seq_along(query_names)) {
    backbone <- rgbif::name_backbone(name = query_names[i])
    if (length(backbone$usageKey) != 1L || is.na(backbone$usageKey)) {
      stop("Unresolved GBIF name: ", query_names[i])
    }
    taxon_lookup$usage_key[i] <- backbone$usageKey
    if (length(backbone$scientificName) == 1L) {
      taxon_lookup$matched_name[i] <- backbone$scientificName
    }
    if (length(backbone$matchType) == 1L) {
      taxon_lookup$match_type[i] <- backbone$matchType
    }
  }
  write.csv(taxon_lookup, file.path(gbif_dir, "gbif_taxon_lookup.csv"),
            row.names = FALSE, fileEncoding = output_encoding)
  if (any(!is.na(taxon_lookup$match_type) & taxon_lookup$match_type != "EXACT")) {
    stop("Review non-exact matches in gbif_taxon_lookup.csv before requesting data.")
  }
  if (any(!nzchar(Sys.getenv(c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL"))))) {
    stop("Set GBIF_USER, GBIF_PWD, and GBIF_EMAIL privately for a new download.")
  }
  gbif_download <- rgbif::occ_download(
    rgbif::pred_in("taxonKey", unique(taxon_lookup$usage_key)),
    rgbif::pred("hasCoordinate", TRUE),
    rgbif::pred("hasGeospatialIssue", FALSE),
    format = "SIMPLE_CSV", 
    user = Sys.getenv("GBIF_USER"), # 
    pwd = Sys.getenv("GBIF_PWD"),   # 
    email = Sys.getenv("GBIF_EMAIL")# 
  )
  gbif_download_key <- as.character(gbif_download)
  writeLines(gbif_download_key, file.path(gbif_dir, "gbif_download_key.txt"))
  rgbif::occ_download_wait(gbif_download)
}

if (length(gbif_download_key) != 1L || !nzchar(gbif_download_key)) {
  stop("Provide one GBIF download key.")
}
writeLines(gbif_download_key, file.path(gbif_dir, "gbif_download_key.txt"))
# Import an already downloaded archive without downloading it again.
archive_file <- file.path(gbif_dir, paste0(gbif_download_key, ".zip"))
if (!file.exists(archive_file)) {
  rgbif::occ_download_get(gbif_download_key, path = gbif_dir, overwrite = FALSE)
}
data <- as.data.frame(rgbif::occ_download_import(
  key = gbif_download_key, path = gbif_dir
))
names(data) <- tolower(names(data))
gbif_metadata <- rgbif::occ_download_meta(gbif_download_key)
saveRDS(gbif_metadata, file.path(gbif_dir, "gbif_download_metadata.rds"))
writeLines(capture.output(print(rgbif::gbif_citation(gbif_metadata))),
           file.path(gbif_dir, "gbif_citation.txt"))

# ---- 3. Clean GBIF coordinates and restrict records ---------------------------

required_gbif_fields <- c(
  "gbifid", "occurrencestatus", "basisofrecord", "year",
  "coordinateprecision", "coordinateuncertaintyinmeters",
  "decimallongitude", "decimallatitude", "specieskey", "datasetkey",
  "species", "scientificname"
)
if (length(setdiff(required_gbif_fields, names(data)))) {
  stop("GBIF archive is missing fields required by the original filtering logic: ",
       paste(setdiff(required_gbif_fields, names(data)), collapse = ", "))
}

# Filter occurrence records and coordinate quality
data$gbifid <- as.character(data$gbifid)
gbif_filter_log <- data.frame(stage = "downloaded", n_records = nrow(data))

data <- data %>%
  filter(occurrencestatus == "PRESENT",
         !basisofrecord %in% c("PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN", "LIVING_SPECIMEN"),
         year >= initial_gbif_year,
         coordinateprecision < coordinate_precision_limit | is.na(coordinateprecision),
         coordinateuncertaintyinmeters < coordinate_uncertainty_limit_m |
           is.na(coordinateuncertaintyinmeters),
         !coordinateuncertaintyinmeters %in% excluded_uncertainty_values,
         is.finite(decimallongitude), is.finite(decimallatitude),
         abs(decimallongitude) <= 180, abs(decimallatitude) <= 90,
         decimallongitude != 0 | decimallatitude != 0)

# Remove spatially suspect and duplicate records
data <- data %>%
  CoordinateCleaner::cc_cen(
    buffer = coordinate_cleaning_buffer_m,
    lon = "decimallongitude",
    lat = "decimallatitude"
  ) %>%
  CoordinateCleaner::cc_cap(
    buffer = coordinate_cleaning_buffer_m,
    lon = "decimallongitude",
    lat = "decimallatitude"
  ) %>%
  CoordinateCleaner::cc_inst(
    buffer = coordinate_cleaning_buffer_m,
    lon = "decimallongitude",
    lat = "decimallatitude"
  ) %>%
  CoordinateCleaner::cc_sea(
    lon = "decimallongitude",
    lat = "decimallatitude"
  ) %>%
  dplyr::distinct(
    decimallongitude, decimallatitude, specieskey, datasetkey,
    .keep_all = TRUE
  ) %>%
  dplyr::mutate(
    gbifid2 = paste0("GBIF", dplyr::row_number())
  ) %>%
  as.data.frame()

gbif_filter_log <- rbind(gbif_filter_log,
                         data.frame(stage = "coordinate_cleaner_and_deduplication", n_records = nrow(data)))
save(data, file = file.path(output_dir, "gbif_occurrences_cleaned.RData"))

tGBIF2 <- data[data$year >= analysis_start_year &
                 data$basisofrecord %in% allowed_record_types, , drop = FALSE]
gbif_filter_log <- rbind(gbif_filter_log,
                         data.frame(stage = "analysis_year_and_record_types", n_records = nrow(tGBIF2)))
write.csv(gbif_filter_log, file.path(output_dir, "gbif_filter_counts.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 4. Match GBIF records to monitoring taxa and combine occurrences ----------
sp.list$gbif_match_name <- vapply(strsplit(sp.list$scientific_name, " +"),
                                  function(x) paste(head(x, 2L), collapse = " "),
                                  character(1))
carp_form <- sp.list$common_name == "\uc774\uc2a4\ub77c\uc5d8\uc789\uc5b4"
sp.list$scientific_name[carp_form] <- "Cyprinus carpio nudus"
sp.list$gbif_match_name[carp_form] <- "carpio nudus"

# Species codes and species-loop order follow the original scientific-name sort.
sp.list2 <- sp.list[order(sp.list$scientific_name), , drop = FALSE]
sp.list2$fcode <- paste0("F", sprintf("%03d", seq_len(nrow(sp.list2))))
row.names(sp.list2) <- NULL
metadata_match <- match(sp.list2$common_name, DB.sp.info$common_name)
if (anyNA(metadata_match)) {
  stop("Species metadata are missing monitoring taxa: ",
       paste(sp.list2$common_name[is.na(metadata_match)], collapse = ", "))
}
if (!all(sp.list2$fcode == DB.sp.info$fcode[metadata_match])) {
  stop("Generated species codes disagree with the supplied metadata. Review the ",
       "taxon list and sort order; do not silently renumber archived species.")
}
write.csv(sp.list2, file.path(output_dir, "species_registry.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

match_parts <- list()
remaining <- seq_len(nrow(tGBIF2))
for (match_field in c("species", "scientificname")) {
  matched_in_pass <- integer()
  for (i in seq_len(nrow(sp.list))) {
    hit <- remaining[grep(sp.list$gbif_match_name[i],
                          as.character(tGBIF2[[match_field]][remaining]),
                          ignore.case = TRUE)]
    if (!length(hit)) next
    match_parts[[length(match_parts) + 1L]] <- data.frame(
      gbif_row = hit, gbifid2 = tGBIF2$gbifid2[hit], gbifid = tGBIF2$gbifid[hit],
      common_name = sp.list$common_name[i], match_field = match_field,
      gbif_species = tGBIF2$species[hit], gbif_scientific_name = tGBIF2$scientificname[hit]
    )
    matched_in_pass <- c(matched_in_pass, hit)
  }
  remaining <- setdiff(remaining, unique(matched_in_pass))
}
if (!length(match_parts)) stop("No GBIF records matched the monitoring taxa.")
gbif_match_audit <- do.call(rbind, match_parts)
write.csv(gbif_match_audit, file.path(output_dir, "gbif_name_matching.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
if (anyDuplicated(gbif_match_audit$gbif_row)) {
  stop("A GBIF record matched multiple taxa. Review gbif_name_matching.csv; ",
       "the script does not resolve ambiguous taxon assignments automatically.")
}
write.csv(tGBIF2[remaining, c("gbifid", "species", "scientificname"), drop = FALSE],
          file.path(output_dir, "gbif_unmatched_records.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

DB.BM2 <- DB.BM[, c("site_id", "common_name", "lon", "lat")]
names(DB.BM2) <- c("siteID", "species", "lon", "lat")
DB.GBIF2 <- data.frame(
  siteID = gbif_match_audit$gbifid2, species = gbif_match_audit$common_name,
  lon = tGBIF2$decimallongitude[gbif_match_audit$gbif_row],
  lat = tGBIF2$decimallatitude[gbif_match_audit$gbif_row]
)
DB.fish <- rbind(DB.BM2, DB.GBIF2)
valid_coordinates <- is.finite(DB.fish$lon) & is.finite(DB.fish$lat) &
  abs(DB.fish$lon) <= 180 & abs(DB.fish$lat) <= 90
write.csv(DB.fish[!valid_coordinates, , drop = FALSE],
          file.path(output_dir, "occurrences_invalid_coordinates.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
DB.fish <- DB.fish[valid_coordinates, , drop = FALSE]
row.names(DB.fish) <- NULL
# Repeated monitoring observations are intentionally not deduplicated here.
save(DB.fish, sp.list2, file = file.path(output_dir, "combined_occurrences.RData"))

# ---- 5. Load baseline environmental rasters and extract values ----------------
climate_files <- list.files(climate_dir, pattern = "\\.(tif|tiff)$",
                            full.names = TRUE, ignore.case = TRUE)
climate_files <- climate_files[grepl("bio[0-9]+", basename(climate_files), ignore.case = TRUE)]
climate_names <- sub(".*?(bio[0-9]+).*", "\\1", tolower(basename(climate_files)), perl = TRUE)
if (length(climate_files) != 19L || anyDuplicated(climate_names) ||
    !setequal(climate_names, paste0("bio", 1:19))) {
  stop("The climate directory must contain one baseline GeoTIFF for each of bio1-bio19.")
}
# Keep the original list.files() ordering; predictor order can affect selection.
Env.CHELSA <- terra::rast(climate_files)
names(Env.CHELSA) <- climate_names
Env.ele <- terra::rast(elevation_file)
if (terra::nlyr(Env.ele) != 1L ||
    !terra::same.crs(Env.CHELSA, "EPSG:4326") ||
    !terra::same.crs(Env.ele, Env.CHELSA)) {
  stop("Provide a single elevation layer and WGS84 rasters with matching CRSs.")
}
Env.ele.c <- terra::crop(Env.ele, Env.CHELSA)
terra::compareGeom(Env.ele.c, Env.CHELSA[[1]], stopOnError = TRUE)
names(Env.ele.c) <- "elev"
Env.maps <- c(Env.ele.c, Env.CHELSA)
env_names <- names(Env.maps)

# Raster values are used as supplied. No undocumented unit conversion is added.
write.csv(data.frame(variable = env_names,
                     source_file = basename(c(elevation_file, climate_files))),
          file.path(output_dir, "environmental_layers.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

occurrence_points <- terra::vect(DB.fish, geom = c("lon", "lat"), crs = "EPSG:4326")
DB.env <- terra::extract(Env.maps, occurrence_points, ID = FALSE)
DB.env <- DB.env[, env_names, drop = FALSE]
DB.fish.env <- cbind(DB.fish, DB.env)
complete_environment <- complete.cases(DB.env) &
  rowSums(!is.finite(as.matrix(DB.env))) == 0L
DB.fish.env2 <- DB.fish.env[complete_environment, , drop = FALSE]
write.csv(DB.fish.env[!complete_environment, c("siteID", "species", "lon", "lat")],
          file.path(output_dir, "occurrences_missing_environment.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
save(DB.fish.env2, sp.list2, env_names,
     file = file.path(output_dir, "occurrence_environment_complete.RData"))

# ---- 6. Generate or reuse the shared global background pool -------------------
# spatial_weights = 0 preserves global sampling without species-specific buffers.
# The requested pool size is approximate, not an enforced final count.
if (generate_background_pool) {
  if (file.exists(background_pool_file)) {
    stop("A background pool already exists. Set generate_background_pool = FALSE ",
         "to reuse it, or use a new output directory for a new run.")
  }
  set.seed(background_seed)
  megaSDM::BackgroundPoints(
    spplist = "fish", envdata = Env.maps, output = background_dir,
    spatial_weights = 0, nbg = background_pool_target,
    method = "Varela", ncores = background_cores
  )
}
if (!file.exists(background_pool_file)) stop("Background pool file not found.")
fish.bg <- read.csv(background_pool_file, stringsAsFactors = FALSE, check.names = FALSE)
# Ignore Species or a CSV row-index column by selecting named fields, never [-1].
if (!all(c("x", "y", env_names) %in% names(fish.bg))) {
  stop("Background pool requires x, y, and all 20 named environmental predictors.")
}
fish.bg <- fish.bg[, c("x", "y", env_names), drop = FALSE]
names(fish.bg)[1:2] <- c("lon", "lat")
if (!all(vapply(fish.bg, is.numeric, logical(1)))) {
  stop("Background coordinates and environmental predictors must be numeric.")
}
fish.bg <- fish.bg[complete.cases(fish.bg) &
                     rowSums(!is.finite(as.matrix(fish.bg))) == 0L, , drop = FALSE]
fish.bg <- fish.bg[fish.bg$lat >= background_latitude_limits[1] &
                     fish.bg$lat <= background_latitude_limits[2] &
                     abs(fish.bg$lon) <= 180, , drop = FALSE]
if (!nrow(fish.bg)) stop("No valid background points remain.")

# ---- 7. Two sampling helpers retained from the original workflow --------------
# PCA environmental-bin sampling, adapted from the supplied VarelaSample helper
# (which cites megaSDM/R/OccurrenceManagement.R by Shipley and collaborators).
# Return source-row indices, not joined coordinate tables. One random record is
# retained per occupied multidimensional bin; selected indices are sorted as in
# the original merge-based return. Unused non-PCA/multiple-bin modes are removed.
environmental_filter_rows <- function(environment, bins = 25L, coverage = 0.95) {
  if (!is.data.frame(environment) || !nrow(environment)) return(integer())
  if (anyNA(environment) || any(!is.finite(as.matrix(environment)))) {
    stop("Environmental filtering requires finite, complete predictor values.")
  }
  if (nrow(environment) == 1L) return(1L)
  variance <- vapply(environment, stats::var, numeric(1))
  environment <- environment[, is.finite(variance) & variance > 0, drop = FALSE]
  # Degenerate case: no environmental distinction is possible.
  if (!ncol(environment)) return(sample.int(nrow(environment), 1L))
  pca <- stats::prcomp(environment, scale. = TRUE)
  cumulative <- summary(pca)$importance[3, ]
  axes <- min(ncol(pca$x), max(2L, which(cumulative > coverage)[1]))
  scores <- as.data.frame(pca$x[, seq_len(axes), drop = FALSE])
  bin_membership <- scores
  for (j in seq_len(ncol(scores))) {
    limits <- range(scores[[j]])
    if (diff(limits) == 0) {
      bin_membership[[j]] <- 1L
    } else {
      bin_membership[[j]] <- ceiling((scores[[j]] - limits[1]) / (diff(limits) / bins))
      bin_membership[[j]][bin_membership[[j]] == 0] <- 1L
    }
  }
  groups <- dplyr::distinct(bin_membership)
  groups$group_id <- seq_len(nrow(groups))
  assigned <- dplyr::left_join(bin_membership, groups, by = names(bin_membership))
  selected <- integer(nrow(groups))
  for (j in seq_len(nrow(groups))) {
    candidates <- which(assigned$group_id == j)
    selected[j] <- candidates[sample.int(length(candidates), 1L)]
  }
  sort(selected)
}

# Sequential geographic thinning, simplified from thin_points_batch_terra.
# WGS84 coordinates only. Every retained pair is strictly farther apart than threshold_m. 
# Checks are applied both across batches and within each batch.
thin_background_points <- function(points, threshold_m, batch_size = 100L,
                                   seed = 910520L) {
  if (!nrow(points)) return(points)
  stopifnot(all(c("lon", "lat") %in% names(points)),
            is.finite(threshold_m), threshold_m >= 0, batch_size >= 1L)
  set.seed(seed)
  points <- points[sample.int(nrow(points)), , drop = FALSE]
  coordinates <- terra::vect(points[, c("lon", "lat")],
                             geom = c("lon", "lat"), crs = "EPSG:4326")
  retained <- integer()
  for (start in seq.int(1L, nrow(points), by = batch_size)) {
    batch <- seq.int(start, min(start + batch_size - 1L, nrow(points)))
    batch_points <- coordinates[batch]
    eligible <- rep(TRUE, length(batch))
    if (length(retained)) {
      distances <- matrix(terra::distance(batch_points, coordinates[retained]),
                          nrow = length(batch))
      eligible <- apply(distances, 1L, min) > threshold_m
    }
    local_retained <- integer()
    if (any(eligible)) {
      within <- matrix(terra::distance(batch_points, batch_points),
                       nrow = length(batch))
      for (j in which(eligible)) {
        if (length(local_retained) &&
            any(within[j, local_retained] <= threshold_m)) next
        local_retained <- c(local_retained, j)
      }
    }
    retained <- c(retained, batch[local_retained])
  }
  result <- points[retained, , drop = FALSE]
  row.names(result) <- NULL
  result
}

# ---- 8. Build species-specific presence-background datasets -------------------
sp.info <- sp.list2
sp.info$N.raw <- as.integer(table(factor(DB.fish$species, levels = sp.info$common_name)))
sp.info$N.complete <- as.integer(table(factor(DB.fish.env2$species, levels = sp.info$common_name)))
sp.info$N.P <- sp.info$N.A <- 0L
sp.info$processing_status <- "fewer_than_two_complete_occurrences"

for (i in seq_len(nrow(sp.list2))) {
  message("Preparing ", sp.list2$fcode[i], ": ", sp.list2$scientific_name[i])
  tP <- DB.fish.env2[DB.fish.env2$species == sp.list2$common_name[i], , drop = FALSE]
  final.DB.P <- tP[FALSE, , drop = FALSE]
  final.DB.A <- fish.bg[FALSE, , drop = FALSE]
  # Preserve the original handling of taxa with <= 1 complete occurrence.
  if (nrow(tP) > 1L) {
    other_presences <- DB.fish.env2[
      DB.fish.env2$species != sp.list2$common_name[i],
      c("lon", "lat", env_names), drop = FALSE
    ]
    # These are candidate BACKGROUND locations, not confirmed absences.
    temp.A <- rbind(fish.bg, other_presences)
    set.seed(6808L + i)
    sampled_A <- environmental_filter_rows(temp.A[, env_names, drop = FALSE],
                                           bins = environmental_bins,
                                           coverage = pca_variance_cutoff)
    tA <- temp.A[sampled_A, , drop = FALSE]
    # Compare backgrounds against ALL complete presences, before presence filtering.
    presence_sf <- sf::st_as_sf(tP, coords = c("lon", "lat"), crs = 4326)
    background_sf <- sf::st_as_sf(tA, coords = c("lon", "lat"), crs = 4326)
    nearby <- sf::st_is_within_distance(background_sf, presence_sf,
                                        dist = minimum_distance_m)
    tA1 <- tA[lengths(nearby) == 0L, , drop = FALSE]
    tA2 <- thin_background_points(tA1, threshold_m = minimum_distance_m,
                                  batch_size = thinning_batch_size, seed = thinning_seed)
    if (nrow(tA2) > maximum_background_per_species) {
      set.seed(123L + i)
      tA2 <- tA2[sample.int(nrow(tA2), maximum_background_per_species), , drop = FALSE]
    }
    final.DB.A <- tA2
    set.seed(919L + i)
    sampled_P <- environmental_filter_rows(tP[, env_names, drop = FALSE],
                                           bins = environmental_bins,
                                           coverage = pca_variance_cutoff)
    final.DB.P <- tP[sampled_P, , drop = FALSE]
    sp.info$processing_status[i] <- "processed"
  }
  sp.info$N.P[i] <- nrow(final.DB.P)
  sp.info$N.A[i] <- nrow(final.DB.A)
  info.sp <- sp.info[i, , drop = FALSE]
  save(info.sp, final.DB.P, final.DB.A,
       file = file.path(species_dir, paste0("presence_background_", info.sp$fcode, ".RData")))
}
write.csv(sp.info, file.path(output_dir, "species_dataset_summary.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 9. Screen domestic distributions using distance to the coastline ---------
# script does not reconstruct FishBase habitat information or manual decisions.
metadata_match <- match(sp.info$common_name, DB.sp.info$common_name)
habitat_fields <- intersect(c("freshwater", "brackish", "saltwater"), names(DB.sp.info))
for (field in habitat_fields) sp.info[[field]] <- DB.sp.info[[field]][metadata_match]
coastal_candidates <- sp.info[sp.info$N.P >= minimum_presence_records, , drop = FALSE]
if (!nrow(coastal_candidates)) stop("No taxa qualify for coastal screening.")
korea_boundary <- terra::vect(boundary_file)
if (!nzchar(terra::crs(korea_boundary))) stop("The boundary file requires a CRS.")
if (!terra::same.crs(korea_boundary, "EPSG:4326")) {
  korea_boundary <- terra::project(korea_boundary, "EPSG:4326")
}
domestic_input <- DB.fish[DB.fish$species %in% coastal_candidates$common_name, , drop = FALSE]
domestic_points <- terra::vect(domestic_input, geom = c("lon", "lat"), crs = "EPSG:4326")
domestic_points <- terra::intersect(korea_boundary, domestic_points)
coastline <- sf::st_read(coastline_file, quiet = TRUE)
if (is.na(sf::st_crs(coastline))) stop("The coastline file requires a CRS.")
coastline <- sf::st_transform(coastline, 4326)
coastline_combined <- sf::st_combine(sf::st_geometry(coastline))

coast_summary <- coastal_candidates[, c("fcode", "common_name", "scientific_name",
                                        habitat_fields), drop = FALSE]
coast_summary$n_domestic_records <- coast_summary$n_clusters <- 0L
coast_summary$distance_mean_m <- coast_summary$distance_median_m <- NA_real_
coast_summary$distance_q70_m <- NA_real_
for (i in seq_len(nrow(coast_summary))) {
  coordinates <- terra::crds(domestic_points[
    domestic_points$species == coast_summary$common_name[i]], df = TRUE)
  coast_summary$n_domestic_records[i] <- nrow(coordinates)
  if (!nrow(coordinates)) next
  if (nrow(coordinates) == 1L) {
    cluster_id <- 1L
  } else {
    distances <- stats::as.dist(geodist::geodist(coordinates, measure = "haversine"))
    cluster_id <- dbscan::dbscan(distances, eps = coastal_cluster_radius_m,
                                 minPts = 1L)$cluster
  }
  # minPts = 1 keeps singleton clusters. Centroids are arithmetic lon/lat means,
  # and each cluster receives equal weight in the coastal-distance summaries.
  centroids <- aggregate(coordinates, by = list(cluster = cluster_id), FUN = mean)
  centroids_sf <- sf::st_as_sf(centroids, coords = c("x", "y"), crs = 4326)
  distance_m <- as.numeric(sf::st_distance(centroids_sf, coastline_combined))
  coast_summary$n_clusters[i] <- length(distance_m)
  coast_summary$distance_mean_m[i] <- mean(distance_m)
  coast_summary$distance_median_m[i] <- median(distance_m)
  coast_summary$distance_q70_m[i] <- unname(stats::quantile(distance_m, probs = 0.70))
}
write.csv(coast_summary, file.path(output_dir, "coastal_screening.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

species_screening <- sp.info
coast_match <- match(species_screening$fcode, coast_summary$fcode)
species_screening$distance_q70_m <- coast_summary$distance_q70_m[coast_match]
species_screening$retained <- species_screening$N.P >= minimum_presence_records &
  !is.na(species_screening$distance_q70_m) &
  species_screening$distance_q70_m >= coastal_distance_cutoff_m
species_screening$screening_result <- "retained"
species_screening$screening_result[species_screening$N.P < minimum_presence_records] <- "insufficient_processed_presences"
species_screening$screening_result[species_screening$N.P >= minimum_presence_records &
                                     is.na(species_screening$distance_q70_m)] <- "no_domestic_coastal_distance"
species_screening$screening_result[!is.na(species_screening$distance_q70_m) &
                                     species_screening$distance_q70_m < coastal_distance_cutoff_m] <- "predominantly_coastal"
write.csv(species_screening, file.path(output_dir, "species_screening.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
freshwater_species <- species_screening[species_screening$retained, , drop = FALSE]
write.csv(freshwater_species, file.path(output_dir, "freshwater_species.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
if (!nrow(freshwater_species)) stop("No taxa passed the species-screening criteria.")

# ---- 10. Select predictors and prepare repeated training/test partitions -------
# Select predictors from pooled presences + backgrounds before partitioning,
# exactly as in the source workflow. This is not fold-wise predictor selection.
# Split files are produced only for the final retained inland-freshwater taxa.
predictor_selection <- list()
partition_counts <- list()
for (i in seq_len(nrow(freshwater_species))) {
  info.sp <- freshwater_species[i, , drop = FALSE]
  source_file <- file.path(species_dir,
                           paste0("presence_background_", info.sp$fcode, ".RData"))
  # Load into an isolated environment so the saved info.sp cannot overwrite the
  # final screening information. This is a read operation, not a workflow wrapper.
  species_data <- new.env(parent = emptyenv())
  load(source_file, envir = species_data)
  final.DB.P <- species_data$final.DB.P
  final.DB.A <- species_data$final.DB.A
  if (!all(final.DB.P$species == info.sp$common_name) ||
      nrow(final.DB.P) < minimum_presence_records ||
      floor(test_fraction * nrow(final.DB.A)) < 1L) {
    stop("Invalid or insufficient presence/background data for ", info.sp$fcode)
  }
  DB.P <- final.DB.P[, env_names, drop = FALSE]
  DB.A <- final.DB.A[, env_names, drop = FALSE]
  DB.P$OC <- 1L
  DB.A$OC <- 0L
  DB.orginal <- rbind(DB.P, DB.A)  # Legacy object spelling retained deliberately.
  predictors <- DB.orginal[, env_names, drop = FALSE]
  predictor_variance <- vapply(predictors, stats::var, numeric(1))
  constant_predictors <- names(predictor_variance)[predictor_variance == 0]
  predictors <- predictors[, predictor_variance > 0, drop = FALSE]
  if (!ncol(predictors)) stop("No varying predictors for ", info.sp$fcode)
  set.seed(variable_selection_seed + match(info.sp$fcode, sp.list2$fcode))
  if (ncol(predictors) > 1L) {
    selection_cor <- usdm::vifcor(predictors, th = correlation_threshold)
    predictors <- usdm::exclude(predictors, selection_cor)
  }
  if (ncol(predictors) > 1L) {
    selection_vif <- usdm::vifstep(predictors, th = vif_threshold)
    predictors <- usdm::exclude(predictors, selection_vif)
  }
  names.used.env <- names(predictors)
  if (!length(names.used.env)) stop("No predictors survived selection for ", info.sp$fcode)
  predictor_selection[[i]] <- data.frame(
    fcode = info.sp$fcode, variable = env_names,
    constant = env_names %in% constant_predictors,
    selected = env_names %in% names.used.env
  )
  DB.target.sp2 <- DB.orginal[, c("OC", names.used.env), drop = FALSE]
  SDM.formula <- paste0("OC ~ ", paste(names.used.env, collapse = " + "))
  t.DB.P <- DB.target.sp2[DB.target.sp2$OC == 1L, , drop = FALSE]
  t.DB.A <- DB.target.sp2[DB.target.sp2$OC == 0L, , drop = FALSE]
  list.DB.train.10R <- list.DB.test.10R <- vector("list", partition_repeats)
  split_row_indices <- vector("list", partition_repeats)
  for (k in seq_len(partition_repeats)) {
    set.seed(partition_seed + k^2)
    test.num.P <- sample.int(nrow(t.DB.P), floor(test_fraction * nrow(t.DB.P)))
    test.num.A <- sample.int(nrow(t.DB.A), floor(test_fraction * nrow(t.DB.A)))
    list.DB.test.10R[[k]] <- rbind(t.DB.P[test.num.P, , drop = FALSE],
                                   t.DB.A[test.num.A, , drop = FALSE])
    list.DB.train.10R[[k]] <- rbind(t.DB.P[-test.num.P, , drop = FALSE],
                                    t.DB.A[-test.num.A, , drop = FALSE])
    # Indices refer to the presence and background tables 
    # in the matching presence_background_Fxxx.RData checkpoint, not to the raw occurrence table.
    split_row_indices[[k]] <- list(test_presence = test.num.P, test_background = test.num.A)
    partition_counts[[length(partition_counts) + 1L]] <- data.frame(
      fcode = info.sp$fcode, replicate = k,
      n_train_presence = nrow(t.DB.P) - length(test.num.P),
      n_test_presence = length(test.num.P),
      n_train_background = nrow(t.DB.A) - length(test.num.A),
      n_test_background = length(test.num.A)
    )
  }
  save(info.sp, DB.orginal, names.used.env, SDM.formula,
       list.DB.train.10R, list.DB.test.10R, split_row_indices,
       file = file.path(model_input_dir, paste0("model_inputs_", info.sp$fcode, ".RData")))
}
write.csv(do.call(rbind, predictor_selection),
          file.path(output_dir, "predictor_selection.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(do.call(rbind, partition_counts),
          file.path(output_dir, "partition_summary.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 11. Record the computational environment ---------------------------------
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
package_versions <- data.frame(package = required_packages,
                               version = vapply(required_packages, function(x) as.character(utils::packageVersion(x)),
                                                character(1)))
write.csv(package_versions, file.path(output_dir, "package_versions.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
message("Data preparation complete. Retained taxa: ", nrow(freshwater_species),
        ". Output directory: ", output_dir)
