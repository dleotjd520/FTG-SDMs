# =============================================================================
# 03. Thermal guild classification
# =============================================================================
# Run the numbered sections sequentially in an interactive R session.
# All settings are defined here; no configuration file or workflow wrapper is used.
# Only the standard R installation is required.

# ---- 0. Paths and analysis settings ------------------------------------------
# Start R in the same project directory used for data preprocessing.
project_dir <- getwd()
preprocessing_dir <- file.path(project_dir, "outputs", "data_preprocessing")
species_list_file <- file.path(preprocessing_dir, "freshwater_species.csv")
species_data_dir <- file.path(preprocessing_dir, "species_datasets")
species_file_template <- "presence_background_%s.RData"  # %s is the fcode.
output_dir <- file.path(project_dir, "outputs", "thermal_guild_classification")

input_encoding <- "UTF-8"  # Matches the preceding preprocessing script's outputs.
output_encoding <- "UTF-8"

# Retain the five temperature predictors used in the source script.
# bio1: annual mean temperature
# bio5: maximum temperature of the warmest month
# bio6: minimum temperature of the coldest month
# bio10: mean temperature of the warmest quarter
# bio11: mean temperature of the coldest quarter
# bio4 and bio7 are not included in this classification.
env.vars <- paste0("bio", c(1, 5, 6, 10, 11))
coord.vars <- c("lon", "lat")
sp.var <- "species"
n_clusters <- 5L

# Preserve the source's explicit aggregation of five clusters into three guilds.
# T1-T5 are ordered by increasing mean species-centroid PC1.
# This is a supplied mapping, not an automatically estimated guild boundary.
# For a different dataset, review cluster.profile before applying this mapping.
guild.map <- data.frame(
  cluster5 = c("T1", "T2", "T3", "T4", "T5"),
  thermal.guild = c(
    "cold-affinity", "cold-affinity", "intermediate-affinity",
    "warm-affinity", "warm-affinity"
  ),
  stringsAsFactors = FALSE
)

# ---- 1. Read the retained freshwater species list -----------------------------
if (!file.exists(species_list_file)) {
  stop("Species list not found: ", species_list_file)
}
if (!dir.exists(species_data_dir)) {
  stop("Species dataset directory not found: ", species_data_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

DB.sp2 <- read.csv(
  species_list_file,
  fileEncoding = input_encoding,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Accept the two legacy metadata headers without translating taxon identifiers.
# Unicode escapes keep the script itself in English/ASCII.
names(DB.sp2)[names(DB.sp2) == "\uad6d\uba85"] <- "common_name"
names(DB.sp2)[names(DB.sp2) == "\ud559\uba85"] <- "scientific_name"

required_fields <- c("fcode", "common_name", "scientific_name")
missing_fields <- setdiff(required_fields, names(DB.sp2))
if (length(missing_fields) > 0L) {
  stop("Species list is missing fields: ", paste(missing_fields, collapse = ", "))
}
if (anyDuplicated(names(DB.sp2)) ||
    anyNA(DB.sp2[, required_fields, drop = FALSE]) ||
    anyDuplicated(DB.sp2$fcode) || anyDuplicated(DB.sp2$common_name)) {
  stop("Species metadata require unique headers, codes, and common-name keys, ",
       "with no missing required values.")
}
for (field in required_fields) {
  if (any(!nzchar(trimws(as.character(DB.sp2[[field]]))))) {
    stop("Empty species metadata values in: ", field)
  }
}
if (nrow(DB.sp2) < n_clusters) {
  stop("At least ", n_clusters, " species are required for the supplied clustering.")
}
# The input is already screened upstream. Do not select taxa again here.
classification_fields <- c("PC1.centroid", "PC2.centroid", "cluster5.raw",
                           "cluster5", "thermal.guild")
if (any(classification_fields %in% names(DB.sp2))) {
  stop("Use the screened species list, not a previous classification output.")
}

species_files <- file.path(
  species_data_dir,
  sprintf(species_file_template, DB.sp2$fcode)
)
if (any(!file.exists(species_files))) {
  stop("Missing species datasets: ",
       paste(basename(species_files[!file.exists(species_files)]), collapse = ", "))
}

# ---- 2. Pool the species-specific background and presence datasets ------------
background_tables <- vector("list", nrow(DB.sp2))
presence_tables <- vector("list", nrow(DB.sp2))

for (i in seq_len(nrow(DB.sp2))) {
  # A fresh environment prevents objects from the preceding species being reused
  # if a file does not contain both required tables.
  species_data <- new.env(parent = emptyenv())
  loaded_objects <- load(species_files[i], envir = species_data)
  if (!all(c("final.DB.A", "final.DB.P") %in% loaded_objects)) {
    stop("Required presence/background objects are missing in: ", species_files[i])
  }
  final.DB.A <- species_data$final.DB.A
  final.DB.P <- species_data$final.DB.P
  if (!is.data.frame(final.DB.A) || !is.data.frame(final.DB.P)) {
    stop("Presence/background objects must be data.frames: ", DB.sp2$fcode[i])
  }
  if (!all(c(coord.vars, env.vars) %in% names(final.DB.A)) ||
      !all(c(sp.var, coord.vars, env.vars) %in% names(final.DB.P))) {
    stop("Required coordinate, species, or temperature fields are missing: ",
         DB.sp2$fcode[i])
  }
  if (nrow(final.DB.A) == 0L || nrow(final.DB.P) == 0L) {
    stop("Empty presence/background dataset for retained species: ", DB.sp2$fcode[i])
  }
  if (anyNA(final.DB.P[[sp.var]]) ||
      !all(as.character(final.DB.P[[sp.var]]) == DB.sp2$common_name[i])) {
    stop("Presence labels do not match the species list: ", DB.sp2$fcode[i])
  }
  
  DB.A.tmp <- final.DB.A[, c(coord.vars, env.vars), drop = FALSE]
  DB.P.tmp <- final.DB.P[, c(sp.var, coord.vars, env.vars), drop = FALSE]
  if (!all(vapply(DB.A.tmp, is.numeric, logical(1))) ||
      !all(vapply(DB.P.tmp[, c(coord.vars, env.vars), drop = FALSE],
                  is.numeric, logical(1)))) {
    stop("Coordinates and temperature predictors must be numeric: ", DB.sp2$fcode[i])
  }
  if (any(!is.finite(as.matrix(DB.A.tmp[, coord.vars, drop = FALSE]))) ||
      any(!is.finite(as.matrix(DB.P.tmp[, coord.vars, drop = FALSE])))) {
    stop("Non-finite coordinates must be resolved during preprocessing: ",
         DB.sp2$fcode[i])
  }
  
  background_tables[[i]] <- DB.A.tmp
  presence_tables[[i]] <- DB.P.tmp
}

# Binding once preserves the input species/record order while avoiding repeated reallocation of the pooled tables.
DB.A.all <- do.call(rbind, background_tables)
DB.P.all <- do.call(rbind, presence_tables)
rm(background_tables, presence_tables, species_data, final.DB.A, final.DB.P,
   DB.A.tmp, DB.P.tmp)

# ---- 3. Remove missing temperatures and duplicate background coordinates -------
input_counts <- data.frame(
  stage = "pooled",
  n_background = nrow(DB.A.all),
  n_presence = nrow(DB.P.all)
)

# Apply the source's complete-case rule to the five temperature predictors only.
DB.A.all <- DB.A.all[complete.cases(DB.A.all[, env.vars, drop = FALSE]), , drop = FALSE]
DB.P.all <- DB.P.all[complete.cases(DB.P.all[, env.vars, drop = FALSE]), , drop = FALSE]
input_counts <- rbind(input_counts, data.frame(
  stage = "complete_temperatures",
  n_background = nrow(DB.A.all), n_presence = nrow(DB.P.all)
))

# Keep the first row at each background coordinate, as in the source.
# Do not deduplicate presence records, including records shared across species.
DB.A.all <- DB.A.all[
  !duplicated(DB.A.all[, coord.vars, drop = FALSE]), , drop = FALSE
]
input_counts <- rbind(input_counts, data.frame(
  stage = "unique_background_coordinates",
  n_background = nrow(DB.A.all), n_presence = nrow(DB.P.all)
))

missing_species <- setdiff(DB.sp2$common_name, unique(DB.P.all[[sp.var]]))
if (length(missing_species) > 0L) {
  stop("No complete temperature records remain for: ",
       paste(missing_species, collapse = ", "))
}
if (nrow(DB.A.all) < 2L ||
    any(!is.finite(as.matrix(DB.A.all[, env.vars, drop = FALSE]))) ||
    any(!is.finite(as.matrix(DB.P.all[, env.vars, drop = FALSE])))) {
  stop("PCA requires at least two background rows and finite temperature values.")
}
background_sd <- vapply(DB.A.all[, env.vars, drop = FALSE], stats::sd, numeric(1))
if (any(!is.finite(background_sd) | background_sd <= 0)) {
  stop("Every temperature predictor must vary in the pooled background data. ",
       "No predictors have been dropped automatically.")
}

# ---- 4. Fit a standardized PCA to the pooled background environment ------------
pca.env <- stats::prcomp(DB.A.all[, env.vars, drop = FALSE], scale. = TRUE)
pc1.A <- stats::predict(pca.env, newdata = DB.A.all[, env.vars, drop = FALSE])[, 1]
pc1_bio1_cor <- stats::cor(pc1.A, DB.A.all$bio1, use = "complete.obs")
if (!is.finite(pc1_bio1_cor)) {
  stop("The PC1-Bio1 correlation is undefined; the thermal direction cannot be set.")
}
sign.pc1 <- ifelse(pc1_bio1_cor < 0, -1, 1)

# Keep pca.env unchanged. Apply sign.pc1 to projected PC1 scores below.
# Save the sign explicitly so later projections can use the same orientation.
pca_importance <- summary(pca.env)$importance
pca_variance <- data.frame(
  component = colnames(pca_importance),
  standard_deviation = as.numeric(pca_importance[1, ]),
  proportion_variance = as.numeric(pca_importance[2, ]),
  cumulative_variance = as.numeric(pca_importance[3, ])
)
pca_loadings <- data.frame(
  variable = rownames(pca.env$rotation),
  PC1_raw = unname(pca.env$rotation[, 1]),
  PC1 = unname(sign.pc1 * pca.env$rotation[, 1]),
  PC2 = unname(pca.env$rotation[, 2])
)

# ---- 5. Project presences and calculate species-level centroids ----------------
score.P <- as.data.frame(
  stats::predict(pca.env, newdata = DB.P.all[, env.vars, drop = FALSE])[, 1:2, drop = FALSE]
)
names(score.P) <- c("PC1", "PC2")
score.P$PC1 <- sign.pc1 * score.P$PC1
score.P$common_name <- DB.P.all[[sp.var]]

species.centroid <- stats::aggregate(
  cbind(PC1, PC2) ~ common_name,
  data = score.P,
  FUN = mean
)
names(species.centroid) <- c("common_name", "PC1.centroid", "PC2.centroid")

# ---- 6. Cluster species centroids using Euclidean distance and Ward.D2 ---------
# Use the first two PC coordinates directly, without rescaling either axis.
X <- species.centroid[, c("PC1.centroid", "PC2.centroid"), drop = FALSE]
# Species-code labels make the saved clustering identifiable without changing
# the ordering or numerical values used in the distance calculation.
rownames(X) <- DB.sp2$fcode[match(species.centroid$common_name, DB.sp2$common_name)]
d <- stats::dist(X, method = "euclidean")
hc <- stats::hclust(d, method = "ward.D2")
species.centroid$cluster5.raw <- unname(stats::cutree(hc, k = n_clusters))

# Order the five clusters by their mean species-centroid PC1, from cold to warm.
cluster.order.tmp <- stats::aggregate(
  PC1.centroid ~ cluster5.raw,
  data = species.centroid,
  FUN = mean
)
cluster.order.tmp <- cluster.order.tmp[order(cluster.order.tmp$PC1.centroid), , drop = FALSE]
raw.cluster.order <- as.character(cluster.order.tmp$cluster5.raw)
raw.to.ordered <- setNames(paste0("T", seq_along(raw.cluster.order)), raw.cluster.order)
species.centroid$cluster5 <- unname(
  raw.to.ordered[as.character(species.centroid$cluster5.raw)]
)

# ---- 7. Summarize species and cluster temperature profiles ---------------------
species.env <- stats::aggregate(
  DB.P.all[, env.vars, drop = FALSE],
  by = list(common_name = DB.P.all[[sp.var]]),
  FUN = mean
)
DB.cluster <- merge(species.centroid, species.env, by = "common_name", all.x = TRUE)

# Each species contributes equally to cluster-level means, regardless of its
# number of occurrence records.
cluster.profile <- stats::aggregate(
  DB.cluster[, c("PC1.centroid", "PC2.centroid", env.vars), drop = FALSE],
  by = list(cluster5 = DB.cluster$cluster5),
  FUN = mean
)
cluster.n <- stats::aggregate(common_name ~ cluster5, data = DB.cluster, FUN = length)
names(cluster.n)[2] <- "n_species"
cluster.profile <- merge(cluster.profile, cluster.n, by = "cluster5")
cluster.profile <- cluster.profile[order(cluster.profile$PC1.centroid), , drop = FALSE]
cluster.profile$PC1_gap_to_next <- c(diff(cluster.profile$PC1.centroid), NA_real_)

# Review cluster.profile before proceeding when using a different dataset.
# The following section retains the source's mapping; it does not optimize it.

# ---- 8. Apply the supplied cluster-to-guild mapping ----------------------------
if (anyDuplicated(guild.map$cluster5) || anyNA(guild.map) ||
    !setequal(guild.map$cluster5, species.centroid$cluster5)) {
  stop("guild.map must assign exactly one guild to every observed cluster.")
}
species.centroid <- merge(species.centroid, guild.map, by = "cluster5", all.x = TRUE)
DB.guild <- merge(
  DB.cluster,
  species.centroid[, c("common_name", "thermal.guild"), drop = FALSE],
  by = "common_name", all.x = TRUE
)

# Each species also contributes equally to guild-level means.
guild.profile <- stats::aggregate(
  DB.guild[, c("PC1.centroid", "PC2.centroid", env.vars), drop = FALSE],
  by = list(thermal.guild = DB.guild$thermal.guild),
  FUN = mean
)
guild.n <- stats::aggregate(common_name ~ thermal.guild, data = DB.guild, FUN = length)
names(guild.n)[2] <- "n_species"
guild.profile <- merge(guild.profile, guild.n, by = "thermal.guild")
guild.profile <- guild.profile[order(guild.profile$PC1.centroid), , drop = FALSE]

# ---- 9. Calculate cluster/guild centroids and their distances ------------------
# Retain the source's numerical summaries without drawing any figures.
cluster.centroid <- stats::aggregate(
  cbind(PC1.centroid, PC2.centroid) ~ cluster5,
  data = species.centroid, FUN = mean
)
guild.centroid <- stats::aggregate(
  cbind(PC1.centroid, PC2.centroid) ~ thermal.guild,
  data = species.centroid, FUN = mean
)
rownames(cluster.centroid) <- cluster.centroid$cluster5
rownames(guild.centroid) <- guild.centroid$thermal.guild
cluster.dist <- stats::dist(
  cluster.centroid[, c("PC1.centroid", "PC2.centroid"), drop = FALSE],
  method = "euclidean"
)
guild.dist <- stats::dist(
  guild.centroid[, c("PC1.centroid", "PC2.centroid"), drop = FALSE],
  method = "euclidean"
)

# ---- 10. Save classification tables and reusable analysis objects --------------
species.centroid2 <- merge(DB.sp2, species.centroid, by = "common_name", all.x = TRUE)
if (nrow(species.centroid2) != nrow(DB.sp2) || anyNA(species.centroid2$thermal.guild)) {
  stop("Every input species must receive exactly one thermal-guild assignment.")
}

write.csv(species.centroid2, file.path(output_dir, "species_thermal_guilds.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(cluster.profile, file.path(output_dir, "thermal_cluster_profiles.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(guild.profile, file.path(output_dir, "thermal_guild_profiles.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(guild.map, file.path(output_dir, "thermal_cluster_guild_mapping.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(pca_variance, file.path(output_dir, "thermal_pca_variance.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(pca_loadings, file.path(output_dir, "thermal_pca_loadings.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(input_counts, file.path(output_dir, "thermal_guild_data_counts.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# Pooled input tables and per-presence scores are not saved again.
# They remain available in this session and can be reconstructed from the inputs.
# As in the source, pca.env itself retains the background PCA scores in pca.env$x.
save(
  DB.sp2, species.centroid, species.centroid2, species.env,
  cluster.profile, guild.profile, cluster.centroid, guild.centroid,
  cluster.dist, guild.dist, hc, raw.to.ordered, guild.map,
  pca.env, sign.pc1, pc1_bio1_cor, pca_variance, pca_loadings,
  env.vars, coord.vars, sp.var, n_clusters, input_counts,
  file = file.path(output_dir, "thermal_guild_classification.RData")
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

message("Thermal guild classification complete: ", nrow(species.centroid2),
        " species. Output directory: ", output_dir)
