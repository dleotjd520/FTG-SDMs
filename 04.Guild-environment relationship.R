# =============================================================================
# 04. Guild-environment relationship
# =============================================================================
# No plotting, SDM fitting, thermal-guild reassignment, or new PDP calculation.
# All settings are below; no configuration file or workflow function is needed.

# ---- 0. Paths, package, and analysis settings --------------------------------
# Start R in the same project directory used by scripts 01-03.
project_dir <- getwd()
preprocessing_dir <- file.path(project_dir, "outputs", "data_preprocessing")
species_list_file <- file.path(preprocessing_dir, "freshwater_species.csv")
species_data_dir <- file.path(preprocessing_dir, "species_datasets")
species_file_template <- "presence_background_%s.RData"  # %s = fcode
guild_file <- file.path(
  project_dir, "outputs", "thermal_guild_classification",
  "species_thermal_guilds.csv"
)
output_dir <- file.path(project_dir, "outputs", "guild_environment_relationship")
input_encoding <- "UTF-8"
output_encoding <- "UTF-8"

env.vars <- c(paste0("bio", c(2, 3, 8, 9, 14, 15, 18, 19)), "elev")
coord.vars <- c("lon", "lat")
guild_levels <- c("cold-affinity", "intermediate-affinity", "warm-affinity")
envfit_seed <- 520L
envfit_permutations <- 999L

# The source aggregates existing PDP CSVs; it does not generate those CSVs.
# Supply archived PDP inputs here. Script 02 does not generate them either.
# Set FALSE to run only the PCA and guild-environment tests (Sections 1-6).
summarise_existing_pdp <- TRUE
pdp_input_dir <- file.path(project_dir, "inputs", "pdp")
pdp_file_pattern <- "^PDP_.*\\.csv$"
pdp_input_encoding <- "EUC-KR"  # Source CSV encoding; use "UTF-8" if appropriate.
used.vars <- env.vars

# Install once in your local R console if necessary:
# install.packages("vegan")
if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("Package 'vegan' is required. Install it before running this script.")
}
library(vegan)

# ---- 1. Read species metadata and existing thermal-guild assignments -----------
if (!file.exists(species_list_file) || !file.exists(guild_file)) {
  stop("Species metadata are missing. Check species_list_file and guild_file.")
}
if (!dir.exists(species_data_dir)) {
  stop("Species dataset directory not found: ", species_data_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

DB.sp2 <- read.csv(
  species_list_file, fileEncoding = input_encoding,
  stringsAsFactors = FALSE, check.names = FALSE
)
fish.info <- read.csv(
  guild_file, fileEncoding = input_encoding,
  stringsAsFactors = FALSE, check.names = FALSE
)

# Accept legacy metadata headers without changing the taxon-name values.
names(DB.sp2)[names(DB.sp2) == "\uad6d\uba85"] <- "common_name"
names(DB.sp2)[names(DB.sp2) == "\ud559\uba85"] <- "scientific_name"
names(fish.info)[names(fish.info) == "\uad6d\uba85"] <- "common_name"
names(fish.info)[names(fish.info) == "\ud559\uba85"] <- "scientific_name"

required_metadata <- c("fcode", "common_name", "scientific_name")
if (!all(required_metadata %in% names(DB.sp2)) ||
    !all(c(required_metadata, "thermal.guild") %in% names(fish.info))) {
  stop("Metadata require fcode, common_name, and scientific_name; ",
       "the guild table also requires thermal.guild.")
}
if (!nrow(DB.sp2) || !nrow(fish.info) ||
    anyDuplicated(names(DB.sp2)) || anyDuplicated(names(fish.info)) ||
    anyNA(DB.sp2[, required_metadata]) ||
    anyNA(fish.info[, c(required_metadata, "thermal.guild")]) ||
    anyDuplicated(DB.sp2$fcode) || anyDuplicated(DB.sp2$common_name) ||
    anyDuplicated(fish.info$fcode) || anyDuplicated(fish.info$common_name)) {
  stop("Metadata must be non-empty with unique headers, codes, and common-name ",
       "keys, and no missing required values.")
}
for (field in required_metadata) {
  DB.sp2[[field]] <- as.character(DB.sp2[[field]])
  fish.info[[field]] <- as.character(fish.info[[field]])
  if (any(!nzchar(trimws(DB.sp2[[field]]))) ||
      any(!nzchar(trimws(fish.info[[field]])))) {
    stop("Empty metadata values in: ", field)
  }
}

guild_index <- match(DB.sp2$fcode, fish.info$fcode)
if (anyNA(guild_index)) {
  stop("Some retained species have no thermal-guild assignment.")
}
if (!all(DB.sp2$common_name == fish.info$common_name[guild_index])) {
  stop("Species codes and common-name keys disagree between the metadata files.")
}
DB.sp2$thermal.guild <- as.character(fish.info$thermal.guild[guild_index])
if (!all(DB.sp2$thermal.guild %in% guild_levels) ||
    length(unique(DB.sp2$thermal.guild)) < 2L) {
  stop("Use the expected guild labels and include at least two represented guilds.")
}

species_files <- file.path(
  species_data_dir, sprintf(species_file_template, DB.sp2$fcode)
)
if (any(!file.exists(species_files))) {
  stop("Missing species datasets: ",
       paste(basename(species_files[!file.exists(species_files)]), collapse = ", "))
}

# ---- 2. Pool processed background and presence records -----------------------
background_tables <- presence_tables <- vector("list", nrow(DB.sp2))
for (i in seq_len(nrow(DB.sp2))) {
  message("Reading species dataset: ", DB.sp2$fcode[i])

  # Isolate loaded objects so a missing object cannot reuse a previous species.
  species_data <- new.env(parent = emptyenv())
  load(species_files[i], envir = species_data)
  if (!all(c("final.DB.A", "final.DB.P") %in% ls(species_data))) {
    stop("Missing final.DB.A or final.DB.P in: ", species_files[i])
  }
  final.DB.A <- species_data$final.DB.A
  final.DB.P <- species_data$final.DB.P
  if (!is.data.frame(final.DB.A) || !is.data.frame(final.DB.P) ||
      !nrow(final.DB.A) || !nrow(final.DB.P) ||
      !all(c(coord.vars, env.vars) %in% names(final.DB.A)) ||
      !all(c("species", coord.vars, env.vars) %in% names(final.DB.P))) {
    stop("Invalid presence/background tables or missing fields for: ",
         DB.sp2$fcode[i])
  }
  if (anyNA(final.DB.P$species) ||
      !all(as.character(final.DB.P$species) == DB.sp2$common_name[i])) {
    stop("Presence species labels do not match the metadata: ", DB.sp2$fcode[i])
  }

  DB.A.tmp <- final.DB.A[, c(coord.vars, env.vars), drop = FALSE]
  DB.P.tmp <- final.DB.P[, c("species", coord.vars, env.vars), drop = FALSE]
  DB.P.tmp$thermal.guild <- DB.sp2$thermal.guild[i]
  background_tables[[i]] <- DB.A.tmp
  presence_tables[[i]] <- DB.P.tmp
}

# Preserve the species-list order and within-file record order.
DB.A.all <- do.call(rbind, background_tables)
DB.P.all <- do.call(rbind, presence_tables)
rm(background_tables, presence_tables, species_data, final.DB.A, final.DB.P,
   DB.A.tmp, DB.P.tmp)

if (!all(vapply(DB.A.all[, c(coord.vars, env.vars)], is.numeric, logical(1))) ||
    !all(vapply(DB.P.all[, c(coord.vars, env.vars)], is.numeric, logical(1)))) {
  stop("Coordinates and environmental predictors must be numeric.")
}

# ---- 3. Remove missing environmental values and duplicate coordinates --------
input_counts <- data.frame(
  stage = "pooled", n_background = nrow(DB.A.all), n_presence = nrow(DB.P.all)
)
DB.A.all <- DB.A.all[
  complete.cases(DB.A.all[, env.vars, drop = FALSE]), , drop = FALSE
]
DB.P.all <- DB.P.all[
  complete.cases(DB.P.all[, env.vars, drop = FALSE]), , drop = FALSE
]
input_counts <- rbind(input_counts, data.frame(
  stage = "complete_environment",
  n_background = nrow(DB.A.all), n_presence = nrow(DB.P.all)
))

# Background: retain the first record at each coordinate.
# Presence: retain the first record for each species-coordinate combination.
# Different species at the same coordinates are retained.
DB.A.all <- DB.A.all[
  !duplicated(DB.A.all[, coord.vars, drop = FALSE]), , drop = FALSE
]
DB.P.all <- DB.P.all[
  !duplicated(DB.P.all[, c("species", coord.vars), drop = FALSE]), , drop = FALSE
]
rownames(DB.A.all) <- rownames(DB.P.all) <- NULL
input_counts <- rbind(input_counts, data.frame(
  stage = "coordinate_deduplication",
  n_background = nrow(DB.A.all), n_presence = nrow(DB.P.all)
))

if (nrow(DB.A.all) < 3L ||
    any(!is.finite(as.matrix(DB.A.all[, c(coord.vars, env.vars)]))) ||
    any(!is.finite(as.matrix(DB.P.all[, c(coord.vars, env.vars)])))) {
  stop("At least three background rows and finite retained coordinates/predictors ",
       "are required. Resolve invalid inputs upstream.")
}
missing_species <- setdiff(DB.sp2$common_name, unique(DB.P.all$species))
if (length(missing_species)) {
  stop("No complete environmental presence records remain for: ",
       paste(missing_species, collapse = ", "))
}
background_sd <- vapply(DB.A.all[, env.vars], stats::sd, numeric(1))
if (any(!is.finite(background_sd) | background_sd <= 0)) {
  stop("A background predictor has zero or invalid variance; no variable ",
       "is removed automatically.")
}

# ---- 4. Environmental PCA, presence projection, and environmental fitting -----
pca.env <- stats::prcomp(DB.A.all[, env.vars, drop = FALSE], scale. = TRUE)
if (length(pca.env$sdev) < 2L || pca.env$sdev[2L] <= 0) {
  stop("The environmental PCA must provide two nonzero components.")
}
# No temperature-based PC1 sign adjustment is applied in this analysis.
presence_scores <- stats::predict(
  pca.env, newdata = DB.P.all[, env.vars, drop = FALSE]
)
sp_centroid <- stats::aggregate(
  presence_scores[, c("PC1", "PC2"), drop = FALSE],
  by = list(species = DB.P.all$species, thermal.guild = DB.P.all$thermal.guild),
  FUN = mean, na.rm = TRUE
)
metadata_index <- match(sp_centroid$species, DB.sp2$common_name)
sp_centroid$fcode <- DB.sp2$fcode[metadata_index]
sp_centroid$scientific_name <- DB.sp2$scientific_name[metadata_index]

# Fit to presence-level scores and predictor values, not to species centroids.
# choices = 1:2 makes the source's default axis selection explicit.
set.seed(envfit_seed)
fit <- vegan::envfit(
  presence_scores, DB.P.all[, env.vars, drop = FALSE],
  permutations = envfit_permutations, choices = c(1, 2)
)
envfit_results <- data.frame(
  variable = rownames(fit$vectors$arrows),
  PC1_direction = fit$vectors$arrows[, 1],
  PC2_direction = fit$vectors$arrows[, 2],
  r_squared = unname(fit$vectors$r),
  p_value = unname(fit$vectors$pvals),
  row.names = NULL
)

pca_importance <- summary(pca.env)$importance
pca_variance <- data.frame(
  component = colnames(pca_importance),
  standard_deviation = pca_importance[1, ],
  variance_fraction = pca_importance[2, ],
  cumulative_variance_fraction = pca_importance[3, ],
  row.names = NULL
)
pca_loadings <- data.frame(
  variable = rownames(pca.env$rotation), pca.env$rotation, row.names = NULL
)
# Row alignment with presence_scores is retained in the saved metadata.
presence_metadata <- DB.P.all[
  , c("species", "thermal.guild", coord.vars), drop = FALSE
]
presence_metadata$fcode <- DB.sp2$fcode[
  match(presence_metadata$species, DB.sp2$common_name)
]

# ---- 5. Species-level environmental means and between-guild tests -------------
# Each species contributes one mean value per variable to the guild comparisons.
species_env_mean <- stats::aggregate(
  DB.P.all[, env.vars, drop = FALSE],
  by = list(species = DB.P.all$species, thermal.guild = DB.P.all$thermal.guild),
  FUN = mean, na.rm = TRUE
)
metadata_index <- match(species_env_mean$species, DB.sp2$common_name)
species_env_mean$fcode <- DB.sp2$fcode[metadata_index]
species_env_mean$scientific_name <- DB.sp2$scientific_name[metadata_index]

kw_rows <- vector("list", length(env.vars))
names(kw_rows) <- env.vars
for (v in env.vars) {
  dat <- species_env_mean[, c("thermal.guild", v)]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  kt <- stats::kruskal.test(dat[[v]] ~ dat$thermal.guild)
  kw_rows[[v]] <- data.frame(
    variable = v, statistic = unname(kt$statistic),
    df = unname(kt$parameter), p.value = kt$p.value
  )
}
kw.results <- do.call(rbind, kw_rows)
rownames(kw.results) <- NULL
kw.results$p.adj <- stats::p.adjust(kw.results$p.value, method = "BH")

species_env_mean$thermal.guild <- factor(
  species_env_mean$thermal.guild, levels = guild_levels
)
pairwise.list <- list()
for (v in env.vars) {
  dat <- species_env_mean[, c("thermal.guild", v)]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  names(dat) <- c("thermal.guild", "value")
  if (length(unique(dat$thermal.guild)) < 2L) next

  pw <- stats::pairwise.wilcox.test(
    x = dat$value, g = dat$thermal.guild,
    p.adjust.method = "BH", exact = FALSE
  )
  pw.df <- as.data.frame(as.table(pw$p.value))
  names(pw.df) <- c("group1", "group2", "p.adj.within.variable")
  pw.df <- pw.df[!is.na(pw.df$p.adj.within.variable), , drop = FALSE]
  pw.df$variable <- rep(v, nrow(pw.df))
  pairwise.list[[v]] <- pw.df[
    , c("variable", "group1", "group2", "p.adj.within.variable"), drop = FALSE
  ]
}
pairwise.results <- do.call(rbind, pairwise.list)
rownames(pairwise.results) <- NULL

# SOURCE REPRODUCTION: apply BH again to the already adjusted values.
# The inherited column name below is retained for compatibility.
# This is NOT one global BH adjustment applied to raw pairwise p-values.
pairwise.results$p.adj.global.BH <- stats::p.adjust(
  pairwise.results$p.adj.within.variable, method = "BH"
)
pairwise.results <- pairwise.results[
  order(pairwise.results$variable, pairwise.results$p.adj.within.variable),
  , drop = FALSE
]
if (any(!is.finite(kw.results$p.value)) || !nrow(pairwise.results)) {
  warning("Some tests have undefined p-values or no finite pairwise results; ",
          "inspect the exported test tables before interpretation.")
}

# ---- 6. Save PCA, environmental-fit, and statistical results ------------------
write.csv(sp_centroid,
          file.path(output_dir, "species_environment_pca_centroids.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(species_env_mean,
          file.path(output_dir, "species_environment_means.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(pca_variance, file.path(output_dir, "environment_pca_variance.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(pca_loadings, file.path(output_dir, "environment_pca_loadings.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(envfit_results, file.path(output_dir, "environment_envfit_results.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(kw.results, file.path(output_dir, "guild_kruskal_wallis_results.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(pairwise.results, file.path(output_dir, "guild_pairwise_wilcoxon_results.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(input_counts, file.path(output_dir, "guild_environment_data_counts.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

save(
  DB.sp2, pca.env, presence_scores, presence_metadata, sp_centroid,
  species_env_mean, fit, envfit_results, kw.results, pairwise.results,
  pca_variance, pca_loadings, input_counts, env.vars, coord.vars, guild_levels,
  envfit_seed, envfit_permutations,
  file = file.path(output_dir, "guild_environment_analysis.RData")
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
message("PCA and guild-environment test results saved: ", output_dir)

# ---- 7. Summarise existing species PDP files by environmental variable --------
# This optional section only reads and aggregates existing numerical PDP data.
# No predictions, smoothing, interpolation, normalisation, or figures are added.

if (summarise_existing_pdp) {
  if (!dir.exists(pdp_input_dir)) {
    stop("PDP input directory not found: ", pdp_input_dir,
         ". Supply existing PDP CSVs or set summarise_existing_pdp <- FALSE. ",
         "PCA and statistical results have already been saved.")
  }
  pdp_files <- list.files(
    pdp_input_dir, pattern = pdp_file_pattern,
    full.names = TRUE, ignore.case = TRUE
  )
  if (!length(pdp_files)) {
    stop("No existing PDP CSVs matched pdp_file_pattern. ",
         "The supplied source does not generate these inputs.")
  }

  # Accepted examples: PDP_F001.csv or PDP_F001_MaxEnt.csv.
  # The second underscore-delimited filename token identifies the species.
  pdp_stems <- tools::file_path_sans_ext(basename(pdp_files))
  pdp_fcodes <- sub("^PDP_([^_]+).*", "\\1", pdp_stems, ignore.case = TRUE)
  if (anyDuplicated(pdp_fcodes) || any(!nzchar(pdp_fcodes)) ||
      any(pdp_fcodes == pdp_stems)) {
    stop("Use exactly one PDP CSV per species, named PDP_<fcode>[_suffix].csv. ",
         "Do not mix alternative model files for the same species.")
  }

  # Preserve the source's scope: all matching files are processed.
  # Extra taxa are flagged, not silently removed; unknown metadata remain NA.
  pdp_input_manifest <- data.frame(
    file = basename(pdp_files), fcode = pdp_fcodes,
    included_in_guild_analysis = pdp_fcodes %in% DB.sp2$fcode,
    n_input_rows = integer(length(pdp_files)),
    stringsAsFactors = FALSE
  )
  missing_pdp_species <- setdiff(DB.sp2$fcode, pdp_fcodes)
  if (length(missing_pdp_species)) {
    warning("No PDP file supplied for guild-analysis species: ",
            paste(missing_pdp_species, collapse = ", "))
  }
  if (any(!pdp_input_manifest$included_in_guild_analysis)) {
    warning("Some PDP files refer to taxa outside the guild analysis. ",
            "All supplied files are retained, as in the source; see the manifest.")
  }

  list.var.PDP <- vector("list", length(used.vars))
  names(list.var.PDP) <- used.vars
  for (v in used.vars) {
    list.var.PDP[[v]] <- data.frame(
      value.x = numeric(), value.y = numeric(), fcode = character()
    )
  }
  pdp_coverage_rows <- list()

  # Read each file once, rather than re-reading it separately for every variable.
  # The within-file arithmetic means and per-variable file order are unchanged.
  for (k in seq_along(pdp_files)) {
    message("Summarising existing PDP file: ", basename(pdp_files[k]))
    total.DB <- read.csv(
      pdp_files[k], fileEncoding = pdp_input_encoding,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    if (!all(c("variable", "value.x", "value.y") %in% names(total.DB)) ||
        anyDuplicated(names(total.DB)) ||
        !is.numeric(total.DB$value.x) || !is.numeric(total.DB$value.y)) {
      stop("PDP CSVs require variable and numeric value.x/value.y columns: ",
           pdp_files[k])
    }
    pdp_input_manifest$n_input_rows[k] <- nrow(total.DB)

    for (v in used.vars) {
      PDP.DB <- total.DB[
        !is.na(total.DB$variable) & total.DB$variable == v,
        c("value.x", "value.y"), drop = FALSE
      ]
      if (any(is.infinite(as.matrix(PDP.DB)))) {
        stop("Infinite PDP values require upstream review: ", pdp_files[k], " / ", v)
      }
      complete_pdp <- complete.cases(PDP.DB)
      pdp_coverage_rows[[length(pdp_coverage_rows) + 1L]] <- data.frame(
        file = basename(pdp_files[k]), fcode = pdp_fcodes[k], variable = v,
        n_input_rows = nrow(PDP.DB), n_complete_rows = sum(complete_pdp),
        n_unique_x = length(unique(PDP.DB$value.x[complete_pdp]))
      )
      # Empty species-variable combinations are recorded, not fabricated.
      if (!any(complete_pdp)) next

      # Original calculation: an unweighted mean at each exact value.x.
      PDP.DB2 <- stats::aggregate(
        value.y ~ value.x, data = PDP.DB, FUN = mean, na.action = stats::na.omit
      )
      PDP.DB2$fcode <- pdp_fcodes[k]
      list.var.PDP[[v]] <- rbind(list.var.PDP[[v]], PDP.DB2)
    }
  }

  pdp_coverage <- do.call(rbind, pdp_coverage_rows)
  empty_pdp_variables <- used.vars[vapply(list.var.PDP, nrow, integer(1)) == 0L]
  if (length(empty_pdp_variables)) {
    warning("No complete PDP data for: ", paste(empty_pdp_variables, collapse = ", "),
            ". Empty entries are retained in list.var.PDP.")
  }
  pdp_long_parts <- vector("list", length(used.vars))
  for (j in seq_along(used.vars)) {
    pdp_part <- list.var.PDP[[j]]
    pdp_part$variable <- rep(used.vars[j], nrow(pdp_part))
    pdp_long_parts[[j]] <- pdp_part[
      , c("variable", "fcode", "value.x", "value.y"), drop = FALSE
    ]
  }
  pdp_summary <- do.call(rbind, pdp_long_parts)
  rownames(pdp_summary) <- NULL
  pdp_metadata_index <- match(pdp_summary$fcode, DB.sp2$fcode)
  pdp_summary$common_name <- DB.sp2$common_name[pdp_metadata_index]
  pdp_summary$scientific_name <- DB.sp2$scientific_name[pdp_metadata_index]
  pdp_summary$thermal.guild <- DB.sp2$thermal.guild[pdp_metadata_index]
  pdp_summary$included_in_guild_analysis <- !is.na(pdp_metadata_index)

  write.csv(pdp_summary, file.path(output_dir, "species_pdp_summary.csv"),
            row.names = FALSE, fileEncoding = output_encoding)
  write.csv(pdp_input_manifest, file.path(output_dir, "pdp_input_manifest.csv"),
            row.names = FALSE, fileEncoding = output_encoding)
  write.csv(pdp_coverage, file.path(output_dir, "pdp_variable_coverage.csv"),
            row.names = FALSE, fileEncoding = output_encoding)
  save(
    list.var.PDP, used.vars, pdp_summary, pdp_input_manifest, pdp_coverage,
    missing_pdp_species, empty_pdp_variables,
    file = file.path(output_dir, "species_pdp_by_variable.RData")
  )
  message("Existing PDP summaries saved: ", output_dir)
}

# ---- 8. Record the R session -------------------------------------------------
# Refresh after the optional PDP section; no analysis or configuration is rerun.
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
message("Guild-environment analysis complete. Output directory: ", output_dir)
