# =============================================================================
# 02. Species distribution modelling
# =============================================================================
# Source workflow: MaxEnt fitting -> species thresholds -> prediction rasters.
# Inputs follow 01_data_collection_and_preprocessing.R. 
# Thermal guild assignments are not predictors and are not required by this script.

# ---- 0. Paths and settings ---------------------------------------------------
project_dir <- getwd()  # Start R in the project root used for preprocessing.
raw_dir <- file.path(project_dir, "data", "raw")
preprocessing_dir <- file.path(project_dir, "outputs", "data_preprocessing")
species_list_file <- file.path(preprocessing_dir, "freshwater_species.csv")
model_input_dir <- file.path(preprocessing_dir, "model_inputs")
species_data_dir <- file.path(preprocessing_dir, "species_datasets")
model_input_template <- "model_inputs_%s.RData"  # %s is the species fcode.
species_data_template <- "presence_background_%s.RData"

output_dir <- file.path(project_dir, "outputs", "species_distribution_modelling")
model_dir <- file.path(output_dir, "models")
threshold_dir <- file.path(output_dir, "thresholds")
map_dir <- file.path(output_dir, "prediction_maps")
threshold_file <- file.path(threshold_dir, "species_thresholds_and_performance.csv")

climate_dirs <- c(
  current = file.path(raw_dir, "CHELSA_1981_2010"),
  nearf = file.path(raw_dir, "CHELSA_2041_2070"),
  farf = file.path(raw_dir, "CHELSA_2071_2100")
)
elevation_file <- file.path(raw_dir, "wc2.1_30s_elev.tif")
boundary_file <- file.path(raw_dir, "korea_boundary_buffer.shp")
water_mask_file <- file.path(raw_dir, "inland_water_mask_1km.tif")
input_encoding <- "UTF-8"
output_encoding <- "UTF-8"
overwrite_output <- FALSE  # Deliberately enable before replacing existing results.

# Model settings retained from the supplied MaxEnt script.
set.num <- 910520L
k.fold.num <- 5L
regmult_bounds <- c(0.1, 2)
bayes_initial_points <- 5L
bayes_iterations <- 5L
bayes_batch_size <- 5L
default_regmult <- 1
regmult_retry_step <- 0.1
maximum_fit_attempts <- 100L  # Safety limit replacing the source's unlimited loop.
optimization_cores <- 5L     # Set to 1 for sequential Bayesian optimization.
prediction_cores <- 5L

# Original domestic-occurrence bounding box; boundaries are excluded.
roi_bounds <- c(xmin = 125, xmax = 130, ymin = 33, ymax = 39)
roi_minimum_coverage <- 0

# IMPORTANT: the source ROI helper passes predicted labels as y_true and observed
# labels as y_pred to MLmetrics::Sensitivity()/Specificity(). This is not ordinary
# TSS. Preserve that calculation by default to avoid silently changing results.
# Use "standard_tss" ONLY after deciding to correct the original method and
# rerunning/rechecking affected thresholds and binary maps. Final performance
# tables always use the conventional observed/predicted argument order.
roi_score_convention <- "source_argument_order"

# Full workflow: current climate plus all six future scenario/period combinations.
scenarios <- rbind(
  data.frame(sce = "normal", period = "current"),
  expand.grid(sce = c("ssp126", "ssp370", "ssp585"),
              period = c("nearf", "farf"), stringsAsFactors = FALSE)
)
scenario_rows <- seq_len(nrow(scenarios))

# ---- 1. Check packages and read the species manifest -------------------------
required_packages <- c("maxnet", "glmnet", "dismo", "caret", "MLmetrics",
                       "ParBayesianOptimization", "doParallel", "foreach", "terra")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install the required R package before continuing: ", package)
  }
}
library(maxnet)
library(terra)

if (!roi_score_convention %in% c("source_argument_order", "standard_tss")) {
  stop("roi_score_convention must be 'source_argument_order' or 'standard_tss'.")
}
if (roi_minimum_coverage < 0 || roi_minimum_coverage > 1 ||
    k.fold.num < 2L || maximum_fit_attempts < 1L ||
    optimization_cores < 1L || prediction_cores < 1L ||
    length(regmult_bounds) != 2L || any(!is.finite(regmult_bounds)) ||
    regmult_bounds[1] <= 0 || regmult_bounds[2] <= regmult_bounds[1]) {
  stop("Invalid model, parallelization, or ROI settings in Section 0.")
}
if (!file.exists(species_list_file)) stop("Species list not found: ", species_list_file)
for (directory in c(output_dir, model_dir, threshold_dir, map_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

species_list <- read.csv(species_list_file, fileEncoding = input_encoding,
                         stringsAsFactors = FALSE, check.names = FALSE)
# Optional legacy-header support. Taxon names themselves are not translated.
names(species_list)[names(species_list) == "\uad6d\uba85"] <- "common_name"
names(species_list)[names(species_list) == "\ud559\uba85"] <- "scientific_name"
required_fields <- c("fcode", "common_name", "scientific_name")
if (!all(required_fields %in% names(species_list)) ||
    anyDuplicated(names(species_list)) || !nrow(species_list)) {
  stop("The species list must contain fcode, common_name, and scientific_name.")
}
if (anyNA(species_list[, required_fields, drop = FALSE]) ||
    anyDuplicated(species_list$fcode)) {
  stop("Species codes must be unique and required metadata must not be missing.")
}
if (any(!grepl("^[A-Za-z0-9_-]+$", species_list$fcode))) {
  stop("Use non-empty, filename-safe species codes (letters, numbers, '_' or '-').")
}
# Use the already screened input list without applying another taxon filter.
# Stable file-code order matches the source's sorted file traversal convention.
species_list <- species_list[order(species_list$fcode), , drop = FALSE]
rownames(species_list) <- NULL
model_input_files <- file.path(model_input_dir,
                               sprintf(model_input_template, species_list$fcode))
species_data_files <- file.path(species_data_dir,
                                sprintf(species_data_template, species_list$fcode))
model_files <- file.path(model_dir,
                         paste0("SDM_output_", species_list$fcode, "_MaxEnt.RData"))
if (any(!file.exists(c(model_input_files, species_data_files)))) {
  missing_files <- c(model_input_files, species_data_files)
  stop("Missing preprocessed inputs: ",
       paste(missing_files[!file.exists(missing_files)], collapse = "; "))
}
write.csv(species_list, file.path(output_dir, "modelled_species.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 2. Shared performance calculation ---------------------------------------
# This small helper replaces the repeated metric definitions in the source files.
# OC = 1 denotes presence; OC = 0 denotes background, not confirmed absence.
evaluate_performance <- function(y_true, y_pred, y_class) {
  if (length(y_true) != length(y_pred) || length(y_true) != length(y_class) ||
      !length(y_true) || anyNA(y_true) || anyNA(y_class) ||
      any(!is.finite(y_pred)) || !all(y_true %in% c(0, 1)) ||
      !all(y_class %in% c(0, 1)) || length(unique(y_true)) != 2L) {
    stop("Performance calculation requires aligned finite data and both OC classes.")
  }
  sensitivity <- MLmetrics::Sensitivity(y_true = y_true, y_pred = y_class, positive = "1")
  specificity <- MLmetrics::Specificity(y_true = y_true, y_pred = y_class, positive = "1")
  tp <- as.numeric(sum(y_true == 1 & y_class == 1))
  tn <- as.numeric(sum(y_true == 0 & y_class == 0))
  fp <- as.numeric(sum(y_true == 0 & y_class == 1))
  fn <- as.numeric(sum(y_true == 1 & y_class == 0))
  denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (denominator == 0) NA_real_ else (tp * tn - fp * fn) / denominator
  metrics <- list(
    ACC = MLmetrics::Accuracy(y_true = y_true, y_pred = y_class),
    AUC = MLmetrics::AUC(y_true = y_true, y_pred = y_pred),
    PRAUC = MLmetrics::PRAUC(y_true = y_true, y_pred = y_pred),
    TSS = sensitivity + specificity - 1,
    F1.1 = MLmetrics::F1_Score(y_true = y_true, y_pred = y_class, positive = "1"),
    Sens.1 = sensitivity,
    Spec.1 = specificity,
    Prec.1 = MLmetrics::Precision(y_true = y_true, y_pred = y_class, positive = "1"),
    Recall.1 = MLmetrics::Recall(y_true = y_true, y_pred = y_class, positive = "1"),
    MCC = mcc
  )
  # Undefined rates are stored as NA rather than being interpreted as zero.
  metrics <- vapply(metrics, function(value) {
    if (length(value) == 1L && is.finite(value)) as.numeric(value) else NA_real_
  }, numeric(1))
  as.data.frame(as.list(metrics), check.names = FALSE)
}

# ---- 3. Fit MaxEnt models for the supplied training/test partitions -----------
for (i in seq_len(nrow(species_list))) {
  info.sp <- species_list[i, , drop = FALSE]
  sp.code <- as.character(info.sp$fcode)
  message("Fitting species ", i, "/", nrow(species_list), ": ", sp.code)
  if (file.exists(model_files[i]) && !overwrite_output) {
    stop("Model output already exists: ", model_files[i],
         ". Skip Section 3 to reuse it, or explicitly enable overwrite_output.")
  }
  model_input <- new.env(parent = emptyenv())
  loaded_objects <- load(model_input_files[i], envir = model_input)
  required_objects <- c("info.sp", "DB.orginal", "names.used.env", "SDM.formula",
                        "list.DB.train.10R", "list.DB.test.10R")
  if (!all(required_objects %in% loaded_objects) ||
      !identical(as.character(model_input$info.sp$fcode), sp.code)) {
    stop("Missing or mismatched model-input objects for ", sp.code)
  }
  DB.orginal <- model_input$DB.orginal  # Original object spelling is retained.
  names.used.env <- model_input$names.used.env
  SDM.formula <- model_input$SDM.formula
  list.DB.train.10R <- model_input$list.DB.train.10R
  list.DB.test.10R <- model_input$list.DB.test.10R
  split_row_indices <- model_input$split_row_indices  # NULL for legacy files.
  n_repeats <- length(list.DB.train.10R)
  if (!n_repeats || length(list.DB.test.10R) != n_repeats ||
      !length(names.used.env) || anyDuplicated(names.used.env)) {
    stop("Invalid partition lists or predictor names for ", sp.code)
  }
  # Ten partitions are expected from the preceding script; no new split is made.
  list.SDM <- vector("list", n_repeats)
  cv_folds <- vector("list", n_repeats)
  optimization_history <- vector("list", n_repeats)
  evaluation_list <- prediction_list <- parameter_list <- vector("list", n_repeats)
  
  for (k in seq_len(n_repeats)) {
    message("  Partition ", k, "/", n_repeats)
    train.DB <- as.data.frame(list.DB.train.10R[[k]])
    test.DB <- as.data.frame(list.DB.test.10R[[k]])
    for (partition_data in list(train.DB, test.DB)) {
      if (!all(c("OC", names.used.env) %in% names(partition_data)) ||
          !all(vapply(partition_data[, c("OC", names.used.env), drop = FALSE],
                      is.numeric, logical(1)))) {
        stop("Missing or non-numeric model fields for ", sp.code, ", partition ", k)
      }
      if (any(!is.finite(as.matrix(partition_data[, c("OC", names.used.env), drop = FALSE]))) ||
          !all(partition_data$OC %in% c(0, 1)) || length(unique(partition_data$OC)) != 2L) {
        stop("Each partition needs finite predictors and both OC classes: ", sp.code)
      }
    }
    set.seed(set.num * k.fold.num + i)
    # Numeric OC is retained exactly as in the source. This is not silently
    # replaced by factor-based stratification; inspect saved folds for rare taxa.
    folds <- caret::createFolds(train.DB$OC, k = k.fold.num, returnTrain = FALSE)
    cv_folds[[k]] <- folds
    
    # Bayesian optimization requires a scoring callback. The workflow itself is
    # still executed by the visible species and partition loops above.
    MAX.opt_fun <- function(regmult) {
      fold_tss <- fold_prauc <- numeric(k.fold.num)
      for (ii in seq_len(k.fold.num)) {
        train_folds <- train.DB[-folds[[ii]], , drop = FALSE]
        valid_fold <- train.DB[folds[[ii]], , drop = FALSE]
        if (length(unique(train_folds$OC)) != 2L || length(unique(valid_fold$OC)) != 2L) {
          stop("A tuning fold lacks presence or background records.")
        }
        set.seed(set.num + k * ii)  # The source method index j was always 1.
        p <- train_folds$OC
        predictors <- train_folds[, names.used.env, drop = FALSE]
        target.model <- maxnet::maxnet(
          p, predictors,
          f = maxnet::maxnet.formula(p, predictors, classes = "default"),
          regmult = regmult, addsamplestobackground = TRUE
        )
        pred.raw <- as.numeric(stats::predict(target.model, valid_fold, type = "cloglog"))
        evaluation <- dismo::evaluate(p = pred.raw[valid_fold$OC == 1],
                                      a = pred.raw[valid_fold$OC == 0])
        cutoff <- as.numeric(dismo::threshold(evaluation, stat = "spec_sens"))
        pred.class <- as.integer(pred.raw >= cutoff)
        score_tss <- tryCatch(
          MLmetrics::Sensitivity(y_true = valid_fold$OC, y_pred = pred.class, positive = "1") +
            MLmetrics::Specificity(y_true = valid_fold$OC, y_pred = pred.class, positive = "1") - 1,
          error = function(e) NA_real_
        )
        score_prauc <- tryCatch(
          MLmetrics::PRAUC(y_true = valid_fold$OC, y_pred = pred.raw),
          error = function(e) NA_real_
        )
        fold_tss[ii] <- if (length(score_tss) == 1L && is.finite(score_tss)) score_tss else 0
        fold_prauc[ii] <- if (length(score_prauc) == 1L && is.finite(score_prauc)) score_prauc else 0
      }
      # Retain the source objective: normalize ACROSS FOLDS within each candidate,
      # sum the two normalized metrics, and average. This is not raw mean TSS.
      if (diff(range(fold_tss)) == 0 || diff(range(fold_prauc)) == 0) {
        stop("The source min-max objective is undefined for a constant fold metric.")
      }
      scaled_tss <- (fold_tss - min(fold_tss)) / diff(range(fold_tss))
      scaled_prauc <- (fold_prauc - min(fold_prauc)) / diff(range(fold_prauc))
      list(Score = mean(scaled_tss + scaled_prauc))
    }
    
    cl <- NULL
    optObj <- tryCatch({
      if (optimization_cores > 1L) {
        cl <- parallel::makeCluster(optimization_cores)
        doParallel::registerDoParallel(cl)
        parallel::clusterExport(
          cl, c("train.DB", "folds", "k.fold.num", "names.used.env", "set.num", "k"),
          envir = environment()
        )
        parallel::clusterEvalQ(cl, {library(maxnet); NULL})
      } else {
        foreach::registerDoSEQ()
      }
      ParBayesianOptimization::bayesOpt(
        FUN = MAX.opt_fun, bounds = list(regmult = regmult_bounds),
        initPoints = bayes_initial_points, iters.n = bayes_iterations,
        iters.k = bayes_batch_size, parallel = optimization_cores > 1L,
        plotProgress = FALSE
      )
    }, error = function(e) e,
    finally = {
      if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    })
    
    optimization_error <- NA_character_
    if (inherits(optObj, "error")) {
      optimization_error <- conditionMessage(optObj)
      warning("Optimization failed for ", sp.code, " / partition ", k,
              "; using the source fallback regmult = ", default_regmult,
              ". Reason: ", optimization_error, call. = FALSE)
      selected_regmult <- default_regmult
      optimization_history[[k]] <- data.frame()
      error_flag <- "Error"
    } else {
      selected_regmult <- as.numeric(ParBayesianOptimization::getBestPars(optObj)$regmult)
      optimization_history[[k]] <- as.data.frame(optObj$scoreSummary)
      error_flag <- NA_character_
    }
    if (length(selected_regmult) != 1L || !is.finite(selected_regmult) || selected_regmult <= 0) {
      stop("Invalid regularization multiplier returned for ", sp.code)
    }
    
    fitted_regmult <- selected_regmult
    now.model <- NULL
    last_fit_error <- NA_character_
    for (fit_attempt in seq_len(maximum_fit_attempts)) {
      set.seed(set.num + k)
      fit_result <- tryCatch(
        maxnet::maxnet(train.DB$OC, train.DB[, names.used.env, drop = FALSE],
                       regmult = fitted_regmult, addsamplestobackground = TRUE),
        error = function(e) e
      )
      if (!inherits(fit_result, "error")) {
        now.model <- fit_result
        break
      }
      last_fit_error <- conditionMessage(fit_result)
      if (fit_attempt < maximum_fit_attempts) fitted_regmult <- fitted_regmult + regmult_retry_step
    }
    if (is.null(now.model)) stop("Final model failed for ", sp.code, ": ", last_fit_error)
    list.SDM[[k]] <- now.model
    
    train.pred.raw <- as.numeric(stats::predict(now.model, train.DB, type = "cloglog"))
    test.pred.raw <- as.numeric(stats::predict(now.model, test.DB, type = "cloglog"))
    if (any(!is.finite(c(train.pred.raw, test.pred.raw)))) {
      stop("Non-finite fitted predictions for ", sp.code, ", partition ", k)
    }
    e.model.train <- dismo::evaluate(p = train.pred.raw[train.DB$OC == 1],
                                     a = train.pred.raw[train.DB$OC == 0])
    e.model.test <- dismo::evaluate(p = test.pred.raw[test.DB$OC == 1],
                                    a = test.pred.raw[test.DB$OC == 0])
    t.model.train <- as.numeric(dismo::threshold(e.model.train, stat = "spec_sens"))
    t.model.test <- as.numeric(dismo::threshold(e.model.test, stat = "spec_sens"))
    # Source convention: both train/test binary metrics use the TEST-derived threshold. 
    # They are not an independent evaluation of threshold selection.
    train.pred.OC <- as.integer(train.pred.raw >= t.model.test)
    test.pred.OC <- as.integer(test.pred.raw >= t.model.test)
    train_metrics <- evaluate_performance(train.DB$OC, train.pred.raw, train.pred.OC)
    test_metrics <- evaluate_performance(test.DB$OC, test.pred.raw, test.pred.OC)
    names(train_metrics) <- paste0("train.", names(train_metrics))
    names(test_metrics) <- paste0("test.", names(test_metrics))
    
    partition_info <- data.frame(
      sp = info.sp$common_name, sp.code = sp.code, SDM.method = "MaxEnt", j = 1L, k = k,
      Error.OC = error_flag, N.ori.train = nrow(train.DB),
      N.ori.train.OC0 = sum(train.DB$OC == 0), N.ori.train.OC1 = sum(train.DB$OC == 1),
      N.ori.test = nrow(test.DB), N.ori.test.OC0 = sum(test.DB$OC == 0),
      N.ori.test.OC1 = sum(test.DB$OC == 1),
      Threshold.train = t.model.train, Threshold.test = t.model.test
    )
    evaluation_list[[k]] <- cbind(partition_info, train_metrics, test_metrics)
    parameter_list[[k]] <- data.frame(
      sp = info.sp$common_name, sp.code = sp.code, SDM.method = "MaxEnt", j = 1L, k = k,
      regmult_selected = selected_regmult, regmult = fitted_regmult,
      error.OC = error_flag, optimization_error = optimization_error,
      fit_attempts = fit_attempt, last_fit_error = last_fit_error
    )
    partition_predictions <- rbind(train.DB[, c("OC", names.used.env), drop = FALSE],
                                   test.DB[, c("OC", names.used.env), drop = FALSE])
    partition_predictions$data.type <- c(rep("train", nrow(train.DB)), rep("test", nrow(test.DB)))
    prediction_list[[k]] <- cbind(
      partition_info[, c("sp", "sp.code", "SDM.method", "j", "k", "Threshold.train", "Threshold.test")],
      partition_predictions,
      Pred.raw = c(train.pred.raw, test.pred.raw), Pred.OC = c(train.pred.OC, test.pred.OC)
    )
  }
  
  SDM.evaluation2 <- do.call(rbind, evaluation_list)
  SDM.pred2 <- do.call(rbind, prediction_list)
  SDM.opt.para2 <- do.call(rbind, parameter_list)
  model_settings <- list(seed = set.num, species_index = i, folds = k.fold.num,
                         regmult_bounds = regmult_bounds, initPoints = bayes_initial_points,
                         iters.n = bayes_iterations, iters.k = bayes_batch_size,
                         optimization_cores = optimization_cores, prediction_type = "cloglog",
                         objective = "within_candidate_fold_minmax_TSS_plus_PRAUC",
                         binary_evaluation_threshold = "test_spec_sens", n_repeats = n_repeats)
  save(info.sp, DB.orginal, names.used.env, SDM.formula,
       list.DB.train.10R, list.DB.test.10R, split_row_indices,
       list.SDM, SDM.evaluation2, SDM.pred2, SDM.opt.para2,
       cv_folds, optimization_history, model_settings, file = model_files[i])
  write.csv(SDM.evaluation2, file.path(model_dir, paste0("model_evaluation_", sp.code, ".csv")),
            row.names = FALSE, fileEncoding = output_encoding)
  write.csv(SDM.pred2, file.path(model_dir, paste0("partition_predictions_", sp.code, ".csv")),
            row.names = FALSE, fileEncoding = output_encoding)
  write.csv(SDM.opt.para2, file.path(model_dir, paste0("model_parameters_", sp.code, ".csv")),
            row.names = FALSE, fileEncoding = output_encoding)
}

# ---- 4. Calculate species thresholds from held-out predictions ----------------
# To resume here, execute Sections 0-2 first and use completed model checkpoints.
if (any(!file.exists(model_files))) stop("Some fitted model files are missing.")
if (file.exists(threshold_file) && !overwrite_output) {
  stop("Threshold summary already exists; skip Section 4 or enable overwrite_output.")
}
threshold_summaries <- vector("list", length(model_files))

for (i in seq_along(model_files)) {
  model_data <- new.env(parent = emptyenv())
  load(model_files[i], envir = model_data)
  info.sp <- model_data$info.sp
  sp.code <- as.character(info.sp$fcode)
  names.used.env <- model_data$names.used.env
  SDM.pred2 <- model_data$SDM.pred2
  SDM.evaluation2 <- model_data$SDM.evaluation2
  if (!identical(sp.code, as.character(species_list$fcode[i])) ||
      anyDuplicated(SDM.evaluation2$k) || anyNA(SDM.evaluation2$test.TSS) ||
      any(!is.finite(SDM.evaluation2$test.TSS))) {
    stop("Invalid species code or replicate TSS weights in ", model_files[i])
  }
  message("Calculating thresholds: ", sp.code)
  occurrence_data <- new.env(parent = emptyenv())
  load(species_data_files[i], envir = occurrence_data)
  final.DB.P <- occurrence_data$final.DB.P
  if (is.null(final.DB.P) ||
      !all(c("siteID", "lon", "lat", names.used.env) %in% names(final.DB.P)) ||
      anyNA(final.DB.P[, c("siteID", "lon", "lat", names.used.env), drop = FALSE])) {
    stop("Incomplete processed presence records for ", sp.code)
  }
  
  DB.pred.test <- SDM.pred2[SDM.pred2$data.type == "test", , drop = FALSE]
  # Preserve source grouping by selected environmental values. 
  # row.ID identifies a predictor profile, NOT necessarily a distinct biological record/location.
  row.unique <- unique(DB.pred.test[, names.used.env, drop = FALSE])
  row.unique$row.ID <- seq_len(nrow(row.unique))
  DB.pred.test2 <- merge(DB.pred.test, row.unique, by = names.used.env)
  response_counts <- aggregate(OC ~ row.ID, DB.pred.test2,
                               function(x) length(unique(x)))
  if (any(response_counts$OC != 1L)) {
    write.csv(response_counts[response_counts$OC != 1L, , drop = FALSE],
              file.path(threshold_dir, paste0("conflicting_profiles_", sp.code, ".csv")),
              row.names = FALSE, fileEncoding = output_encoding)
    stop("A predictor profile has both OC labels for ", sp.code,
         ". The source chose the first label; this version stops for review.")
  }
  y_by_id <- aggregate(OC ~ row.ID, DB.pred.test2, function(x) x[1])
  # Raw test TSS weights are intentionally retained here, unlike map weights.
  if (any(SDM.evaluation2$test.TSS < 0)) {
    warning("Negative test TSS retained for source-style threshold weighting: ",
            sp.code, call. = FALSE)
  }
  oof_w_tbl <- merge(DB.pred.test2, SDM.evaluation2[, c("k", "test.TSS")],
                     by = "k", all.x = TRUE)
  if (anyNA(oof_w_tbl$test.TSS)) stop("Unmatched TSS weights for ", sp.code)
  num <- aggregate(oof_w_tbl$Pred.raw * oof_w_tbl$test.TSS,
                   by = list(row.ID = oof_w_tbl$row.ID), FUN = sum, na.rm = TRUE)
  den <- aggregate(oof_w_tbl$test.TSS,
                   by = list(row.ID = oof_w_tbl$row.ID), FUN = sum, na.rm = TRUE)
  names(num)[2] <- "num"
  names(den)[2] <- "den"
  oof_w <- merge(y_by_id, merge(num, den, by = "row.ID"), by = "row.ID")
  oof_w$p_oof_w <- with(oof_w, ifelse(den > 0, num / den, NA_real_))
  # Do not silently substitute uniform/non-negative weights or discard profiles.
  if (any(!is.finite(oof_w$p_oof_w)) || any(oof_w$p_oof_w < 0 | oof_w$p_oof_w > 1)) {
    write.csv(oof_w, file.path(threshold_dir, paste0("invalid_weighted_predictions_", sp.code, ".csv")),
              row.names = FALSE, fileEncoding = output_encoding)
    stop("Undefined or out-of-range source-weighted predictions for ", sp.code,
         ". Review raw TSS weights before choosing a revised weighting rule.")
  }
  tev <- dismo::evaluate(p = oof_w$p_oof_w[oof_w$OC == 1],
                         a = oof_w$p_oof_w[oof_w$OC == 0])
  final_cutoff <- as.numeric(dismo::threshold(tev, stat = "spec_sens"))
  
  # Original ROI matching: domestic processed presences joined to ALL train/test
  # predictions by selected predictors, then weighted within siteID. 
  # This is not an independent held-out domestic evaluation or an administrative-border test.
  DB.korea <- final.DB.P[
    final.DB.P$lon > roi_bounds["xmin"] & final.DB.P$lon < roi_bounds["xmax"] &
      final.DB.P$lat > roi_bounds["ymin"] & final.DB.P$lat < roi_bounds["ymax"], , drop = FALSE]
  DB.K.thr <- numeric(0)
  final_cutoff2 <- final_cutoff
  roi_status <- "no_domestic_records"
  roi_candidates <- data.frame(threshold = numeric(), coverage = numeric(),
                               selection_score = numeric())
  if (nrow(DB.korea)) {
    DB.k2 <- merge(DB.korea, SDM.pred2[, c(names.used.env, "k", "Pred.raw")],
                   by = names.used.env)
    DB.k3 <- merge(DB.k2, SDM.evaluation2[, c("k", "test.TSS")], by = "k")
    if (!nrow(DB.k3) || !all(unique(DB.korea$siteID) %in% DB.k3$siteID)) {
      stop("Domestic presence records could not all be matched for ", sp.code)
    }
    DB.K.thr <- vapply(split(DB.k3, DB.k3$siteID), function(x) {
      sum(x$Pred.raw * x$test.TSS) / sum(x$test.TSS)
    }, numeric(1))
    if (any(!is.finite(DB.K.thr)) || any(DB.K.thr < 0 | DB.K.thr > 1)) {
      stop("Invalid TSS-weighted domestic predictions for ", sp.code)
    }
    roi_status <- "unadjusted"
    if (!any(DB.K.thr >= final_cutoff)) {
      lo <- min(DB.K.thr)
      hi <- max(DB.K.thr)
      if (lo == hi) {
        final_cutoff2 <- lo
        roi_status <- "adjusted_constant_domestic_prediction"
      } else {
        candidates <- sort(unique(oof_w$p_oof_w[oof_w$p_oof_w >= lo & oof_w$p_oof_w <= hi]))
        if (length(candidates) < 100L) {
          candidates <- sort(unique(c(candidates, seq(lo, hi, length.out = 200L))))
        }
        coverage <- vapply(candidates, function(cutoff) mean(DB.K.thr >= cutoff), numeric(1))
        scores <- rep(NA_real_, length(candidates))
        for (candidate_index in which(coverage >= roi_minimum_coverage)) {
          yhat <- as.integer(oof_w$p_oof_w >= candidates[candidate_index])
          # Explicitly distinguish literal-source reproduction from a correction.
          if (roi_score_convention == "source_argument_order") {
            y_reference <- yhat
            y_prediction <- oof_w$OC
          } else {
            y_reference <- oof_w$OC
            y_prediction <- yhat
          }
          score <- MLmetrics::Sensitivity(y_true = y_reference, y_pred = y_prediction, positive = "1") +
            MLmetrics::Specificity(y_true = y_reference, y_pred = y_prediction, positive = "1") - 1
          if (length(score) == 1L && is.finite(score)) scores[candidate_index] <- score
        }
        if (!any(is.finite(scores))) {
          stop("All ROI candidate scores are undefined for ", sp.code,
               "; inspect the scoring convention rather than assigning an arbitrary threshold.")
        }
        final_cutoff2 <- candidates[which.max(scores)]
        roi_status <- "adjusted_within_domestic_range"
        roi_candidates <- data.frame(threshold = candidates, coverage = coverage,
                                     selection_score = scores)
      }
    }
  }
  
  DF.final1 <- evaluate_performance(oof_w$OC, oof_w$p_oof_w,
                                    as.integer(oof_w$p_oof_w >= final_cutoff))
  DF.final2 <- evaluate_performance(oof_w$OC, oof_w$p_oof_w,
                                    as.integer(oof_w$p_oof_w >= final_cutoff2))
  names(DF.final1) <- paste0("Ens.", names(DF.final1))
  names(DF.final2) <- paste0("ROI.", names(DF.final2))
  threshold_summaries[[i]] <- cbind(
    info.sp,
    data.frame(thr1.Ensemble = final_cutoff, thr2.ROI = final_cutoff2,
               roi_status = roi_status, roi_score_convention = roi_score_convention,
               n_oof_predictor_profiles = nrow(oof_w), n_domestic_records = nrow(DB.korea),
               n_domestic_sites = length(DB.K.thr),
               domestic_coverage_before = if (length(DB.K.thr)) mean(DB.K.thr >= final_cutoff) else NA_real_,
               domestic_coverage_after = if (length(DB.K.thr)) mean(DB.K.thr >= final_cutoff2) else NA_real_),
    DF.final1, DF.final2
  )
  save(info.sp, names.used.env, row.unique, oof_w, DB.K.thr,
       final_cutoff, final_cutoff2, roi_status, roi_score_convention, roi_candidates,
       file = file.path(threshold_dir, paste0("threshold_details_", sp.code, ".RData")))
  write.csv(oof_w, file.path(threshold_dir, paste0("weighted_test_predictions_", sp.code, ".csv")),
            row.names = FALSE, fileEncoding = output_encoding)
}
Info.Fishes <- do.call(rbind, threshold_summaries)
write.csv(Info.Fishes, threshold_file, row.names = FALSE, fileEncoding = output_encoding)

# ---- 5. Prepare spatial inputs and the raster-prediction callback -------------
# For mapping-only resumption, execute Sections 0-1 and then Sections 5-7.
# Existing fitted models and the completed threshold summary are required.
required_spatial_files <- c(elevation_file, boundary_file, water_mask_file, threshold_file)
if (any(!file.exists(c(required_spatial_files, model_files)))) {
  stop("Fitted models, thresholds, elevation, boundary, and water mask are required for mapping.")
}
if (!length(scenario_rows) || anyNA(scenario_rows) || anyDuplicated(scenario_rows) ||
    !all(scenario_rows %in% seq_len(nrow(scenarios)))) {
  stop("scenario_rows must select unique existing rows of scenarios.")
}
needed_climate_dirs <- climate_dirs[unique(scenarios$period[scenario_rows])]
if (any(!dir.exists(needed_climate_dirs))) stop("One or more selected climate directories are missing.")

f.thr <- read.csv(threshold_file, fileEncoding = output_encoding, check.names = FALSE,
                  stringsAsFactors = FALSE)
if (!all(c("fcode", "thr2.ROI") %in% names(f.thr)) || anyDuplicated(f.thr$fcode)) {
  stop("Threshold table requires one fcode and thr2.ROI per species.")
}
Env.ele <- terra::rast(elevation_file)
if (terra::nlyr(Env.ele) != 1L) stop("The elevation file must contain one layer.")
names(Env.ele) <- "elev"
kor.boundary <- terra::vect(boundary_file)
mask.water <- terra::rast(water_mask_file)
if (terra::nlyr(mask.water) != 1L) stop("The inland-water mask must contain one layer.")
water_values <- terra::freq(mask.water)
if (is.null(water_values) || !nrow(water_values)) stop("The inland-water mask has no valid cells.")
if (any(water_values$value == 0, na.rm = TRUE)) {
  stop("The source polygon mask includes every non-NA cell. Supply a water-only mask",
       " with non-water cells set to NA; do not supply an unexplained 0/1 mask.")
}
mask.wv <- terra::as.polygons(mask.water, aggregate = TRUE, na.rm = TRUE)
# No basin-specific dispersal restriction is applied in these three source files.
# Such restrictions belong to the subsequent spatial-analysis script.

predict_maxnet_cloglog <- function(model, newdata, ...) {
  as.numeric(stats::predict(model, newdata, type = "cloglog", ...))
}

# ---- 6. Generate current and future ensemble prediction rasters ---------------
map_manifest <- list()
ensemble_weight_tables <- list()
prediction_grid <- NULL
for (cmt in scenario_rows) {
  scenario <- scenarios[cmt, , drop = FALSE]
  scenario_name <- paste(scenario$period, scenario$sce, sep = "_")
  message("Predicting climate scenario: ", scenario_name)
  climate_files <- list.files(climate_dirs[[scenario$period]],
                              pattern = "\\.tif(f)?$", full.names = TRUE, ignore.case = TRUE)
  climate_files <- climate_files[grepl("bio[0-9]+", basename(climate_files), ignore.case = TRUE)]
  if (scenario$sce != "normal") {
    climate_files <- climate_files[grepl(scenario$sce, basename(climate_files), ignore.case = TRUE)]
  }
  # Require one file for each predictor. No GCM averaging, resampling, or silent
  # choice among alternative climate products is performed.
  climate_labels <- sub(".*?(bio[0-9]+).*", "\\1", tolower(basename(climate_files)), perl = TRUE)
  climate_numbers <- suppressWarnings(as.integer(sub("bio", "", climate_labels, fixed = TRUE)))
  if (length(climate_files) != 19L || anyNA(climate_numbers) ||
      !setequal(climate_numbers, seq_len(19L)) || anyDuplicated(climate_numbers)) {
    stop("Need exactly bio1-bio19 for ", scenario_name,
         ". Keep only one selected climate product/GCM per scenario and period.")
  }
  file_order <- order(climate_numbers)
  Env.CHELSA <- terra::rast(climate_files[file_order])
  if (terra::nlyr(Env.CHELSA) != 19L) stop("Climate inputs must be single-layer GeoTIFFs.")
  names(Env.CHELSA) <- paste0("bio", seq_len(19L))
  Env.ele.c <- terra::crop(Env.ele, Env.CHELSA)
  terra::compareGeom(Env.ele.c, Env.CHELSA[[1]], stopOnError = TRUE)
  Env.maps <- c(Env.ele.c, Env.CHELSA)
  boundary_use <- terra::project(kor.boundary, terra::crs(Env.maps))
  water_use <- terra::project(mask.wv, terra::crs(Env.maps))
  tmap <- terra::crop(Env.maps, boundary_use)
  tmap <- terra::mask(tmap, boundary_use)
  Env.maps.kor <- terra::mask(tmap, tmap[["elev"]])
  Env.maps.water <- terra::mask(Env.maps.kor, water_use)
  if (is.null(prediction_grid)) {
    prediction_grid <- terra::rast(Env.maps.water[[1]])
  } else {
    terra::compareGeom(prediction_grid, Env.maps.water[[1]], stopOnError = TRUE)
  }
  rm(Env.CHELSA, Env.ele.c, Env.maps, tmap, Env.maps.kor)
  scenario_output_dir <- file.path(map_dir, scenario_name)
  dir.create(scenario_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (ii in seq_along(model_files)) {
    model_data <- new.env(parent = emptyenv())
    load(model_files[ii], envir = model_data)
    info.sp <- model_data$info.sp
    sp.code <- as.character(info.sp$fcode)
    list.SDM <- model_data$list.SDM
    t.ev <- model_data$SDM.evaluation2
    names.used.env <- model_data$names.used.env
    if (!identical(sp.code, as.character(species_list$fcode[ii])) ||
        !all(names.used.env %in% names(Env.maps.water)) ||
        anyDuplicated(t.ev$k) || !setequal(t.ev$k, seq_along(list.SDM))) {
      stop("Model, predictor, or replicate mismatch for ", sp.code)
    }
    # Align weights explicitly with model-list order; keep ALL fitted replicates.
    model.select.num <- seq_along(list.SDM)
    t.ev <- t.ev[match(model.select.num, t.ev$k), , drop = FALSE]
    w <- t.ev$test.TSS
    if (any(is.infinite(w))) stop("Infinite map weights for ", sp.code)
    w[is.na(w) | w < 0] <- 0
    uniform_fallback <- sum(w) == 0
    w <- if (uniform_fallback) rep(1 / length(w), length(w)) else w / sum(w)
    final_cutoff <- f.thr$thr2.ROI[match(sp.code, f.thr$fcode)]
    if (length(final_cutoff) != 1L || !is.finite(final_cutoff)) {
      stop("Missing or invalid mapping threshold for ", sp.code)
    }
    continuous_file <- file.path(scenario_output_dir,
                                 paste0("suitability_", sp.code, "_MaxEnt_", scenario_name, ".tif"))
    binary_file <- file.path(scenario_output_dir,
                             paste0("binary_", sp.code, "_MaxEnt_", scenario_name, ".tif"))
    sd_file <- file.path(scenario_output_dir,
                         paste0("prediction_sd_", sp.code, "_MaxEnt_", scenario_name, ".tif"))
    if (any(file.exists(c(continuous_file, binary_file, sd_file))) && !overwrite_output) {
      stop("Prediction rasters already exist for ", sp.code, " / ", scenario_name)
    }
    predicted.map.list <- vector("list", length(list.SDM))
    for (k in model.select.num) {
      # All 20 input layers are retained in the NA mask, as in the source.
      predicted.map.list[[k]] <- terra::predict(
        Env.maps.water, list.SDM[[k]], fun = predict_maxnet_cloglog,
        na.rm = TRUE, cores = prediction_cores, cpkgs = "maxnet"
      )
    }
    stk <- terra::rast(predicted.map.list)
    ensemble_prob <- terra::app(stk, fun = function(x, weights) {
      if (is.null(dim(x))) sum(x * weights) else as.numeric(x %*% weights)
    }, weights = w)
    uncertainty_sd <- terra::app(stk, fun = "sd")  # Unweighted sample SD; not a CI.
    binary_map <- (ensemble_prob >= final_cutoff) * 1
    names(ensemble_prob) <- "suitability"
    names(binary_map) <- "suitable_habitat"
    names(uncertainty_sd) <- "prediction_sd"
    terra::writeRaster(ensemble_prob, continuous_file, overwrite = overwrite_output)
    terra::writeRaster(binary_map, binary_file, overwrite = overwrite_output)
    terra::writeRaster(uncertainty_sd, sd_file, overwrite = overwrite_output)
    map_manifest[[length(map_manifest) + 1L]] <- data.frame(
      fcode = sp.code, period = scenario$period, scenario = scenario$sce,
      n_models = length(list.SDM), threshold = final_cutoff,
      suitability_file = continuous_file, binary_file = binary_file, sd_file = sd_file
    )
    ensemble_weight_tables[[length(ensemble_weight_tables) + 1L]] <- data.frame(
      fcode = sp.code, climate = scenario_name, k = model.select.num,
      test.TSS = t.ev$test.TSS, ensemble_weight = w, uniform_fallback = uniform_fallback
    )
    rm(predicted.map.list, stk, ensemble_prob, binary_map, uncertainty_sd, model_data)
  }
  rm(Env.maps.water)
}
# These manifests describe the scenarios selected in the current execution.
write.csv(do.call(rbind, map_manifest), file.path(output_dir, "prediction_map_manifest.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
write.csv(do.call(rbind, ensemble_weight_tables), file.path(output_dir, "ensemble_weights.csv"),
          row.names = FALSE, fileEncoding = output_encoding)

# ---- 7. Record the computational environment ---------------------------------
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
package_versions <- data.frame(
  package = required_packages,
  version = vapply(required_packages, function(x) as.character(utils::packageVersion(x)), character(1))
)
write.csv(package_versions, file.path(output_dir, "package_versions.csv"),
          row.names = FALSE, fileEncoding = output_encoding)
message("Species distribution modelling complete. Outputs: ", output_dir)
