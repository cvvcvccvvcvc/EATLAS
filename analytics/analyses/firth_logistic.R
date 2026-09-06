suppressPackageStartupMessages(library(logistf))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("usage: firth_logistic.R <cohort.tsv> <specs.tsv> <results.tsv> <versions.tsv> <workers>")
}

cohort <- read.delim(args[[1]], check.names = FALSE, stringsAsFactors = FALSE)
specs <- read.delim(args[[2]], check.names = FALSE, stringsAsFactors = FALSE)
workers <- suppressWarnings(as.integer(args[[5]]))
if (is.na(workers) || workers < 1) {
  stop("workers must be a positive integer")
}
workers <- min(workers, nrow(specs))

matches_variant_type <- function(values, selector) {
  if (selector == "all") return(rep(TRUE, length(values)))
  if (selector == "indel") return(values %in% c("insertion", "deletion"))
  values == selector
}

matches_consequence <- function(values, selector) {
  if (selector == "all") return(rep(TRUE, length(values)))
  pattern <- paste0("(^|\\|)", selector, "(\\||$)")
  grepl(pattern, values)
}

matches_target_context <- function(values, selector) {
  if (selector == "all") return(rep(TRUE, length(values)))
  values == selector
}

fit_one <- function(spec) {
  keep <- cohort[[spec$eligibility_column]] == 1 &
    matches_variant_type(cohort$variant_subtype, spec$variant_type) &
    matches_target_context(cohort$target_context, spec$target_context) &
    matches_consequence(cohort$consequence_groups, spec$consequence)
  model_data <- data.frame(
    benign = cohort$benign[keep],
    ALT_observed = cohort[[spec$observation_column]][keep],
    score = cohort$score[keep]
  )
  model_warnings <- character()
  clean_message <- function(value) gsub("[\t\r\n]+", " ", value)
  tryCatch({
    fit <- withCallingHandlers(logistf(
      benign ~ ALT_observed + splines::ns(score, df = 3),
      data = model_data,
      firth = TRUE,
      pl = TRUE,
      plconf = 2
    ), warning = function(warning) {
      model_warnings <<- c(model_warnings, clean_message(conditionMessage(warning)))
      invokeRestart("muffleWarning")
    })
    coefficient_index <- match("ALT_observed", names(fit$coefficients))
    values <- c(
      exp(fit$coefficients[[coefficient_index]]),
      exp(fit$ci.lower[[coefficient_index]]),
      exp(fit$ci.upper[[coefficient_index]]),
      fit$prob[[coefficient_index]]
    )
    nonconverged <- any(grepl(
      "not converg|non.?converg|maximum number of iterations",
      model_warnings,
      ignore.case = TRUE
    ))
    invalid <- any(!is.finite(values)) || any(values[1:3] <= 0) ||
      values[[2]] > values[[3]] || values[[4]] < 0 || values[[4]] > 1
    if (invalid) {
      model_warnings <- c(
        model_warnings,
        "Firth returned invalid effect, confidence limits, or p-value."
      )
    }
    usable <- !nonconverged && !invalid
    data.frame(
      analysis_id = spec$analysis_id,
      odds_ratio = if (usable) values[[1]] else NA_real_,
      ci_low = if (usable) values[[2]] else NA_real_,
      ci_high = if (usable) values[[3]] else NA_real_,
      plr_p = if (usable) values[[4]] else NA_real_,
      status = if (!usable) {
        "not_estimable"
      } else if (length(model_warnings)) {
        "estimated_warning"
      } else {
        "estimated"
      },
      reason = paste(unique(model_warnings), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }, error = function(error) {
    data.frame(
      analysis_id = spec$analysis_id,
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      plr_p = NA_real_,
      status = "not_estimable",
      reason = paste(
        c(
          unique(model_warnings),
          paste0("Firth model failed: ", clean_message(conditionMessage(error)))
        ),
        collapse = "; "
      ),
      stringsAsFactors = FALSE
    )
  })
}

fit_indices <- seq_len(nrow(specs))
if (workers > 1 && .Platform$OS.type != "windows") {
  fitted <- parallel::mclapply(
    fit_indices,
    function(index) fit_one(specs[index, ]),
    mc.cores = workers,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  fitted <- lapply(fit_indices, function(index) fit_one(specs[index, ]))
}
results <- do.call(rbind, fitted)
write.table(results, args[[3]], sep = "\t", row.names = FALSE, quote = FALSE, na = "")

versions <- data.frame(
  component = c("R", "logistf"),
  version = c(as.character(getRversion()), as.character(packageVersion("logistf"))),
  stringsAsFactors = FALSE
)
write.table(versions, args[[4]], sep = "\t", row.names = FALSE, quote = FALSE)
