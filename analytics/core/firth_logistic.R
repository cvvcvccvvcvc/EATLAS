suppressPackageStartupMessages(library(logistf))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("usage: firth_logistic.R <cohort.tsv> <specs.tsv> <results.tsv> <versions.tsv>")
}

cohort <- read.delim(args[[1]], check.names = FALSE, stringsAsFactors = FALSE)
specs <- read.delim(args[[2]], check.names = FALSE, stringsAsFactors = FALSE)

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
  tryCatch({
    fit <- logistf(
      benign ~ ALT_observed + splines::ns(score, df = 3),
      data = model_data,
      firth = TRUE,
      pl = TRUE,
      plconf = 2
    )
    coefficient_index <- match("ALT_observed", names(fit$coefficients))
    data.frame(
      analysis_id = spec$analysis_id,
      odds_ratio = exp(fit$coefficients[[coefficient_index]]),
      ci_low = exp(fit$ci.lower[[coefficient_index]]),
      ci_high = exp(fit$ci.upper[[coefficient_index]]),
      plr_p = fit$prob[[coefficient_index]],
      status = "estimated",
      reason = "",
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
      reason = paste0("Firth model failed: ", conditionMessage(error)),
      stringsAsFactors = FALSE
    )
  })
}

results <- do.call(rbind, lapply(seq_len(nrow(specs)), function(index) fit_one(specs[index, ])))
write.table(results, args[[3]], sep = "\t", row.names = FALSE, quote = FALSE, na = "")

versions <- data.frame(
  component = c("R", "logistf"),
  version = c(as.character(getRversion()), as.character(packageVersion("logistf"))),
  stringsAsFactors = FALSE
)
write.table(versions, args[[4]], sep = "\t", row.names = FALSE, quote = FALSE)
