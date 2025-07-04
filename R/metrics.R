calibration_metrics <- function(res, prob1_prior = 0.5) {
  bp <- binary_probabilities_from_stats(res$stats)
  if(length(prob1_prior) == 1) {
    t_res <- t.test(bp$prob, mu = prob1_prior)
  } else if(all(prob1_prior %in% c(0,1))) {
    stopifnot(length(prob1_prior) == length(res))
    t_res <- t.test(bp$prob, prob1_prior)
  } else {
    stop("Invalid prob1_prior")
  }
  miscalibration_stats <- miscalibration_resampling_stats(bp$prob, bp$simulated_value)
  reliability_diag <- my_reliability_diag(bp)

  max_rank <- unique(res$stats$max_rank)
  stopifnot(length(max_rank) == 1)

  n_sims <- length(res)
  log_gammas <- res$stats |> group_by(variable) |>
    filter(all(!is.na(rank))) |>
    summarise(log_gamma = log_gamma_stat(rank, !!max_rank))


  gamma_limit <- SBC:::adjust_gamma(N = n_sims, L = 1, K = max_rank + 1)
  max_ecdf_diff <- max(abs(qbinom(c(gamma_limit / 2, 1 - gamma_limit), size = n_sims, prob = 0.5) / n_sims - 0.5))

  structure(
    list(t  = t_res,
         n_sims = n_sims,
         miscalibration_stats = miscalibration_stats,
         log_gammas = log_gammas,
         max_ecdf_diff = max_ecdf_diff,
         log_gamma_limit = log(gamma_limit),
         prob1_prior = prob1_prior,
         reliability_diag = reliability_diag),
    class = "calibration_metrics"
  )
}

#' @export
print.calibration_metrics <- function(m) {
  cat(m$n_sims, " simulations\n")
  cat("DAP - ", m$t$method, " - H0: ", m$t$null.value, "\n95% CI for difference: [", m$t$conf.int[1] - m$t$null.value, ", ", m$t$conf.int[2] - m$t$null.value,"], p = ", m$t$p.value,  "\n", sep = "")
  print(m$miscalibration_stats)
  cat("Log gamma range: [", min(m$log_gammas$log_gamma), ", ", max(m$log_gammas$log_gamma), "]\n",
      "Q95% under null: ", m$log_gamma_limit, ", ", sum(m$log_gammas$log_gamma < m$log_gamma_limit), "/", nrow(m$log_gammas), " (", scales::percent(mean(m$log_gammas$log_gamma < m$log_gamma_limit)), ") below threshold.\n", sep = "")
  cat("Sensitive to eCDF difference up to ", m$max_ecdf_diff, "\n")
  print(m$reliability_diag)
}


report_success_metrics <- function(m, dap_digits = 3, miscalib_digits = 4, ecdf_diff_digits = 3, tex = FALSE) {
  my_format <- function(x, digits) {
    format(round(x, digits), scientific = FALSE)
  }
  text <- paste0("95% CI for DAP difference from prior: ", my_format(m$t$conf.int[1] - m$t$null.value, dap_digits), " -- ", my_format(m$t$conf.int[2] - m$t$null.value, dap_digits),"; ",
      "miscalibration: ", my_format(m$miscalibration_stats$observed, miscalib_digits), ", ", scales::percent(1 - m$miscalibration_stats$alpha), " quantile under null: ", my_format(m$miscalibration_stats$rejection_limit, miscalib_digits),"; ",
      "SBC sensitive to eCDF difference up to ", my_format(m$max_ecdf_diff, ecdf_diff_digits),
      sep = "")

  if(tex) {
    text <- gsub("%", "\\%", text, fixed = TRUE)
  }
  cat(text)
}

success_metrics_for_table <- function(m, dap_digits = 4, miscalib_digits = 4, ecdf_diff_digits = 3) {
  my_format <- function(x, digits) {
    format(round(x, digits), scientific = FALSE)
  }
  stopifnot(m$miscalibration_stats$alpha == 0.05)
  data.frame(check.names = FALSE,
    "#Sims" = m$n_sims,
    "DAP 95% CI" = paste0(my_format(m$t$conf.int[1] - m$t$null.value, dap_digits), " -- ", my_format(m$t$conf.int[2] - m$t$null.value, dap_digits)),
    "Miscalibration" = my_format(m$miscalibration_stats$observed, miscalib_digits),
    "Miscalibration Q95%" = my_format(m$miscalibration_stats$rejection_limit, miscalib_digits),
    "SBC sensitivivity" = my_format(m$max_ecdf_diff, ecdf_diff_digits)
  ) |> tibble::remove_rownames()
}
