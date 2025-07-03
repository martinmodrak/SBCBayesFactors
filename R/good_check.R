good_check <- function(probs, true_model, prior_prob1 = 0.5) {
  bf01s <- ((1 - probs) / probs) / ((1 - prior_prob1) / prior_prob1)


  bf01s_given_h1 <- bf01s[true_model == 1]
  bf10s_given_h0 <- 1/bf01s[true_model == 0]

  list(
    "01" = t.test(bf01s_given_h1, mu = 1),
    "10" = t.test(bf10s_given_h0, mu = 1),
    "combined" = t.test(c(bf01s_given_h1, bf10s_given_h0), mu = 1)
  )
}

good_check_okada <- function(probs, true_model, prior_prob1 = 0.5) {
  bf01s <- ((1 - probs) / probs) / ((1 - prior_prob1) / prior_prob1)

  bf01s_frac_given_h0 <- bf01s[true_model == 0]^-0.5
  bf01s_frac_given_h1 <- bf01s[true_model == 1]^0.5

  t.test(bf01s_frac_given_h0, bf01s_frac_given_h1)
}


run_good_sims <- function(n_sims, n_runs, m0, m1 = 0, step = 10, start = 1,
                          simulate_from = 1, numerator = 0, moment = 1) {
  if(simulate_from == 1) {
    rng_def <- m1
  } else if(simulate_from == 0) {
    rng_def <- m0
  } else {
    stop("simulate_from")
  }

  if(rng_def == "cauchy") {
      rng <- rcauchy
  } else {
      rng <- \(...) {rnorm(..., mean = rng_def)}
  }


  y <- matrix(rng(n_sims * n_runs), nrow = n_sims)

  if(m0 == "cauchy") {
    ll0 <- dcauchy(y, log = TRUE)
  } else {
    ll0 <- dnorm(y, mean = m0, log = TRUE)
  }
  if(m1 == "cauchy") {
    ll1 <- dcauchy(y, log = TRUE)
  } else {
    ll1 <- dnorm(y, mean = m1, log = TRUE)
  }
  if(numerator == 0) {
    bf <- exp(ll0 - ll1)
  } else if (numerator == 1) {
    bf <- exp(ll1 - ll0)
  }

  bf_moment <- bf^moment
  cmeans <- apply(bf_moment, MARGIN = 2, FUN = \(x) cumsum(x) / (1:length(x)))
  colnames(cmeans) <- 1:n_runs
  sims <- seq(start, n_sims, by = step)
  cmeans <- cmeans[sims,]
  rownames(cmeans) <- sims
  attr(cmeans, "raw_bf_moment") <- bf_moment

  cmeans
}

plot_good_example <- function(sims_matrix,
                              linealpha = 0.3,
                              correct_value = 1,
                              correct_range = 0.1,
                              ylim = NULL) {

  sims <- as.integer(rownames(sims_matrix))
  plot_df <- as.data.frame(sims_matrix) |> mutate(sim = sims) |>
    tidyr::pivot_longer(c(everything(), - all_of("sim")), names_to = "run", values_to = "mean_bf")

  if(!is.null(correct_value)) {
    correct_geom <- geom_hline(yintercept = correct_value, color = "darkblue")
  } else {
    correct_geom <- NULL
  }

  if(!is.null(correct_value) && !is.null(correct_range)) {
    correct_rect_geom <- geom_rect(xmin = min(sims), xmax = max(sims), ymin = correct_value - correct_range, ymax = correct_value + correct_range, fill = "lightblue")
    expand <- expand_limits(y = c(correct_value - correct_range, correct_value + correct_range))
  } else {
    correct_rect_geom <- NULL
    expand <- NULL
  }

  plot_df |> ggplot(aes(x = sim, y = mean_bf, group = run)) +
    correct_rect_geom +
    correct_geom +
    geom_line(alpha = linealpha) +
    scale_y_continuous("Cumulative mean BF", limits = ylim) +
    #scale_y_log10("Cumulative mean BF") +
    scale_x_continuous("Simulations") +
    expand
}
