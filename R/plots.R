my_reliability_diag <- function(binary_probs, region.method = "resampling") {
  rel_diag <- reliabilitydiag::reliabilitydiag(
    x = binary_probs$prob,
    y = binary_probs$simulated_value,
    region.level = 0.95,
    region.method = region.method,
    n.boot = 1000
  )



  res <- ggplot2::autoplot(
    rel_diag,
    params_CEPsegment = list(linewidth = 0.5, color = "red")
  )

  # attr(res, "max_diff") <- rel_diag$x$regions |>
  #   mutate(max_diff = pmax(abs(x - lower), abs(x - upper))) |>
  #   pull(max_diff) |>
  #   max()

  res
}

compare_binary_probs <- function(bp1, bp2, scale = c("prob", "BF"), lab1 = waiver(), lab2 = waiver(), title = NULL, reference = NULL) {
  stopifnot(identical(bp1$sim_id, bp2$sim_id))
  stopifnot(identical(bp1$variable, bp2$variable))
  stopifnot(identical(bp1$simulated_value, bp2$simulated_value))

  bp_comparison <- bp1 |>
    select(sim_id, variable, simulated_value, prob) |>
    rename(prob1 = prob) |>
    mutate(prob2 = bp2$prob)

  scale <- match.arg(scale)
  if(scale == "prob") {
    base_plot <-
    bp_comparison |> ggplot() + aes(x = prob1, y = prob2) +#, color = after_stat(ndensity))+
      scale_x_continuous(lab1) +
      scale_y_continuous(lab2)
  } else {
    stopifnot(!is.null(reference))
    stopifnot(reference %in% bp1$variable)
    base_plot <- bp_comparison |> group_by(sim_id) |>
      mutate(BF1 = pmax(1e-15, prob1 / prob1[variable == reference]),
             BF2 = pmax(1e-15, prob2 / prob2[variable == reference])) |>
      ungroup() |>
      filter(variable != reference) |>
      ggplot() + aes(x = BF1, y = BF2) + #, color = after_stat(ndensity))+
      scale_x_log10(lab1) +
      scale_y_log10(lab2)
  }
  base_plot  + geom_abline(color = "orangered", intercept = 0, slope = 1) +
    #ggpointdensity::geom_pointdensity(adjust = 4)  +
    geom_point(alpha = 0.03) +
    facet_wrap(~variable, nrow = 1) + coord_fixed() + scale_color_viridis_c("relative density") +  ggtitle(title)

}
