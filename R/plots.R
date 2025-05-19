my_reliability_diag <- function(binary_probs, region.method = "resampling") {
  rel_diag <- reliabilitydiag::reliabilitydiag(
    x = binary_probs$prob,
    y = binary_probs$simulated_value,
    region.level = 0.95,
    region.method = region.method,
    n.boot = 1000
  )



  res <- autoplot(
    rel_diag,
    params_CEPsegment = list(linewidth = 0.5, color = "red")
  )

  # attr(res, "max_diff") <- rel_diag$x$regions |>
  #   mutate(max_diff = pmax(abs(x - lower), abs(x - upper))) |>
  #   pull(max_diff) |>
  #   max()

  res
}
