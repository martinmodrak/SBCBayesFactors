my_reliability_diag <- function(binary_probs, region.method = "resampling") {
  autoplot(
    reliabilitydiag::reliabilitydiag(
      x = binary_probs$prob,
      y = binary_probs$simulated_value,
      region.level = 0.95,
      region.method = region.method,
      n.boot = 1000
    ),
    params_CEPsegment = list(linewidth = 0.5, color = "red")
  )
}
