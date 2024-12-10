my_reliability_diag <- function(binary_probs) {
  plot(reliabilitydiag::reliabilitydiag(x = binary_probs$prob, y = binary_probs$simulated_value, region.level = 0.95, region.method = "resampling", n.boot = 1000))
}
