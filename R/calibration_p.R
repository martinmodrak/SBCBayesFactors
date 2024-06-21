calibration_p <- function(x, y, adjust.method = "holm") {
  ys <- split(y, x)
  ns <- lengths(ys)
  ys <- sapply(ys, sum)
  probs <- as.numeric(names(ys))

  n <- length(ys)
  all_ps <- rep(NA_real_, n*(n + 1))
  p_index <- 1
  for(j in 1:n) {
    csum_ys <- cumsum(ys[j:n])
    csum_ns <- cumsum(ns[j:n])
    p_low <- pbinom(csum_ys, size = csum_ns, prob = probs[j])
    p_high <- pbinom(csum_ys - 1, size = csum_ns, prob = probs[j:n], lower.tail = FALSE)

    n_new_ps <- n - j + 1
    all_ps[p_index:(p_index + n_new_ps - 1)] <- p_low
    p_index <- p_index + n_new_ps
    all_ps[p_index:(p_index + n_new_ps - 1)] <- p_high
    p_index <- p_index + n_new_ps
  }
  stopifnot(all(!is.na(all_ps)))
  stopifnot(length(all_ps) == n * (n + 1))
  if(adjust.method == "my") {
    min(1, min(all_ps) * n * (n + 1))
  } else {
    min(p.adjust(all_ps, method = adjust.method))
  }
}
