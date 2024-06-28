#' Calls sinib::psinib, but handles some extra edge cases
my_psinib <- function(q, size, prob, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(length(size) == length(prob))
  stopifnot(all(prob >= 0) & all(prob <= 1))

  size <- size[prob > 0]
  prob <- prob[prob > 0]

  q <- q - sum(size[prob == 1])
  size <- size[prob < 1]
  prob <- prob[prob < 1]
  if(length(size) == 1) {
    return(pbinom(q, size, prob, lower.tail = lower.tail, log.p = log.p))
  } else {
    return(sinib::psinib(q, size, prob, lower.tail = lower.tail, log.p = log.p))
  }
}

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
  min(p.adjust(all_ps, method = adjust.method))
}

calibration_p_saddle <- function(x, y, adjust.method = "holm", round = TRUE) {
  #cat(paste0("x <- c(", paste(x, collapse = ", "), "); y <- c(", paste(y, collapse = ", "), ");\n"))
  if(round) {
    digits <- 2
    d10 <- 10^digits

    # TODO: this is opposite the original code where
    # they floor the upr, check which is correct
    x_upr <- ceiling(x * d10)/d10
    ys_upr <- split(y, x_upr)
    probs_upr <- as.numeric(names(ys_upr))
    ns_upr <- as.integer(lengths(ys_upr))
    ys_upr <- as.integer(sapply(ys_upr, sum))

    x_lwr <- floor(x * d10)/d10
    ys_lwr <- split(y, x_lwr)
    probs_lwr <- as.numeric(names(ys_lwr))
    ns_lwr <- as.integer(lengths(ys_lwr))
    ys_lwr <- as.integer(sapply(ys_lwr, sum))

  } else {
    ys <- split(y, x)
    probs <- as.numeric(names(ys))
    ns <- as.integer(lengths(ys))
    ys <- as.integer(sapply(ys, sum))

    ys_lwr <- ys
    probs_lwr <- probs
    ns_lwr <- ns

    ys_upr <- ys
    probs_upr <- probs
    ns_upr <- ns

  }

  stopifnot(length(ys_lwr) == length(ys_upr))
  n <- length(ys_lwr)
  all_ps <- rep(NA_real_, n*(n + 1))
  p_index <- 1
  for(j in 1:n) {
    for(k in 1:j) {
      all_ps[p_index]     <- my_psinib(sum(ys_lwr[k:j]), size = ns_lwr[k:j], prob = probs_lwr[k:j])
      all_ps[p_index + 1] <- my_psinib(sum(ys_upr[k:j]) - 1L, size = ns_upr[k:j], prob = probs_upr[k:j], lower.tail = FALSE)
      p_index <- p_index + 2
    }
  }
  stopifnot(all(!is.na(all_ps)))
  stopifnot(length(all_ps) == n * (n + 1))
  min(p.adjust(all_ps, method = adjust.method))
}

brier_resampling_p <- function(x,y, B = 1000) {
  actual_brier <- sum((x-y)^2)
  brier_null <- replicate(B, {
    yrep <- rbinom(length(x), size = 1, prob = x)
    sum((x - yrep)^2)
  })
  #print(hist(brier_null))
  mean(actual_brier <= brier_null)
}

brier_resampling_history <- function(x,y, B = 1000) {
  actual_brier <- cumsum((x-y)^2)
  brier_null <- replicate(B, {
    yrep <- rbinom(length(x), size = 1, prob = x)
    cumsum((x - yrep)^2)
  })
  #print(hist(brier_null))
  rowMeans(sweep(brier_null, MARGIN = 1, STATS = actual_brier, FUN = `>=`))
}

miscalibration_single <- function(x,y) {
  ord <- order(x, -y)
  x <- x[ord]
  y <- y[ord]
  #CEP_pav <- stats::isoreg(y)$yf
  CEP_pav <- monotone::monotone(y)
  #Using brier score
  Sc <- mean((CEP_pav - y)^2)
  mean((x - y) ^2) - Sc
}

# Faster reimplementation from https://www.pnas.org/doi/full/10.1073/pnas.2016191118#sec-4
# and the reliabilitydiag package
# TODO check licensing
miscalibration_resampling_p <- function(x,y, B = 1000) {
  actual_miscalibration <- miscalibration_single(x,y)
  misc_null <- replicate(B, {
    yrep <- rbinom(length(x), size = 1, prob = x)
    miscalibration_single(x, yrep)
  })
  max(mean(actual_miscalibration <= misc_null), 0.5/B)
}

CEP_bins <- function(x, y) {
  ord <- order(x, -y)
  x <- x[ord]
  y <- y[ord]
  CEP_pav <- monotone::monotone(y)
  bins <- rle(CEP_pav)
  break_locations <- cumsum(head(bins$lengths, -1))
  return(x[break_locations])
}

hosmer_lemeshow_p <- function(x,y, breaks = CEP_bins(x,y), log.p = FALSE) {
  x_cut <- cut(x, breaks = c(0, breaks, 1.001), labels = 1:(length(breaks) + 1), include.lowest = TRUE)
  stopifnot(all(!is.na(x_cut)))
  observed <- xtabs(cbind(y0 = 1 - y, y1 = y) ~ x_cut)
  expected <- xtabs(cbind(yhat0 = 1 - x, yhat1 = x) ~
                      x_cut)
  chisq <- sum((observed - expected)^2/expected)
  G <- length(unique(x_cut))
  return(pchisq(chisq, G - 2, lower.tail = FALSE, log.p = log.p))
}
