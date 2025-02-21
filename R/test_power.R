test_func_DAP_t <- function(probs, mu) {
  if(sd(probs) == 0) {
    if(abs(mean(probs) - mu) < 1e-8) {
      return(data.frame(p = 1, ci_low = mu - 1e-8, ci_high = mu + 1e-8))
    } else {
      return(data.frame(p = NA_real_, ci_low = mean(probs) - 1e-8, ci_high = mean(probs) + 1e-8))
    }
  } else {
      t <- t.test(probs, mu = mu)
      return(data.frame(p = t$p.value, ci_low = t$conf.int[1], ci_high = t$conf.int[2]))
  }
}

test_func_DAP_bayesT <- function(probs, mu) {
  uniform_sd <- sqrt(1/12) # Assume alternative matching the variance of the uniform distribution
    # Handling some edge cases typically with few sims where the Bayes factors
    # don't play nicely
    if(sd(probs) == 0) {
      if(abs(mean(probs) - mu) < 1e-8) {
        return(data.frame(p = 1, ci_low = mu - 1e-8, ci_high = mu + 1e-8))
      } else {
        return(data.frame(p = NA_real_, ci_low = mean(probs) - 1e-8, ci_high = mean(probs) + 1e-8))
      }
    } else if(sd(probs) < 0.0001) {
      if(abs(mean(probs) - expected) < 0.001) {
        return(data.frame(p = 1, ci_low = mu - 0.001, ci_high = mu + 0.001))
      } else {
        return(data.frame(p = 2 * pnorm(-abs(mean(probs) - mu), 0, 0.001), ci_low = mean(probs) - 0.001, ci_high = mu + 0.001))
      }
    } else {
      bf <- tryCatch(BayesFactor::ttestBF(probs, mu = mu, rscale = uniform_sd),
                     error = \(e) { message("ttestBF failed for c(",paste0(probs, collapse = ", "), ")");  NA_real_ })
      nullConnection <- file(nullfile(), open = "w")
      bf_posterior <- tryCatch(
        { sink(file = nullConnection, type = "message")
          BayesFactor::ttestBF(probs, mu = mu, rscale = uniform_sd, posterior = TRUE, iterations = 4000)@data$y} ,
        finally = {sink(type = "message"); close(nullConnection)},
                     error = \(e) { message("ttestBF failed for c(",paste0(probs, collapse = ", "), ")");  NA_real_ })
      return(data.frame(p = plogis(-bf@bayesFactor[1, "bf"]), ci_low = quantile(bf_posterior, 0.025), ci_high = quantile(bf_posterior, 0.975)))
    }
}

test_func_DAP_gaffke <- function(probs, mu, B = 4000) {
   m <- gaffke_m(probs, B = B)
   ci <- gaffke_ci_from_m(m)
   data.frame(p = gaffke_p_from_m(m, mu, B = B), ci_low = ci[1], ci_high = ci[2])
}

test_DAP_power <- function(probs, size, test_func, mu = 0.5, B = 2000) {
   res_list <- future.apply::future_replicate(
     B,
     simplify = FALSE,
     future.globals = c("probs", "size", "mu", "test_func"),
     expr = {
       p_to_test <- sample(probs, size = size)
       test_func(probs = p_to_test, mu = mu)
   })
   res_df <- do.call(rbind, res_list)
   res_df
}
