test_func_DAP_t <- function(probs, mu) {
  t <- t.test(probs, mu = mu)
  data.frame(p = t$p.value, ci_low = t$conf.int[1], ci_high = t$conf.int[2])
}

test_DAP_power <- function(probs, size, test_func, mu = 0.5, B = 100) {
   res_list <- future.apply::future_replicate(B, simplify = FALSE, {
     p_to_test <- sample(probs, size = size)
     test_func(probs = probs, mu = mu)
   })
   res_df <- do.call(rbind, res_list)
   res_df |> mutate(ci_low = pmax(ci_low, 0), ci_high = pmin(ci_high, 1),
                    ci_width = ci_high - ci_low) |>
     summarise(power = mean(p <= 0.05), mean_ci_width = mean(ci_width), sd_ci_width = sd(ci_width))
}
