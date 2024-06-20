good_check <- function(true_model, p_M1) {
  #TODO log version
  p_wrong <- if_else(true_model == 0, p_M1, 1 - p_M1)

  bf_wrong <- p_wrong / (1 - p_wrong)
  bf_correct <- 1/bf_wrong

  list(
    mean_bf_wrong = mean(bf_wrong),
    e_bf_correct = mean(bf_correct),
    e2_bf_wrong = mean(bf_wrong)^2 + var(bf_wrong),
    moment_diff = mean(bf_correct) - mean(bf_wrong)^2 - var(bf_wrong)
  )
}
