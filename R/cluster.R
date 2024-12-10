setup_cluster <- function() {
  mirai::daemons(min(parallelly::availableCores(), 40), dispatcher = FALSE, force = FALSE)
  future::plan(future.mirai::mirai_cluster)
  options(SBC.generator_chunk_size = 500)
}
