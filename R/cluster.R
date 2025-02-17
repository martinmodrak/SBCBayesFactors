setup_cluster <- function(max_local_workers = Inf, max_server_workers = 140, force = FALSE) {
  is_local <- parallelly::availableCores() < 80
  if(is_local) {
    n_workers <- min(parallelly::availableCores(), max_local_workers)
  } else {
    n_workers <- min(parallelly::availableCores(), max_server_workers)
  }
  mirai::daemons(n_workers, dispatcher = FALSE, force = force)
  future::plan(future.mirai::mirai_cluster)
  options(SBC.generator_chunk_size = 500)
}
