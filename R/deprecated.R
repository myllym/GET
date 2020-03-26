#' @export
global_envelope_test2d <- function(..., image_sets=image_sets) {
  warning("global_envelope_test2d is deprecated. Please use global_envelope_test.")
  global_envelope_test(..., curve_sets=image_sets)
}
#' @export
graph.fanova2d <- function(..., image_set=image_set) {
  warning("graph.fanova2d is deprecated. Please use graph.fanova.")
  graph.fanova(..., curve_set=image_set)
}
#' @export
graph.flm2d <- function(..., image_sets=image_sets) {
  warning("graph.flm2d is deprecated. Please use graph.flm.")
  graph.flm(..., curve_sets=image_sets)
}
#' @export
frank.fanova2d <- function(..., image_set=image_set) {
  warning("frank.fanova2d is deprecated. Please use frank.fanova.")
  frank.fanova(..., curve_set=image_set)
}
#' @export
frank.flm2d <- function(..., image_sets=image_sets) {
  warning("frank.flm2d is deprecated. Please use frank.flm.")
  frank.flm(..., curve_sets=image_sets)
}
