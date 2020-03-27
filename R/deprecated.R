#' @export
central_region2d <- function(image_sets, ...) {
  warning("central_region2d is deprecated. Please use central_region.")
  central_region(curve_sets=image_sets, ...)
}
#' @export
global_envelope_test2d <- function(image_sets, ...) {
  warning("global_envelope_test2d is deprecated. Please use global_envelope_test.")
  global_envelope_test(curve_sets=image_sets, ...)
}
#' @export
graph.fanova2d <- function(nsim, image_set, groups, ...) {
  warning("graph.fanova2d is deprecated. Please use graph.fanova.")
  graph.fanova(nsim=nsim, curve_set=image_set, groups=groups, ..., n.aver=1)
}
#' @export
frank.fanova2d <- function(nsim, image_set, groups, ...) {
  warning("frank.fanova2d is deprecated. Please use frank.fanova.")
  frank.fanova(nsim=nsim, curve_set=image_set, groups=groups, ...)
}
#' @export
graph.flm2d <- function(nsim, formula.full, formula.reduced, image_sets, factors = NULL, ...) {
  warning("graph.flm2d is deprecated. Please use graph.flm.")
  graph.flm(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
            curve_sets=image_sets, factors=factors, ...)
}
#' @export
frank.flm2d <- function(nsim, formula.full, formula.reduced, image_sets, factors = NULL, ...) {
  warning("frank.flm2d is deprecated. Please use frank.flm.")
  frank.flm(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
            curve_sets=image_sets, factors=factors, ...)
}
