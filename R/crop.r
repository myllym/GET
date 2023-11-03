#' Crop the curves
#'
#' Crop the curves to a certain interval, or crop missing and infinite argument
#' values from the curves
#'
#'
#' The curves can be cropped to a certain interval defined by the arguments r_min and r_max.
#' Also the argument values of the sets of curves which have missing or infinite
#' values for any of the curves can be removed from the set (\code{allfinite = TRUE}).
#' The interval should generally be chosen carefully for classical deviation tests.
#' @param curve_set A \code{\link{curve_set}} object, or
#' an \code{envelope} object of \pkg{spatstat}. If an envelope object is given,
#' it must contain the summary functions from the simulated patterns which can be
#' achieved by setting savefuns = TRUE when calling the \code{envelope} function.
#' @param allfinite Logical. TRUE means that the argument values where any of the
#' curves have missing or infinite values are removed. FALSE means that only
#' \code{r_min} and \code{r_max} apply.
#' @param r_min The minimum radius to include.
#' @param r_max The maximum radius to include.
#' @return A curve_set object containing the cropped summary functions and
#'   the cropped radius vector.
#' @export
crop_curves <- function(curve_set, allfinite = TRUE, r_min = NULL, r_max = NULL) {
  if(!is.null(r_min) | !is.null(r_max)) if(!is.vector(curve_set$r))
    stop("curve_set$r is not a vector: r_min and r_max cannot be used.")
  curve_set <- convert_envelope(curve_set, allfinite=FALSE, verbose=FALSE)

  n_r_min <- length(r_min)
  if(n_r_min > 0L && (n_r_min != 1L || !is.finite(r_min))) {
    stop('r_min must be a finite scalar value or NULL.')
  }

  n_r_max <- length(r_max)
  if(n_r_max > 0L && (n_r_max != 1L || !is.finite(r_max))) {
    stop('r_max must be a finite scalar value or NULL.')
  }

  r <- curve_set[['r']]

  if(n_r_min == 1L) {
    if(n_r_max == 1L) {
      if(r_min >= r_max) {
        stop('r_min must be smaller than r_max.')
      }
      cut_idx <- which(r >= r_min & r <= r_max)
    }
    else {
      cut_idx <- which(r >= r_min)
    }
  } else {
    if(n_r_max == 1L) {
      cut_idx <- which(r <= r_max)
    }
    else {
      cut_idx <- seq_along(r)
    }
  }

  if(length(cut_idx) < 1L) {
    stop('r_min and r_max cropped everything away.')
  }

  curve_set <- curve_set[cut_idx, ]

  if(allfinite) {
    argmissing <- curve_set_argmissing(curve_set)
    if(sum(!argmissing) < 1L) {
      stop('allfinite cropped everything away.')
    }
    curve_set <- curve_set[!argmissing, ]
  }

  check_curve_set_content(curve_set, allfinite=allfinite)
  curve_set
}
