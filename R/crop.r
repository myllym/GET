#' Crop the curves to a certain interval
#'
#' Crop the curves to a certain interval
#'
#'
#' The curves can be cropped to a certain interval defined by the arguments r_min and r_max.
#' The interval should generally be chosen carefully for classical deviation tests.
#' @param curve_set A curve_set (see \code{\link{create_curve_set}}) or
#' an \code{envelope} object of \pkg{spatstat}. If an envelope object is given,
#' it must contain the summary functions from the simulated patterns which can be
#' achieved by setting savefuns = TRUE when calling the \code{envelope} function.
#' @param r_min The minimum radius to include.
#' @param r_max The maximum radius to include.
#' @return A curve_set object containing the cropped summary functions and
#'   the cropped radius vector.
#' @export
crop_curves <- function(curve_set, r_min = NULL, r_max = NULL) {
  if(!is.null(r_min) | !is.null(r_max)) if(!is.vector(curve_set$r)) stop("curve_set$r is not a vector: r_min and r_max cannot be used.")
  curve_set <- convert_envelope(curve_set, allow_Inf_values = TRUE)

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
      return(check_curve_set_content(curve_set, allow_Inf_values = FALSE))
    }
  }

  if(length(cut_idx) < 1L) {
    stop('r_min and r_max cropped everything away.')
  }

  curve_set[['r']] <- r[cut_idx]
  curve_set[['funcs']] <- curve_set[['funcs']][cut_idx, , drop = FALSE]
  theo <- curve_set[['theo']]
  if(!is.null(theo)) curve_set[['theo']] <- theo[cut_idx]
  check_curve_set_content(curve_set, allow_Inf_values = FALSE)
  curve_set
}
