#' Scale summary functions.
#'
#' @inheritParams convert_envelope
#' @param scaling The name of the scaling to use. Options include 'none',
#'   'env', 'envdir' and 'scale'. 'env' is default.
#' @param ... Further arguments passed to the chosen scaling function.
#' @export
scale_curves <- function(curve_set, scaling = 'env', ...) {
    curve_set <- convert_envelope(curve_set)

    possible_scalings <- c('env', 'envdir', 'scale', 'none')
    if (length(scaling) != 1L || !(scaling %in% possible_scalings)) {
        stop('scaling must be one of the following: ',
             paste(possible_scalings, collapse = ','))
    }
    scaler <- switch(scaling,
                     env = env_scaling,
                     envdir = envdir_scaling,
                     scale = scale_scaling,
                     none = identity)

    scaled_set <- scaler(curve_set, ...)
    scaled_set
}

#' Turns a divisor into a coeff.
#'
#' Takes the inverse of the input and replaces non-finite values with 1.
divisor_to_coeff <- function(x) {
    y <- 1 / x
    y[!is.finite(y)] <- 1
    y
}

#' Multiply by a coefficient.
#'
#' @inheritParams convert_envelope
#' @param coeff The coefficient vector, often of the length of one curve.
weigh_curves <- function(curve_set, coeff) {
    curve_set[['obs']] <- coeff * curve_set[['obs']]
    curve_set[['sim_m']] <- coeff * curve_set[['sim_m']]
    if (length(curve_set[['theo']]) > 0L) {
        curve_set[['theo']] <- coeff * curve_set[['theo']]
    }
    curve_set
}

#' Check for an increasing two-element vector of probabilities.
check_probs <- function(probs) {
    # Leave further validity checking of probs and type to quantile.
    if (length(probs) != 2L || !all(is.finite(probs)) ||
        probs[1] >= probs[2]) {
        stop('probs must be a two-element vector with both values finite ',
             'and within [0, 1]. The first value must be smaller than the ',
             'second.')
    }
}

#' Env scaling.
#'
#' @inheritParams convert_envelope
#' @param probs A two-element vector containing the lower and upper
#'   quantiles for the envelope, in that order and on the interval [0, 1].
#'   The default values are 0.025 and 0.975 as in the article by Møller and
#'   Berthelsen.
#' @param type The type of algorithm to use for the
#'   \code{\link[stats]{quantile}} function. Type 8 is the default.
#' @param ... Further arguments passed to quantile.
#' @return A scaled curve_set.
#' @references J. Møller and K. K. Berthelsen, “Transforming spatial point
#'   processes into Poisson processes using random superposition,” Advances
#'   in Applied Probability, vol. 44, no. 1, pp. 42–62, 2012.
env_scaling <- function(curve_set, probs = c(0.025, 0.975), type = 8, ...) {
    check_probs(probs)

    # Dimensions: 2, r_idx.
    quant_m <- apply(curve_set[['sim_m']], 1, quantile, probs = probs,
                     type = type, ...)
    # upper - lower
    div <- as.vector(diff(quant_m))
    res <- weigh_curves(curve_set, divisor_to_coeff(div))
    res
}

#' Weigh a matrix or a vector with two different coeffs depending on which
#' side of middle curve each element is.
#'
#' Used by \code{\link{envdir_scaling}}.
weigh_both_sides <- function(x, upper_coeff, lower_coeff) {
    if (is.matrix(x)) {
        dims <- dim(x)
        y <- matrix(0, nrow = dims[1], ncol = dims[2],
                    dimnames = dimnames(x))
    } else if (is.vector(x)) {
        y <- numeric(length(x))
        names(y) <- names(x)
    } else {
        stop('x must be either a matrix or a vector.')
    }

    upper_idx <- x > 0
    lower_or_equal_idx <- !upper_idx

    y[upper_idx] <- upper_coeff * x[upper_idx]
    y[lower_or_equal_idx] <- lower_coeff * x[lower_or_equal_idx]
    y
}

#' EnvDir scaling.
#'
#' @details Notice that this scaling is only defined for residuals.
#'
#' @inheritParams convert_envelope
#' @inheritParams env_scaling
#' @return A scaled curve_set.
envdir_scaling <- function(curve_set, probs = c(0.025, 0.975), type = 8,
                           ...) {
    check_probs(probs)
    check_residualness(curve_set)

    # Dimensions: 2, r_idx.
    quant_m <- apply(curve_set[['sim_m']], 1, quantile, probs = probs,
                     type = type, ...)
    abs_coeff <- divisor_to_coeff(abs(quant_m))
    lower_coeff <- abs_coeff[1, , drop = TRUE]
    upper_coeff <- abs_coeff[2, , drop = TRUE]

    res <- with(curve_set, list(r = r,
                                obs = weigh_both_sides(obs, upper_coeff,
                                                       lower_coeff),
                                sim_m = weigh_both_sides(sim_m, upper_coeff,
                                                         lower_coeff)))
    res[['is_residual']] <- TRUE

    res <- create_curve_set(res)
    res
}

#' Scale with the standard deviation.
#'
#' @inheritParams convert_envelope
#' @return A scaled curve_set.
scale_scaling <- function(curve_set) {
    sim_sd <- apply(curve_set[['sim_m']], 1, sd)
    res <- weigh_curves(curve_set, divisor_to_coeff(sim_sd))
    res
}
