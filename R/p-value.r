#' Estimate p-value.
#'
#' @param x The first argument.
#' @param ... Additional arguments.
#' @export
#' @seealso \code{\link{estimate_p_value.default}}
estimate_p_value <- function (x, ...) UseMethod('estimate_p_value')

# FIXME: Do we need to consider NA values at some point in life? Our
# methods should not produce them but if we export this function, they
# should be handled.
# FIXME: How about both two-sided and one-sided?
# FIXME: north_note_2002 uses the term anti-conservative for the incorrect p-value.
#' Estimate p-value.
#'
#' Estimates the p-value of the given observation for the given set of Monte
#' Carlo samples. User can choose which method is used to treat possible
#' tied values.
#'
#' @usage \method{estimate_p_value}{default}(obs, sim_vec, ties = 'midrank', ...)
#'
#' @param obs The data sample. A scalar real value. Must not be
#'   NULL.
#' @param sim_vec The Monte Carlo samples. A vector of real values.
#'   Must not be NULL.
#' @param ties The method to treat tied values with. If one or more of the
#'   elements of sim_vec are equal to obs, how should the rank of obs be
#'   determined? For 'conservative' the resulting p-value will be the
#'   highest possible. For 'liberal' the p-value will be the lowest
#'   possible. For 'random' the rank of the obs within the tied values is
#'   uniformly sampled so that the resulting p-value is at most the
#'   conservative option and at least the liberal option. For 'midrank'
#'   the mid-rank within the tied values is taken. 'midrank' is the default.
#' @return The p-value estimate. A scalar real value between 0 and 1.
#'
#' @references Hájek & Šidák & Sen. Theory of Rank Tests. 1999. ff. 130.
#' @export
estimate_p_value.default <- function(obs, sim_vec, ties = 'midrank') {
    if (length(obs) != 1L || !is.finite(obs) || !is.numeric(obs)) {
        stop('obs must be a scalar finite real value.')
    }
    n_sim <- length(sim_vec)
    if (n_sim < 1L || !all(is.finite(sim_vec)) ||
        !all(is.numeric(sim_vec))) {
        stop('sim_vec must have at least one element and the ',
             'elements must be finite and real.')
    }
    possible_ties <- c('midrank', 'random', 'conservative', 'liberal')
    if (length(ties) != 1L || !(ties %in% possible_ties)) {
        stop('ties must be exactly one of the following: ',
             paste(possible_ties, collapse=', '))
    }

    # FIXME: Ugly redundancy.
    if (ties %in% 'conservative') {
        n_less <- sum(sim_vec < obs)
    } else if (ties %in% 'liberal') {
        n_less <- sum(sim_vec <= obs)
    } else if (ties %in% 'midrank') {
        n_smaller <- sum(sim_vec < obs)
        n_equal <- sum(sim_vec == obs)
        n_less <- n_smaller + n_equal / 2L
    } else if (ties %in% 'random') {
        n_smaller <- sum(sim_vec < obs)
        n_equal <- sum(sim_vec == obs)
        n_rand_smaller <- sample(seq(0L, n_equal, by = 1L), size = 1L)
        n_less <- n_smaller + n_rand_smaller
    } else {
        stop('ties argument not recognized: ', ties)
    }
    n_all <- n_sim + 1L
    p_estimate <- 1 - n_less / n_all
    p_estimate
}

#' Estimate p-value.
#'
#' Check \code{\link{estimate_p_value.default}} for details.
#'
#' @usage \method{estimate_p_value}{deviation_set}(deviation_set, ...)
#'
#' @param deviation_set A deviation_set object. Use deviation() to get one.
#' @param ... Additional parameters to be passed to \code{\link{estimate_p_value.default}}.
#' @return The p-value.
estimate_p_value.deviation_set <- function(deviation_set, ...) {
    check_deviation_set(deviation_set)
    p <- with(deviation_set,
              estimate_p_value.default(obs, sim, ...))
    p
}
