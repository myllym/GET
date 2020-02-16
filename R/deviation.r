# Deviation of curves.
#
# @param curve_set A residual curve_set object. Can be obtained by using
#   residual().
# @param measure The deviation measure to use. Default is 'max'. Must be
#   one of the following: 'max', 'int' or 'int2'.
# @param ... Arguments to be passed to the measure function, if applicable.
# @return A deviation_set object. The list has two elements: obs and sim.
#   obs is scalar while sim is a vector with at least one element.
#' @importFrom methods is
deviation <- function(curve_set, measure = 'max', ...) {
  possible_measures <- c('max', 'int', 'int2')

  # deviation() should not accept envelope objects as they are not in
  # residual form.
  if(!methods::is(curve_set, 'curve_set')) {
    stop('curve_set must have class "curve_set".')
  }
  check_curve_set_content(curve_set)
  check_residualness(curve_set)

  if(length(measure) != 1L || !(measure %in% possible_measures)) {
    stop('measure must be exactly one of the following: ',
         paste(possible_measures, collapse = ', '))
  }
  data_and_sim_curves <- data_and_sim_curves(curve_set)

  switch(measure,
         'max' = {
           res <- apply(abs(data_and_sim_curves), 1, max)
         },
         'int' = {
           res <- rowSums(abs(data_and_sim_curves))
         },
         'int2' = {
           res <- rowSums((data_and_sim_curves)^2)
         })
  res
}

# Check the content validity of a potential deviation_set object.
#
# @param deviation_set A potential deviation_set object.
check_deviation_set_content <- function(deviation_set) {
  possible_names <- c('obs', 'sim')

  n <- length(deviation_set)
  if(n < 1L) {
    stop('deviation_set must have some elements.')
  }
  if(!is.list(deviation_set)) {
    stop('deviation_set must be a list.')
  }

  name_vec <- names(deviation_set)
  if(length(name_vec) != n) {
    stop('Every element of deviation_set must be named.')
  }
  if(!all(name_vec %in% possible_names)) {
    stop('deviation_set should contain exactly these elements: ',
         paste(possible_names, collapse = ', '))
  }

  obs <- deviation_set[['obs']]
  if(length(obs) != 1L) {
    stop('deviation_set[["obs"]] must be a scalar.')
  }
  if(!is.vector(obs)) {
    stop('deviation_set[["obs"]] must be a vector.')
  }
  if(!all(is.numeric(obs)) || !all(is.finite(obs))) {
    stop('deviation_set[["obs"]] must have a finite numeric value.')
  }

  sim <- deviation_set[['sim']]
  if(length(sim) < 1L) {
    stop('deviation_set[["sim"]] must have at least one value.')
  }
  if(!is.vector(sim)) {
    stop('deviation_set[["sim"]] must be a vector.')
  }
  if(!all(is.numeric(sim)) || !all(is.finite(sim))) {
    stop('deviation_set[["sim"]] must have only finite numeric values.')
  }

  deviation_set
}

# Check the object.
#
# @param deviation_set A potential deviation_set object.
check_deviation_set <- function(deviation_set) {
  if(!inherits(deviation_set, 'deviation_set')) {
    stop('deviation_set must have class "deviation_set".')
  }
  check_deviation_set_content(deviation_set)
  deviation_set
}

# Create a deviation set out of a list in the right form.
#
# @param deviation_set The list to be turned into a deviation_set object.
create_deviation_set <- function(deviation_set) {
  check_deviation_set_content(deviation_set)
  class(deviation_set) <- 'deviation_set'
  deviation_set
}

# Check class.
#
# @param x An object to be checked.
is.deviation_set <- function(x) inherits(x, 'deviation_set')
