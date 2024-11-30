# Deviation of curves.
#
# @param curve_set A residual curve_set object. Can be obtained by using
#   residual().
# @param measure The deviation measure to use. Default is 'max'. Must be
#   one of the following: 'max', 'int' or 'int2'.
# @param alternative Only used for measure = 'max'!
# @return A deviation_set object. The list has two elements: obs and sim.
#   obs is scalar while sim is a vector with at least one element.
deviation <- function(curve_set, measure = 'max',
                      alternative = c("two.sided", "greater", "less")) {
  possible_measures <- c('max', 'int', 'int2')
  alternative <- match.arg(alternative)

  # deviation() should not accept envelope objects as they are not in
  # residual form.
  if(!inherits(curve_set, 'curve_set')) {
    stop('curve_set must have class "curve_set".')
  }
  check_curve_set_content(curve_set)
  check_residualness(curve_set)

  if(length(measure) != 1L || !(measure %in% possible_measures)) {
    stop('measure must be exactly one of the following: ',
         paste(possible_measures, collapse = ', '))
  }
  all_curves <- data_and_sim_curves(curve_set)

  switch(measure,
         'max' = {
           switch(alternative,
                  two.sided = {
                    res <- apply(abs(all_curves), 1, max)
                  },
                  greater = {
                    res <- apply(all_curves, 1, max)
                  },
                  less = {
                    res <- apply(-all_curves, 1, max)
                  })
         },
         'int' = {
           res <- rowSums(abs(all_curves))
         },
         'int2' = {
           res <- rowSums((all_curves)^2)
         })
  res
}
