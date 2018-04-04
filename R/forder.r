#' Functional ordering
#'
#' Calculates different measures for ordering the functions (or vectors)
#' from the most extreme to least extreme one
#'
#'
#' Given a \code{curve_set} (see \code{\link{create_curve_set}} for how to create such an object)
#' or an \code{\link[spatstat]{envelope}} object,
#' which contains both the data curve (or function or vector) \eqn{T_1(r)}{T_1(r)} and
#' the simulated curves \eqn{T_2(r),\dots,T_{s+1}(r)}{T_2(r),...,T_(s+1)(r)},
#' the functions are ordered by one of the following measures.
#'
#' First the curves in the curve set are ranked from the most extreme one to the least extreme one
#' either by using the extreme ranks R_i and/or the extreme rank lengths \eqn{R_i^{\text{erl}}}{Rerl_i}.
#' The option \code{measure = "rank"}, the ordering is based on the extreme ranks R_i.
#' The extreme rank is defined as the minimum of pointwise ranks of the curve \eqn{T_i(r)}{T_i(r)},
#' where the pointwise rank is the rank of the value of the curve for a specific r-value among the
#' corresponding values of the s other curves such that the lowest ranks correspond to the most extreme
#' values of the curves. How the pointwise ranks are determined exactly depends on the whether a
#' one-sided (\code{alternative} is "less" or "greater") or the two-sided test (\code{alternative="two.sided"}) is
#' chosen, for details see Mrkvička et al. (2017, page 1241) or Mrkvička et al. (2018, page 6).
#'
#' The extreme ranks can contain many ties, for which reason Myllymäki et al. (2017) proposed the
#' extreme rank length ordering. This ordering can be obtained by specifying \code{measure = "erl"}.
#' Considering the vector of pointwise ordered ranks
#' \eqn{\mathbf{R}_i}{RP_i} of the ith curve, the extreme rank length measure
#' \eqn{R_i^{\text{erl}}}{Rerl_i}
#' is equal to
#' \deqn{R_i^{\text{erl}} = \frac{1}{s+1}\sum_{j=1}^{s+1} \1(\mathbf{R}_j \prec \mathbf{R}_i)}{Rerl_i = \sum_{j=1}^{s} 1(RP_j "<" RP_i) / (s + 1)}
#' where \eqn{\mathbf{R}_j \prec \mathbf{R}_i}{RP_j "<" RP_i} if and only if
#' there exists \eqn{n\leq d}{n<=d} such that for the first k, \eqn{k<n}{k<n}, pointwise ordered
#' ranks of \eqn{\mathbf{R}_j}{RP_j} and \eqn{\mathbf{R}_i}{RP_i} are equal and the n'th rank of
#' \eqn{\mathbf{R}_j}{RP_j} is smaller than that of \eqn{\mathbf{R}_i}{RP_i}.
#'
#' Further options for the \code{measure} argument are \code{max}, \code{int} and \code{int2}
#' which can be used together with \code{scaling}. See the help in \code{\link{deviation_test}}
#' for these options of measures (measure) and scalings.
#'
#' @return A vector containing one of the above mentioned measures k for each of the functions
#' in the curve set. If the component \code{obs} in the curve set is a vector, then its measure
#' will be the first component (named 'obs') in the returned vector.
#'
#' @inheritParams deviation_test
#' @param curve_set A \code{curve_set} object.
#' @export
forder <- function(curve_set, r_min = NULL, r_max = NULL,
                    measure = 'erl', scaling = 'qdir',
                    alternative=c("two.sided", "less", "greater"),
                    use_theo = TRUE) {
  possible_measures <- c('rank', 'erl', 'max', 'int', 'int2')
  if(!(measure %in% possible_measures)) stop("Unreasonable measure argument!\n")

  curve_set <- crop_curves(curve_set, r_min = r_min, r_max = r_max)

  if(measure %in% c('max', 'int', 'int2')) {
    curve_set <- residual(curve_set, use_theo = use_theo)
    curve_set <- scale_curves(curve_set, scaling = scaling)
    data_and_sim_curves <- data_and_sim_curves(curve_set)
    curve_set_tmp <- create_curve_set(list(r=curve_set[['r']],
                                           obs=data_and_sim_curves[1,],
                                           sim_m=t(data_and_sim_curves[-1,]),
                                           is_residual=curve_set[['is_residual']]))
    devs <- deviation(curve_set_tmp, measure = measure)
    distance <- c(devs$obs, devs$sim)
  }
  else {
    alternative <- match.arg(alternative)
    data_and_sim_curves <- data_and_sim_curves(curve_set)
    Nfunc <- dim(data_and_sim_curves)[1]
    nr <- length(curve_set[['r']])

    loranks <- apply(data_and_sim_curves, MARGIN=2, FUN=rank, ties.method = "average")
    hiranks <- Nfunc+1-loranks

    switch(alternative,
           "two.sided" = {
             allranks <- pmin(loranks, hiranks)
           },
           "less" = {
             allranks <- loranks
           },
           "greater" = {
             allranks <- hiranks
           })

    switch(measure,
           rank = {
             distance <- apply(allranks, MARGIN=1, FUN=min) # extreme ranks R_i
           },
           erl = {
             sortranks <- apply(allranks, 1, sort) # curves now represented as columns
             lexo_values <- do.call("order", split(sortranks, row(sortranks))) # indices! of the functions from the most extreme to least extreme one
             newranks <- 1:Nfunc
             distance <- newranks[order(lexo_values)] / Nfunc # ordering of the functions by the extreme rank lengths
           })
  }
  names(distance) <- rownames(data_and_sim_curves)
  distance
}
