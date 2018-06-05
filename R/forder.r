# Continuous pointwise ranks
#
# Calculate continuous pointwise ranks of the curves from the largest (smallest rank)
# to the smallest (largest rank).
cont_pointwise_hiranks <- function(data_and_sim_curves) {
  Nfunc <- dim(data_and_sim_curves)[1]
  nr <- dim(data_and_sim_curves)[2]
  RRR <- array(0, c(Nfunc, nr))
  for(i in 1:nr) {
    y = data_and_sim_curves[,i]
    ordery = order(y, decreasing = TRUE) # index of function values at "i" from largest to smallest
    for(j in 2:(Nfunc-1)) {
      if(y[ordery[j-1]] == y[ordery[j+1]])
        RRR[ordery[j],i] = j-1
      else
        RRR[ordery[j],i] = j-1+(y[ordery[j-1]]-y[ordery[j]])/(y[ordery[j-1]]-y[ordery[j+1]])
    }
    RRR[ordery[1],i] = exp(-(y[ordery[1]]-y[ordery[2]])/(y[ordery[2]]-y[ordery[Nfunc]]))
    RRR[ordery[Nfunc],i] = Nfunc
    while(j <= Nfunc-2) {
      k = 1
      if(y[ordery[j]] == y[ordery[j+2]]) {
        k = 3; S = 3*j+3
        while(j+k <= Nfunc && y[ordery[j]] == y[ordery[j+k]]) { k = k+1; S = S+j+k }
        for(t in j:(j+k-1)) { RRR[ordery[t],i] = S/k }
      }
      j = j+k
    }
  }
  RRR
}

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
#' the functions are ordered from the most extreme one to the least extreme one
#' by one of the following measures (specified by the argument \code{measure}):
#'
#' \itemize{
#'  \item \code{'rank'}: extreme rank (Myllymäki et al., 2017).
#' The extreme rank \eqn{R_i}{R_i} is defined as the minimum of pointwise ranks of the curve
#' \eqn{T_i(r)}{T_i(r)}, where the pointwise rank is the rank of the value of the curve for a
#' specific r-value among the corresponding values of the s other curves such that the lowest
#' ranks correspond to the most extreme values of the curves. How the pointwise ranks are determined
#' exactly depends on the whether a one-sided (\code{alternative} is "less" or "greater") or the
#' two-sided test (\code{alternative="two.sided"}) is chosen, for details see
#' Mrkvička et al. (2017, page 1241) or Mrkvička et al. (2018, page 6).
#'  \item \code{'erl'}: extreme rank length extreme rank lengths were proposed by Myllymäki et al. (2017)
#'  as a refinement to extreme ranks, because the extreme ranks can contain many ties.
#'  Considering the vector of pointwise ordered ranks \eqn{\mathbf{R}_i}{RP_i} of the ith curve,
#'  the extreme rank length measure \eqn{R_i^{\text{erl}}}{Rerl_i} is equal to
#' \deqn{R_i^{\text{erl}} = \frac{1}{s+1}\sum_{j=1}^{s+1} \1(\mathbf{R}_j \prec \mathbf{R}_i)}{Rerl_i = \sum_{j=1}^{s} 1(RP_j "<" RP_i) / (s + 1)}
#' where \eqn{\mathbf{R}_j \prec \mathbf{R}_i}{RP_j "<" RP_i} if and only if
#' there exists \eqn{n\leq d}{n<=d} such that for the first k, \eqn{k<n}{k<n}, pointwise ordered
#' ranks of \eqn{\mathbf{R}_j}{RP_j} and \eqn{\mathbf{R}_i}{RP_i} are equal and the n'th rank of
#' \eqn{\mathbf{R}_j}{RP_j} is smaller than that of \eqn{\mathbf{R}_i}{RP_i}.
#'  \item \code{'cont'}:
#'  \item \code{'area'}:
#'  \item \code{'max'} and \code{'int'} and \code{'int2'}:
#' Further options for the \code{measure} argument that can be used together with \code{scaling}.
#' See the help in \code{\link{deviation_test}} for these options of \code{measure} and \code{scaling}.
#' These measures are largest for the most extreme functions and smallest for the least extreme ones.
#' The arguments \code{use_theo} and \code{probs} are relevant for these measures only (otherwise ignored).
#' }
#'
#' @return A vector containing one of the above mentioned measures k for each of the functions
#' in the curve set. If the component \code{obs} in the curve set is a vector, then its measure
#' will be the first component (named 'obs') in the returned vector.
#'
#' @inheritParams crop_curves
#' @param measure The measure to use to order the functions from the most extreme to the least extreme
#' one. Must be one of the following: 'rank', 'erl', 'cont', 'area', 'max', 'int', 'int2'. Default is 'erl'.
#' @param scaling The name of the scaling to use if measure is 'max', 'int' or 'int2'.
#' Options include 'none', 'q', 'qdir' and 'st', where 'qdir' is the default.
#' @param alternative A character string specifying the alternative hypothesis.
#' Must be one of the following: "two.sided" (default), "less" or "greater".
#' The last two options only available for \code{type = 'rank'} and \code{type = 'erl'}.
#' @param use_theo Logical. When calculating the measures 'max', 'int', 'int2',
#'  should the theoretical function from \code{curve_set} be used (if 'theo' provided),
#'  see \code{\link{deviation_test}}.
#' @param probs A two-element vector containing the lower and upper
#'   quantiles for the measure 'q' or 'qdir', in that order and on the interval [0, 1].
#'   The default values are 0.025 and 0.975, suggested by Myllymäki et al. (2015, 2017).
#' @export
#' @examples
#' if(requireNamespace("fda", quietly = TRUE)) {
#'   curve_set <- create_curve_set(list(r = as.numeric(row.names(fda::growth$hgtf)),
#'                                      obs = fda::growth$hgtf))
#'   plot(curve_set, ylab="height")
#'   forder(curve_set, measure = "max", scaling="qdir")
#'   forder(curve_set, measure = "rank")
#'   forder(curve_set, measure = "erl")
#' }
forder <- function(curve_set, r_min = NULL, r_max = NULL,
                    measure = 'erl', scaling = 'qdir',
                    alternative=c("two.sided", "less", "greater"),
                    use_theo = TRUE, probs = c(0.025, 0.975)) {
  possible_measures <- c('rank', 'erl', 'cont', 'area', 'max', 'int', 'int2')
  if(!(measure %in% possible_measures)) stop("Unreasonable measure argument!\n")

  curve_set <- crop_curves(curve_set, r_min = r_min, r_max = r_max)

  if(measure %in% c('max', 'int', 'int2')) {
    curve_set <- residual(curve_set, use_theo = use_theo)
    if(scaling %in% c('q', 'qdir')) curve_set <- scale_curves(curve_set, scaling = scaling, probs = probs)
    else curve_set <- scale_curves(curve_set, scaling = scaling)
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

    # calculate pointwise ranks
    if(measure %in% c('rank', 'erl')) {
      loranks <- apply(data_and_sim_curves, MARGIN=2, FUN=rank, ties.method = "average")
      hiranks <- Nfunc+1-loranks
    }
    else { # 'cont', 'area'
      hiranks <- cont_pointwise_hiranks(data_and_sim_curves)
      loranks <- cont_pointwise_hiranks(-data_and_sim_curves)
    }
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

    # calculate measures from the pointwise ranks
    switch(measure,
           rank = {
             distance <- apply(allranks, MARGIN=1, FUN=min) # extreme ranks R_i
           },
           erl = {
             sortranks <- apply(allranks, 1, sort) # curves now represented as columns
             lexo_values <- do.call("order", split(sortranks, row(sortranks))) # indices! of the functions from the most extreme to least extreme one
             newranks <- 1:Nfunc
             distance <- newranks[order(lexo_values)] / Nfunc # ordering of the functions by the extreme rank lengths
           },
           cont = {
             distance <- apply(allranks, MARGIN=1, FUN=min)
           },
           area = {
             RRRm <- apply(allranks, MARGIN=1, FUN=min)
             RRRm <- ceiling(RRRm) # = R_i
             distance <- array(0, Nfunc)
             for(j in 1:Nfunc) distance[j] <- sum(allranks[j, allranks[j,] <= RRRm[j]] - RRRm[j]+1) + RRRm[j]*nr
           })
  }
  names(distance) <- rownames(data_and_sim_curves)
  distance
}
