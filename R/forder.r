# Continuous pointwise ranks
#
# Calculate continuous pointwise ranks of the curves from the largest (smallest rank)
# to the smallest (largest rank) for the data in the vector y (corresponding to functions at an argument value (r)).
contrank <- function(y, ordery=order(y, decreasing=TRUE)) {
  Nfunc <- length(y)
  y <- y[ordery]
  ties <- y[1:(Nfunc-2)] == y[3:Nfunc] # same as j = 2:(Nfunc-1); ties <- y[j-1] == y[j+1]
  RR <- 0:(Nfunc-1) # Initialize the vector with j-1 (values for the case of ties)
  RR[1] <- exp(-(y[1]-y[2])/(y[2]-y[Nfunc]))
  if(!any(ties)) {
    RR[2:(Nfunc-1)] <- 1:(Nfunc-2)+(y[1:(Nfunc-2)]-y[2:(Nfunc-1)])/(y[1:(Nfunc-2)]-y[3:Nfunc])
    # same as j <- 2:(Nfunc-1); RR[j] <- j-1+(y[j-1]-y[j])/(y[j-1]-y[j+1])
  }
  else { # The case of some ties
    j <- (2:(Nfunc-1))[!ties]
    jm1 <- j-1
    RR[j] <- jm1+(y[jm1]-y[j])/(y[jm1]-y[j+1])
    # Treat ties
    j <- 1
    while(j <= Nfunc-2) {
      k <- 1
      if(ties[j]) {
        k <- 3; S <- 3*j+3
        while(j+k <= Nfunc && y[j] == y[j+k]) { k <- k+1; S <- S+j+k }
        for(t in j:(j+k-1)) { RR[t] <- S/k }
      }
      j <- j+k
    }
  }
  RR[ordery] <- RR
  RR
}

# @return A function to calculate pointwise rank in different situations given by the measure and the alternative
find_calc_pointwiserank <- function(measure, alternative) {
  if(measure %in% c('rank', 'erl')) {
    avrank <- function(x) { rank(x, ties.method = "average") }
    switch(alternative,
           "two.sided" = {
             function(x) {
               loranks <- avrank(x)
               hiranks <- length(x)+1-loranks
               pmin(loranks, hiranks)
             }
           },
           "less" = {
             function(x) { avrank(x) }
           },
           "greater" = {
             function(x) { length(x)+1-avrank(x) }
           })
  }
  else if(measure %in% c('cont', 'area')) {
    switch(alternative,
           "two.sided" = {
             function(y) {
               ordery <- order(y, decreasing = TRUE)
               hiranks <- contrank(y, ordery)
               loranks <- contrank(-y, rev(ordery))
               pmin(loranks, hiranks)
             }
           },
           "less" = {
             function(y) { contrank(-y) }
           },
           "greater" = {
             function(y) { contrank(y) }
           })
  } else { stop("Internal error in GET.")}
}

# Functionality for functional ordering based on a curve set
individual_forder <- function(curve_set,
                              measure = 'erl', scaling = 'qdir',
                              alternative=c("two.sided", "less", "greater"),
                              use_theo = TRUE,
                              probs = c(0.025, 0.975), quantile.type = 7) {
  possible_measures <- c('rank', 'erl', 'cont', 'area', 'max', 'int', 'int2')
  if(!(measure %in% possible_measures)) stop("Unreasonable measure argument!\n")

  curve_set <- convert_envelope(curve_set)

  if(measure %in% c('max', 'int', 'int2')) {
    curve_set <- residual(curve_set, use_theo = use_theo)
    curve_set <- scale_curves(curve_set, scaling = scaling, probs = probs, type = quantile.type)
    distance <- deviation(curve_set, measure = measure)
  }
  else {
    alternative <- match.arg(alternative)
    data_and_sim_curves <- data_and_sim_curves(curve_set)
    Nfunc <- dim(data_and_sim_curves)[1]
    nr <- curve_set_narg(curve_set)

    # Calculate pointwise ranks for each argument value (r)
    calc_pointwiserank <- find_calc_pointwiserank(measure, alternative)
    for(i in 1:nr) {
      data_and_sim_curves[,i] <- calc_pointwiserank(data_and_sim_curves[,i]) # overwriting curves by their ranks
    }
    allranks <- data_and_sim_curves

    # Calculate measures from the pointwise ranks
    applyfuncs <- function(f) {
      sapply(1:nrow(allranks), function(i) { f(allranks[i,]) })
    }
    switch(measure,
           rank = {
             distance <- applyfuncs(min) # extreme ranks R_i
           },
           erl = {
             sortranks <- applyfuncs(sort)
             lexo_values <- do.call("order", split(sortranks, row(sortranks))) # indices! of the functions from the most extreme to least extreme one
             newranks <- 1:Nfunc
             distance <- newranks[order(lexo_values)] / Nfunc # ordering of the functions by the extreme rank lengths
           },
           cont = {
             distance <- applyfuncs(min)
             if(alternative == "two.sided")
               distance <- distance / ceiling(Nfunc/2)
             else
               distance <- distance / (Nfunc-1)
           },
           area = {
             distance <- applyfuncs(function(x) { sum(pmin(ceiling(min(x)), x)) })
             if(alternative == "two.sided")
               distance <- distance / (ceiling(Nfunc/2)*nr)
             else
               distance <- distance / ((Nfunc-1)*nr)
           })
  }
  names(distance) <- rownames(data_and_sim_curves)
  list(distance = distance, measure = measure)
}

# Functionality for functional ordering based on several curve sets
combined_forder <- function(curve_sets, ...) {
  ntests <- length(curve_sets)
  curve_sets <- check_curve_set_dimensions(curve_sets)

  # 1) First stage: Calculate the functional orderings individually for each curve_set
  res_ls <- lapply(curve_sets, FUN = function(x) { individual_forder(x, ...) })

  # 2) Second stage: ERL ordering
  # Create a curve_set for the ERL test
  k_ls <- lapply(res_ls, FUN = function(x) x$distance)
  k_mat <- do.call(rbind, k_ls, quote=FALSE)
  curve_set_u <- create_curve_set(list(r=1:ntests, obs=k_mat))
  # Construct the one-sided ERL central region
  if(res_ls[[1]]$measure %in% c('max', 'int', 'int2')) alt2 <- "greater"
  else alt2 <- "less"
  individual_forder(curve_set_u, measure="erl", alternative=alt2)
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
#' by one of the following measures (specified by the argument \code{measure}).
#' Note that \code{'erl'}, \code{'cont'} and \code{'area'} were proposed as a refinement to
#' the extreme ranks \code{'rank'}, because the extreme ranks can contain many ties.
#' All of these completely non-parametric measures are smallest for the most extreme functions
#' and largest for the least extreme ones,
#' whereas the deviation measures (\code{'max'}, \code{'int'} and \code{'int2'}) obtain largest values
#' for the most extreme functions.
#' \itemize{
#'  \item \code{'rank'}: extreme rank (Myllymäki et al., 2017).
#' The extreme rank \eqn{R_i}{R_i} is defined as the minimum of pointwise ranks of the curve
#' \eqn{T_i(r)}{T_i(r)}, where the pointwise rank is the rank of the value of the curve for a
#' specific r-value among the corresponding values of the s other curves such that the lowest
#' ranks correspond to the most extreme values of the curves. How the pointwise ranks are determined
#' exactly depends on the whether a one-sided (\code{alternative} is "less" or "greater") or the
#' two-sided test (\code{alternative="two.sided"}) is chosen, for details see
#' Mrkvička et al. (2017, page 1241) or Mrkvička et al. (2018, page 6).
#'  \item \code{'erl'}: extreme rank length (Myllymäki et al., 2017).
#'  Considering the vector of pointwise ordered ranks \eqn{\mathbf{R}_i}{RP_i} of the ith curve,
#'  the extreme rank length measure \eqn{R_i^{erl}}{Rerl_i} is equal to
#' \deqn{R_i^{erl} = \frac{1}{s+1}\sum_{j=1}^{s+1} \mathbf{1}(\mathbf{R}_j "<" \mathbf{R}_i)}{Rerl_i = \sum_{j=1}^{s} 1(RP_j "<" RP_i) / (s + 1)}
#' where \eqn{\mathbf{R}_j "<" \mathbf{R}_i}{RP_j "<" RP_i} if and only if
#' there exists \eqn{n\leq d}{n<=d} such that for the first k, \eqn{k<n}{k<n}, pointwise ordered
#' ranks of \eqn{\mathbf{R}_j}{RP_j} and \eqn{\mathbf{R}_i}{RP_i} are equal and the n'th rank of
#' \eqn{\mathbf{R}_j}{RP_j} is smaller than that of \eqn{\mathbf{R}_i}{RP_i}.
#'  \item \code{'cont'}: continuous rank (Hahn, 2015; Mrkvička et al., 2019)
#' based on minimum of continuous pointwise ranks
#'  \item \code{'area'}: area rank (Mrkvička et al., 2019) based on area between continuous
#'  pointwise ranks and minimum pointwise ranks for those argument (r) values for which pointwise
#'  ranks achieve the minimum (it is a combination of erl and cont)
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
#' @param curve_sets A \code{curve_set} object or a list of \code{curve_set} objects.
#' @param measure The measure to use to order the functions from the most extreme to the least extreme
#' one. Must be one of the following: 'rank', 'erl', 'cont', 'area', 'max', 'int', 'int2'. Default is 'erl'.
#' @param scaling The name of the scaling to use if measure is 'max', 'int' or 'int2'.
#' Options include 'none', 'q', 'qdir' and 'st', where 'qdir' is the default.
#' @param alternative A character string specifying the alternative hypothesis.
#' Must be one of the following: "two.sided" (default), "less" or "greater".
#' The last two options only available for types \code{'rank'}, \code{'erl'},
#' \code{'cont'} and \code{'area'}.
#' @param use_theo Logical. When calculating the measures 'max', 'int', 'int2',
#'  should the theoretical function from \code{curve_set} be used (if 'theo' provided),
#'  see \code{\link{deviation_test}}.
#' @param probs A two-element vector containing the lower and upper
#'   quantiles for the measure 'q' or 'qdir', in that order and on the interval [0, 1].
#'   The default values are 0.025 and 0.975, suggested by Myllymäki et al. (2015, 2017).
#' @param quantile.type As type argument of \code{\link[stats]{quantile}}, how to
#' calculate quantiles for 'q' or 'qdir'.
#' @export
#' @references
#' Hahn U (2015). “A note on simultaneous Monte Carlo tests.” Technical report, Centre for Stochastic Geometry and advanced Bioimaging, Aarhus University.
#'
#' Mrkvička, T., Hahn, U. and Myllymäki, M. (2018). A one-way ANOVA test for functional data with graphical interpretation. arXiv:1612.03608 [stat.ME]
#'
#' Mrkvička, T., Myllymäki, M. and Narisetty, N. N. (2019) New methods for multiple testing in permutation inference for the general linear model. arXiv:1906.09004 [stat.ME]
#'
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2015). Deviation test construction and power comparison for marked spatial point patterns. Spatial Statistics 11: 19-34. doi: 10.1016/j.spasta.2014.11.004
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#' @examples
#' if(requireNamespace("fda", quietly = TRUE)) {
#'   # Consider ordering of the girls in the Berkeley Growth Study data
#'   # available from the R package fda, see ?growth, according to their
#'   # annual heights or/and changes within years.
#'   # First create sets of curves (vectors), for raw heights and
#'   # for the differences within the years
#'   years <- paste(1:18)
#'   curves <- fda::growth[['hgtf']][years,]
#'   cset1 <- create_curve_set(list(r = as.numeric(years),
#'                                  obs = curves))
#'   plot(cset1, ylab="Height")
#'   cset2 <- create_curve_set(list(r = as.numeric(years[-1]),
#'                                  obs = curves[-1,] - curves[-nrow(curves),]))
#'   plot(cset2)
#'
#'   # Order the girls from most extreme one to the least extreme one, below using the 'area' measure
#'   # a) according to their heights
#'   forder(cset1, measure = 'area')
#'   # Print the 10 most extreme girl indices
#'   order(forder(cset1, measure = 'area'))[1:10]
#'   # b) according to the changes (print indices)
#'   order(forder(cset2, measure = 'area'))[1:10]
#'   # c) simultaneously with respect to heights and changes (print indices)
#'   csets <- list(Height = cset1, Change = cset2)
#'   order(forder(csets, measure = 'area'))[1:10]
#' }
forder <- function(curve_sets, measure = 'erl', scaling = 'qdir',
                   alternative=c("two.sided", "less", "greater"),
                   use_theo = TRUE, probs = c(0.025, 0.975), quantile.type = 7) {
  if(class(curve_sets)[1] == "list") {
    res <- combined_forder(curve_sets,
                           measure = measure, scaling = scaling,
                           alternative = alternative,
                           use_theo = use_theo,
                           probs = probs, quantile.type = quantile.type)
  }
  else {
    res <- individual_forder(curve_sets,
                             measure = measure, scaling = scaling,
                             alternative = alternative,
                             use_theo = use_theo,
                             probs = probs, quantile.type = quantile.type)
  }
  res$distance
}
