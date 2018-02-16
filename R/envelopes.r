#' Central region / Global envelope
#'
#' Provides central regions or global envelopes or confidence bands
#'
#'
#' Given a \code{curve_set} (see \code{\link{create_curve_set}} for how to create such an object)
#' or an \code{\link[spatstat]{envelope}} object, the function \code{central_region}
#' construcst a central region, i.e. a global envelope, from the given set of functions (or vectors).
#' There are two options for the functions that the \code{curve_set} can contain:
#' \itemize{
#'  \item If the component \code{obs} of the \code{curve_set} is a matrix,
#' then it is assumed that all the functions are data/observed. In this case,
#' the component \code{sim_m} of the \code{curve_set} (which can be then NULL)
#' is ignored and the central region constructed from the functions given in \code{obs}.
#'  \item If the component \code{obs} is a vector, then \code{sim_m} should be provided as well
#' and it is assumed to contain simulated functions (obtained, e.g., from some model or by permutation).
#' Then the central region is constructed from the functions given in \code{sim_m}.
#' }
#' Thus the \code{curve_set} contains functions (or vectors)
#' \eqn{T_1(r),\dots,T_s(r)}{T_1(r),...,T_s(r)}.
#' In the case of one observed function only,
#' the data function is considered to be \eqn{T_1(r)}{T_1(r)}.
#'
#' Generally an envelope is a band bounded by the vectors (or functions) 
#' \eqn{T_{\text{low}}}{T_lo} and \eqn{T_{\text{hi}}}{T_hi}.
#' A \eqn{100(1-\alpha)}{100(1-alpha)}\% or 100*coverage\% global envelope is a set
#' \eqn{(T_{\text{low}}, T_{\text{hi}})}{(T_lo, T_hi)} of envelope vectors
#' such that the probability that \eqn{T_i}{T_i} falls outside this envelope
#' in any of the d points of the vector \eqn{T_i}{T_i} is less or equal to \eqn{\alpha}{alpha}.
#'
#' The central regions, i.e. global envelopes, can be constructed based on different measures
#' that order the functions from the most extreme one to the least extreme one.
#' The measure can be chosen with the argument \code{type} and the options are given
#' in the following with further description of the global envelope related to the measure:
#' \itemize{
#'  \item \code{'rank'}: extreme rank. Then the global envelope is the global rank envelope
#' proposed by Myllymäki et al. (2017).
#'  \item \code{'erl'}: extreme rank length (Myllymäki et al.,2017).
#' Then the envelope is the global rank envelope based on the extreme rank length ordering.
#' This envelope is constructed as the convex hull of the functions which have extreme rank
#' length measure \eqn{R_i^{\text{erl}}}{Rerl_i}
#' that is larger or equal to the critical \eqn{\alpha}{alpha} level of the extreme rank
#' length measure (Mrkvička et al., 2018).
#'  \item \code{'qdir'}: the directional quantile maximum absolute deviation (MAD) measure
#' (Myllymäki et al., 2015).
#' The directional quantile envelope test (Myllymäki et al., 2017),
#' which takes into account the unequal variances of the test function T(r) for
#' different distances r and is also protected against asymmetry of T(r).
#'  \item \code{'st'}: the studentized MAD measure (Myllymäki et al., 2015).
#' The studentised envelope test (Myllymäki et al., 2017), which takes into account the unequal
#' variances of the test function T(r) for different distances r.
#'  \item \code{'unscaled'}: the unscaled MAD measure. The unscaled envelope test (Ripley, 1981),
#' which leads to envelopes with constant width. It corresponds to the classical
#' maximum deviation test without scaling. This test suffers from unequal variance
#' of T(r) over the distances r and from the asymmetry of distribution of T(r).
#' We recommend to use the other alternatives instead. This unscaled global envelope is
#' provided for reference.
#' }
#' We note that the global envelopes \code{'rank'} and \code{'erl'}
#' are completely non-parametric tests and thus protected against the unequal variances and asymmetry
#' mentioned above.
#' See detailed description of all the measures in \code{\link{forder}}.
#'
#' For each curve in the curve_set, both the data curve and the simulations,
#' an above mention measure k is determined. If savedevs = TRUE, then the
#' measure values \eqn{k_1, k_2, ..., k_{s+1}}{k_1, k_2, ..., k_(s+1)} are
#' returned in the component 'k', where k[1] is the value for the data.
#' Based on the chosen measure, the central region, i.e. the global envelope, is constructed
#' on the chosen interval of argument values (the functions in the \code{curve_set} are assumed
#' to be given on this interval only).
#'
#' @references
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' Mrkvička, T., Hahn, U. and Myllymäki, M. (2018). A one-way ANOVA test for functional data with graphical interpretation. arXiv:1612.03608 [stat.ME]
#'
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2015). Deviation test construction and power comparison for marked spatial point patterns. Spatial Statistics 11: 19-34. doi: 10.1016/j.spasta.2014.11.004
#'
#' Ripley, B.D. (1981). Spatial statistics. Wiley, New Jersey.
#'
#' @param curve_set A curve_set (see \code{\link{create_curve_set}}) or
#' an \code{\link[spatstat]{envelope}} object. If an envelope object is given,
#' it must contain the summary functions from the simulated patterns which can be
#' achieved by setting savefuns = TRUE when calling \code{\link[spatstat]{envelope}}.
#' @param coverage A number between 0 and 1. The 100*coverage\% central region will be calculated.
#' @param savedevs Logical. Should the measure values k_i, i=1,...,s, be returned? Default: FALSE.
#' @param alternative A character string specifying the alternative hypothesis.
#' Must be one of the following: "two.sided" (default), "less" or "greater".
#' The last two options only available for \code{type = 'rank'} and \code{type = 'erl'}.
#' @param type The type of the global envelope with current options for 'rank', 'erl',
#' 'qdir', 'st' and 'unscaled'. See details.
#' @param ... Additional parameters to be passed to \code{\link[stats]{quantile}}
#' (through \code{\link{forder}}) in the case \code{type = 'qdir'}.
#' @param probs A two-element vector containing the lower and upper
#'   quantiles for the 'qdir' envelope, in that order and on the interval [0, 1].
#'   The default values are 0.025 and 0.975, suggested by Myllymäki et al. (2015, 2017).
#' @return An object of class "envelope" and "fv" (see \code{\link[spatstat]{fv.object}}),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item r = the vector of values of the argument r at which the test was made
#' \item obs = the data function, if there is only one data function. Otherwise not existing.
#' \item lo = the lower envelope based on the simulated functions
#' \item hi = the upper envelope based on the simulated functions
#' \item central = If the curve_set (or envelope object) contains a component 'theo',
#'       then this function is used as the central curve and returned in this component.
#'       Otherwise, the central_curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'       Used for visualization only.
#' }
#' Additionally, the return value has attributes
#' \itemize{
#'   \item method = The name of the envelope test ("Rank envelope test" for the rank envelope test)
#'   \item alternative = The alternative specified in the function call.
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#'   \item p_interval = The p-value interval [p_liberal, p_conservative].
#'   \item ties = As the argument \code{ties}.
#'   \item k_alpha = The value of k corresponding to the 100(1-alpha)\% global envelope.
#'   \item k = Global rank values (for type='rank') or extreme rank lengths (for type='erl').
#'   k[1] is the value for the data pattern. Returned only if savedevs = TRUE.
#'   \item call = The call of the function.
#' }
#' and a punch of attributes for the "fv" object type, see \code{\link[spatstat]{fv}}.
#' @export
central_region <- function(curve_set, coverage=0.95, savedevs=FALSE,
                           alternative=c("two.sided", "less", "greater"),
                           type="rank", probs = c(0.025, 0.975), ...) {
  if(coverage < 0 | coverage > 1) stop("Unreasonable value of alpha.")
  alpha <- 1 - coverage
  if(!is.logical(savedevs)) cat("savedevs should be logical. Using the default FALSE.")
  alternative <- match.arg(alternative)
  if(type %in% c("qdir", "st", "unscaled") && alternative != "two.sided") {
    warning("For qdir, st and unscaled envelopes only the two.sided alternative is valid.\n")
    alternative <- "two.sided"
  }

  picked_attr <- pick_attributes(curve_set, alternative=alternative) # saving for attributes / plotting purposes
  curve_set <- convert_envelope(curve_set)

  if(!(type %in% c("rank", "erl", "qdir", "st", "unscaled")))
    stop("No such a type for global envelope.\n")

  # Measures for functional ordering
  measure <- type
  scaling <- ""
  switch(type,
         qdir = {
           measure <- "max"
           scaling <- "qdir"
         },
         st = {
           measure <- "max"
           scaling = "st"
         },
         unscaled = {
           measure <- "max"
           scaling = "none"
         })
  distance <- forder(curve_set, measure = measure, scaling = scaling, alternative = alternative, ...)

  data_and_sim_curves <- data_and_sim_curves(curve_set) # all the functions
  Nfunc <- length(distance) # Number of functions
  nr <- length(curve_set[['r']])
  # Define the central curve T_0
  T_0 <- get_T_0(curve_set)

  #-- Global envelopes
  switch(type,
         rank = {
           #-- the 100(1-alpha)% global rank envelope
           distancesorted <- sort(distance, decreasing=TRUE)
           kalpha <- distancesorted[floor((1-alpha)*(Nfunc))]
           LB <- array(0, nr)
           UB <- array(0, nr)
           for(i in 1:nr){
             Hod <- sort(data_and_sim_curves[,i])
             LB[i]<- Hod[kalpha]
             UB[i]<- Hod[Nfunc-kalpha+1]
           }
         },
         erl = {
           #-- the 100(1-alpha)% global ERL envelope
           distance_lexo_sorted <- sort(distance, decreasing=TRUE)
           kalpha <- distance_lexo_sorted[floor((1-alpha)*Nfunc)]
           curves_for_envelope <- data_and_sim_curves[which(distance >= kalpha),]
           LB <- apply(curves_for_envelope, MARGIN=2, FUN=min)
           UB <- apply(curves_for_envelope, MARGIN=2, FUN=max)
         },
         qdir = {
           curve_set <- residual(curve_set, use_theo = TRUE)
           quant_m <- curve_set_quant(curve_set, probs = probs, ...)
           #-- the 100(1-alpha)% global directional quantile envelope
           distancesorted <- sort(distance)
           kalpha <- distancesorted[floor((1-alpha)*Nfunc)]
           LB <- T_0 - kalpha*abs(quant_m[1,])
           UB <- T_0 + kalpha*abs(quant_m[2,])
         },
         st = {
           sdX <- curve_set_sd(curve_set)
           #-- calculate the 100(1-alpha)% global studentized envelope
           distancesorted <- sort(distance)
           kalpha <- distancesorted[floor((1-alpha)*Nfunc)]
           LB <- T_0 - kalpha*sdX
           UB <- T_0 + kalpha*sdX
         },
         unscaled = {
           #-- calculate the 100(1-alpha)% global unscaled envelope
           distancesorted <- sort(distance)
           kalpha <- distancesorted[floor((1-alpha)*Nfunc)]
           LB <- T_0 - kalpha
           UB <- T_0 + kalpha
         })

  switch(alternative,
         "two.sided" = {},
         "less" = { UB <- Inf },
         "greater" = { LB <- -Inf })

  if(is.vector(curve_set[['obs']]))
    res <- structure(data.frame(r=curve_set[['r']], obs=curve_set[['obs']], central=T_0, lo=LB, hi=UB),
                     class = c("envelope", "fv", "data.frame"))
  else
    res <- structure(data.frame(r=curve_set[['r']], central=T_0, lo=LB, hi=UB),
                     class = c("envelope", "fv", "data.frame"))
  attr(res, "method") <- "Global envelope"
  attr(res, "type") <- type
  attr(res, "alternative") <- alternative
  attr(res, "critical_alpha") <- kalpha
  if(savedevs) attr(res, "measure") <- distance
  # for fv
  attr(res, "fname") <- picked_attr$fname
  attr(res, "argu") <- "r"
  attr(res, "valu") <- "central"
  attr(res, "ylab") <- picked_attr$ylab
  attr(res, "fmla") <- ". ~ r"
  attr(res, "alim") <- c(min(curve_set[['r']]), max(curve_set[['r']]))
  attr(res, "labl") <- picked_attr$labl
  attr(res, "desc") <- picked_attr$desc
  #attr(res, "unitname") <- "unit / units"
  attr(res, "shade") <- c("lo", "hi")
  attr(res, "call") <- match.call()
  res
}

#' The rank envelope test
#'
#' The rank envelope test, p-values and global envelopes
#'
#'
#' The rank envelope test is a completely non-parametric test, which provides
#' the 100(1-alpha)\% global envelope for the chosen test function T(r) on
#' the chosen interval of distances and associated p-values.
#'
#' Given a \code{curve_set} (see \code{\link{create_curve_set}} for how to create such an object)
#' or an \code{\link[spatstat]{envelope}} object,
#' which contains both the data curve (or function or vector) \eqn{T_1(r)}{T_1(r)} and
#' the simulated curves \eqn{T_2(r),\dots,T_{s+1}(r)}{T_2(r),...,T_(s+1)(r)},
#' the test is carried out as described in the following sections that describe
#' ordering of the functions, p-values and global envelopes.
#'
#' @section Ranking of the curves:
#' First the curves in the curve set are ranked from the most extreme one to the least extreme one
#' either by using the extreme ranks R_i and/or the extreme rank lengths \eqn{R_i^{\text{erl}}}{Rerl_i}.
#' The extreme rank is defined as the minimum of pointwise ranks of the curve \eqn{T_i(r)}{T_i(r)},
#' where the pointwise rank is the rank of the value of the curve for a specific r-value among the
#' corresponding values of the s other curves such that the lowest ranks correspond to the most extreme
#' values of the curves. How the pointwise ranks are determined exactly depends on the whether a
#' one-sided (\code{alternative} is "less" or "greater") or the two-sided test (\code{alternative="two.sided"}) is
#' chosen, for details see Mrkvička et al. (2017, page 1241) or Mrkvička et al. (2018, page 6).
#' The extreme ranks can contain many ties, for which reason Myllymäki et al. (2017) proposed the
#' extreme rank length ordering. Considering the vector of pointwise ordered ranks
#' \eqn{\mathbf{R}_i}{RP_i} of the ith curve, the extreme rank length measure is equal to
#' \deqn{R_i^{\text{erl}} = \frac{1}{s+1}\sum_{j=1}^{s+1} \1(\mathbf{R}_j \prec \mathbf{R}_i)}{Rerl_i = \sum_{j=1}^{s} 1(RP_j "<" RP_i) / (s + 1)}
#' where \eqn{\mathbf{R}_j \prec \mathbf{R}_i}{RP_j "<" RP_i} if and only if
#' there exists \eqn{n\leq d}{n<=d} such that for the first k, \eqn{k<n}{k<n}, pointwise ordered
#' ranks of \eqn{\mathbf{R}_j}{RP_j} and \eqn{\mathbf{R}_i}{RP_i} are equal and the n'th rank of
#' \eqn{\mathbf{R}_j}{RP_j} is smaller than that of \eqn{\mathbf{R}_i}{RP_i}.
#'
#' For each curve in the curve_set, both the data curve and the simulations,
#' an above mention measure k is determined. If savedevs = TRUE, then the
#' measure values \eqn{k_1, k_2, ..., k_{s+1}}{k_1, k_2, ..., k_(s+1)} are
#' returned in the component 'k', where k[1] is the value for the data.
#'
#' @section P-values:
#' In the case \code{type="rank"}, based on the extreme ranks k_i, i=1, ..., s+1,
#' the p-interval is calculated. Because the extreme ranks contain ties, there is not just
#' one p-value. The p-interval is given by the most liberal and the most conservative p-value
#' estimate. This interval is by default plotted for the object returned by the rank_envelope function.
#' Also a single p-value is calculated and returned in the attribute \code{p}.
#' By default this single p-value is the mid-rank p-value, but another option can be used by
#' specifying \code{ties} argument.
#'
#' If the case \code{type="erl"}, the (single) p-value based on the extreme rank length ordering
#' of the functions is calculated and returned in the attribute \code{p}.
#'
#' @section Global envelope:
#' The 100(1-alpha)\% global envelope is provided in addition to the p-values. 
#' If \code{type="rank"} then the envelope is the global rank envelope proposed by
#' Myllymäki et al. (2017).
#' If \code{type="erl"} then the envelope is the global rank envelope based on the
#' extreme rank length ordering. This envelope is constructed as the convex hull of
#' the functions which have extreme rank length measure \eqn{R_i^{\text{erl}}}{Rerl_i}
#' that is larger or equal to the critical \eqn{\alpha}{alpha} level of the extreme rank
#' length measure (Mrkvička et al., 2018).
#'
#' @section Number of simulations:
#' The extreme rank length ordering test (\code{type="erl"}) allows in principle a lower numbe
#' of simulations to be used than the test based on extreme ranks (\code{type="rank"}).
#' However, we recommend some thousands of simulations in any case to achieve a good power
#' and repeatability of the test.
#'
#' @references
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' Mrkvička, T., Myllymäki, M. and Hahn, U. (2017). Multiple Monte Carlo testing, with applications in spatial point processes. Statistics & Computing 27 (5): 1239-1255. doi: 10.1007/s11222-016-9683-9
#'
#' Mrkvička, T., Hahn, U. and Myllymäki, M. (2018). A one-way ANOVA test for functional data with graphical interpretation. arXiv:1612.03608 [stat.ME]
#'
#' @param curve_set A curve_set (see \code{\link{create_curve_set}}) or an \code{\link[spatstat]{envelope}}
#'  object. If an envelope object is given, it must contain the summary
#'  functions from the simulated patterns which can be achieved by setting
#'  savefuns = TRUE when calling \code{\link[spatstat]{envelope}}.
#' @param alpha The significance level. The 100(1-alpha)\% global envelope will be calculated.
#' @param savedevs Logical. Should the global rank values k_i, i=1,...,nsim+1 be returned? Default: FALSE.
#' @param alternative A character string specifying the alternative hypothesis. Must be one of the following:
#'         "two.sided" (default), "less" or "greater".
#' @param type The type of the global envelope with current options for "rank" and "erl".
#' If "rank", the global rank envelope accompanied by the p-interval is given (Myllymäki et al., 2017).
#' If "erl", the global rank envelope based on extreme rank lengths accompanied by the extreme rank
#' length p-value is given (Myllymäki et al., 2017, Mrkvicka et al., 2018). See details and additional
#' sections thereafter.
#' @param lexo Obsolete. Use type instead.
#' @param ties The method to obtain a unique p-value when type = "rank".
#' Possible values are 'midrank', 'random', 'conservative', 'liberal' and 'erl'.
#' For 'conservative' the resulting p-value will be the highest possible.
#' For 'liberal' the p-value will be the lowest possible.
#' For 'random' the rank of the obs within the tied values is uniformly sampled so that the resulting
#' p-value is at most the conservative option and at least the liberal option.
#' For 'midrank' the mid-rank within the tied values is taken.
#' For 'erl' the extreme rank length p-value is calculated.
#' The default is 'midrank'.
#' @return An object of class "envelope_test", "envelope" and "fv" (see \code{\link[spatstat]{fv.object}}),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item r = the vector of values of the argument r at which the test was made
#' \item obs = values of the test function for the data point pattern
#' \item lo = the lower envelope based on the simulated functions
#' \item hi = the upper envelope based on the simulated functions
#' \item central = If the curve_set (or envelope object) contains a component 'theo',
#'       then this function is used as the central curve and returned in this component.
#'       Otherwise, the central_curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'       Used for visualization only.
#' }
#' Additionally, the return value has attributes
#' \itemize{
#'   \item method = The name of the envelope test ("Rank envelope test" for the rank envelope test)
#'   \item alternative = The alternative specified in the function call.
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#'   \item p_interval = The p-value interval [p_liberal, p_conservative].
#'   \item ties = As the argument \code{ties}.
#'   \item k_alpha = The value of k corresponding to the 100(1-alpha)\% global envelope.
#'   \item k = Global rank values (for type="rank") or extreme rank lengths (for type="erl").
#'   k[1] is the value for the data pattern. Returned only if savedevs = TRUE.
#'   \item call = The call of the function.
#' }
#' and a punch of attributes for the "fv" object type, see \code{\link[spatstat]{fv}}.
#' @export
#' @seealso \code{\link{random_labelling}}, \code{\link{plot.envelope_test}}
#' @examples
#'
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' require(spatstat)
#' pp <- unmark(spruces)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=2499, savefuns=TRUE, correction="translate")
#' # The rank envelope test
#' res <- rank_envelope(env)
#' # Plot the result.
#' # - The central curve is now obtained from env[['theo']], which is the
#' # value of the L-function under the null hypothesis (L(r) = r).
#' plot(res)
#' # or (requires R library ggplot2)
#' plot(res, plot_style="ggplot2")
#'
#' ## Advanced use:
#' # Choose the interval of distances [r_min, r_max] (at the same time create a curve_set from 'env')
#' curve_set <- crop_curves(env, r_min = 1, r_max = 7)
#' # For better visualisation, take the L(r)-r function
#' curve_set <- residual(curve_set, use_theo = TRUE)
#' # Do the rank envelope test
#' res <- rank_envelope(curve_set); plot(res, plot_style="ggplot2")
#'
#' ## Random labeling test
#' #----------------------
#' # requires library 'marksummary'
#' mpp <- spruces
#' # 1) Perform simulations under the random labelling hypothesis and calculate
#' # the test function T(r) for the data pattern (mpp) and each simulation.
#' # The command below specifies that the test function is T(r) = \hat{L}_m(r),
#' # which is an estimator of the mark-weighted L function, L_m(r),
#' # with translational edge correction (default).
#' # The random_labelling function returns the centred functions \hat{L}_m(r)-T_0(r),
#' # where T_0(r) = \hat{L}(r) is the unmarked L function.
#' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=2499, r_min=1.5, r_max=9.5)
#' # 2) Do the rank envelope test
#' res <- rank_envelope(curve_set)
#' # 3) Plot the test result
#' plot(res, plot_style="ggplot2", ylab=expression(italic(L[m](r)-L(r))))
#'
#' # Make the test using instead the test function T(r) = \hat{L}_mm(r);
#' # which is an estimator of the mark-weighted L function, L_mm(r),
#' # with translational edge correction (default).
#' curve_set <- random_labelling(mpp, mtf_name = 'mm', nsim=2499, r_min=1.5, r_max=9.5)
#' res <- rank_envelope(curve_set)
#' plot(res, plot_style="ggplot2", ylab=expression(italic(L[mm](r)-L(r))))
#'
#' ## Goodness-of-fit test (typically conservative, see dg.global_envelope for adjusted tests)
#' #-----------------------------------------------
#' pp <- unmark(spruces)
#' # Minimum distance between points in the pattern
#' min(nndist(pp))
#' # Fit a model
#' fittedmodel <- ppm(pp, interaction=Hardcore(hc=1)) # Hardcore process
#'
#' \dontrun{
#' # Simulating Gibbs process by 'envelope' is slow, because it uses the MCMC algorithm
#' #env <- envelope(fittedmodel, fun="Jest", nsim=999, savefuns=TRUE,
#'                  correction="none", r=seq(0, 4, length=500))
#'
#' # Using direct algorihm can be faster, because the perfect simulation is used here.
#' simulations <- NULL
#' for(j in 1:2499) {
#'    simulations[[j]] <- rHardcore(beta=exp(fittedmodel$coef[1]),
#'                                  R = fittedmodel$interaction$par$hc,
#'                                  W = pp$window);
#'    if(j%%10==0) cat(j, "...", sep="")
#' }
#' env <- envelope(pp, simulate=simulations, fun="Jest", nsim=length(simulations),
#'                 savefuns=TRUE, correction="none", r=seq(0, 4, length=500))
#' curve_set <- crop_curves(env, r_min = 1, r_max = 3.5)
#' res <- rank_envelope(curve_set); plot(res, plot_style="ggplot2")
#' }
#'
#' ## A test based on a low dimensional random vector
#' #-------------------------------------------------
#' # Let us generate some example data.
#' X <- matrix(c(-1.6,1.6),1,2) # data pattern X=(X_1,X_2)
#' if(requireNamespace("mvtnorm", quietly = TRUE)) {
#'   Y <- mvtnorm::rmvnorm(200,c(0,0),matrix(c(1,0.5,0.5,1),2,2)) # simulations
#'   plot(Y, xlim=c(min(X[,1],Y[,1]), max(X[,1],Y[,1])), ylim=c(min(X[,2],Y[,2]), max(X[,2],Y[,2])))
#'   points(X, col=2)
#'
#'   # Test the null hypothesis is that X is from the distribution of Y's (or if it is an outlier).
#'
#'   # Case 1. The test vector is (X_1, X_2)
#'   cset1 <- create_curve_set(list(r=1:2, obs=as.vector(X), sim_m=t(Y)))
#'   res1 <- rank_envelope(cset1)
#'   plot(res1)
#'
#'   # Case 2. The test vector is (X_1, X_2, (X_1-mean(Y_1))*(X_2-mean(Y_2))).
#'   t3 <- function(x, y) { (x[,1]-mean(y[,1]))*(x[,2]-mean(y[,2])) }
#'   cset2 <- create_curve_set(list(r=1:3, obs=c(X[,1],X[,2],t3(X,Y)), sim_m=rbind(t(Y), t3(Y,Y))))
#'   res2 <- rank_envelope(cset2)
#'   plot(res2)
#' }
rank_envelope <- function(curve_set, alpha=0.05, savedevs=FALSE,
                          alternative=c("two.sided", "less", "greater"),
                          type="rank", lexo=NULL, ties) {
    if(alpha < 0 | alpha > 1) stop("Unreasonable value of alpha.")
    if(!is.logical(savedevs)) cat("savedevs should be logical. Using the default FALSE.")
    alternative <- match.arg(alternative)

    picked_attr <- pick_attributes(curve_set, alternative=alternative) # saving for attributes / plotting purposes
    curve_set <- convert_envelope(curve_set)

    if(!(type %in% c("rank", "erl"))) stop("No such a type for global envelope.\n")
    if(is.logical(lexo) && lexo) type <- "erl"
    # The type of the p-value
    if(missing(ties)) ties <- p_value_ties_default()
    possible_ties <- c('midrank', 'random', 'conservative', 'liberal', 'erl')
    if(!(ties %in% possible_ties)) stop("Unreasonable ties argument!\n")

    # data_curve = the vector of test function values for data
    # sim_curves = matrix where each row contains test function values of a simulation under null hypothesis
    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    nr <- length(curve_set$r)
    # Define the central curve T_0
    T_0 <- get_T_0(curve_set)

    data_and_sim_curves <- rbind(data_curve, sim_curves)
    loranks <- apply(data_and_sim_curves, MARGIN=2, FUN=rank, ties.method = "average")
    hiranks <- Nsim+2-loranks
    # k:
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

    #-- the ERL p-value
    if(type == "erl" | ties == "erl") { # rank the curves by lexical ordering
        # order ranks within each curve
        sortranks <- apply(allranks, 1, sort) # curves now represented as columns
        lexo_values <- do.call("order", split(sortranks, row(sortranks))) # indices! of the functions from the most extreme to least extreme one
        newranks <- 1:(Nsim+1)
        distance_lexo <- newranks[order(lexo_values)] # ordering of the functions by the extreme rank counts
        #-- calculate the p-value
        u_lexo <- -distance_lexo
        p <- estimate_p_value(x=u_lexo[1], sim_vec=u_lexo[-1], ties="conservative")
    }
    #-- the p-interval (based on extreme ranks) and global envelopes
    switch(type,
           rank = {
             distance <- apply(allranks, MARGIN=1, FUN=min) # extreme ranks R_i
             u <- -distance
             #-- p-interval
             p_low <- estimate_p_value(x=u[1], sim_vec=u[-1], ties='liberal')
             p_upp <- estimate_p_value(x=u[1], sim_vec=u[-1], ties='conservative')
             if(ties != "erl") p <- estimate_p_value(x=u[1], sim_vec=u[-1], ties=ties) # Note: case ties=="erl" calculated above
             #-- the 100(1-alpha)% global rank envelope
             distancesorted <- sort(distance, decreasing=TRUE)
             kalpha <- distancesorted[floor((1-alpha)*(Nsim+1))]
             LB <- array(0, nr);
             UB <- array(0, nr);
             for(i in 1:nr){
               Hod <- sort(data_and_sim_curves[,i])
               LB[i]<- Hod[kalpha];
               UB[i]<- Hod[Nsim+1-kalpha+1];
             }
           },
           erl = {
             #-- the 100(1-alpha)% global ERL envelope
             distance_lexo_sorted <- sort(distance_lexo, decreasing=TRUE)
             kalpha_lexo <- distance_lexo_sorted[floor((1-alpha)*(Nsim+1))]
             curves_for_envelope <- data_and_sim_curves[which(distance_lexo >= kalpha_lexo),]
             LB <- apply(curves_for_envelope, MARGIN=2, FUN=min)
             UB <- apply(curves_for_envelope, MARGIN=2, FUN=max)
           })

    switch(alternative,
            "two.sided" = {},
            "less" = { UB <- Inf },
            "greater" = { LB <- -Inf })

    res <- structure(data.frame(r=curve_set[['r']], obs=data_curve, central=T_0, lo=LB, hi=UB),
                     class = c("envelope_test", "envelope", "fv", "data.frame"))
    attr(res, "method") <- "Rank envelope test"
    attr(res, "alternative") <- alternative
    attr(res, "p") <- p
    switch(type,
           rank = {
             attr(res, "k_alpha") <- kalpha
             attr(res, "k") <- distance
             attr(res, "p_interval") <- c(p_low, p_upp)
             attr(res, "ties") <- ties
           },
           erl = {
             attr(res, "k_alpha") <- kalpha_lexo
             attr(res, "k") <- distance_lexo / (Nsim+1)
             attr(res, "p_interval") <- NULL
             attr(res, "ties") <- "extreme rank length"
           })
    # for fv
    attr(res, "fname") <- picked_attr$fname
    attr(res, "argu") <- "r"
    attr(res, "valu") <- "obs"
    attr(res, "ylab") <- picked_attr$ylab
    attr(res, "fmla") <- ". ~ r"
    attr(res, "alim") <- c(min(curve_set[['r']]), max(curve_set[['r']]))
    attr(res, "labl") <- picked_attr$labl
    attr(res, "desc") <- picked_attr$desc
    #attr(res, "unitname") <- "unit / units"
    attr(res, "shade") <- c("lo", "hi")
    attr(res, "call") <- match.call()
    res
}

#' Print method for the class 'envelope_test'
#' @usage \method{print}{envelope_test}(x, ...)
#'
#' @param x an 'envelope_test' object
#' @param ... Ignored.
#'
#' @method print envelope_test
#' @export
print.envelope_test <- function(x, ...) {
    cat(attr(x, "method"), "\n",
        " p-value of the test: ", attr(x, "p"), sep="")
    if(!is.null(attr(x, "ties"))) cat(" (ties method: ", attr(x, "ties"), ")\n", sep="")
    else cat("\n")
    if(!is.null(attr(x, "p_interval")))
        cat(" p-interval         : (", attr(x, "p_interval")[1], ", ", attr(x, "p_interval")[2],")\n", sep="")
}

#' Plot method for the class 'envelope_test'
#' @usage \method{plot}{envelope_test}(x, plot_style="basic", base_size=15, dotplot=length(x$r)<10,
#'                                      main, ylim, xlab, ylab, use_ggplot2, ...)
#'
#' @param x an 'envelope_test' object
#' @param plot_style One of the following "basic", "fv" or "ggplot2".
#' The option "basic" (default) offers a very basic global envelope plot.
#' The option "fv" utilizes the plot routines of the function value table \code{\link[spatstat]{fv.object}}.
#' For "ggplot2", a plot with a coloured envelope ribbon is provided. Requires R library ggplot2.
#' The option "fv" is currently only available for tests with one test function, whereas the other true allow
#' also tests with several tests functions.
#' @param base_size Base font size, to be passed to theme style when \code{plot_style = "ggplot2"}.
#' @param dotplot Logical. If TRUE, then instead of envelopes a dot plot is done.
#' Suitable for low dimensional test vectors. Only applicable if \code{plot_style} is "basic".
#' Default: TRUE if the dimension is less than 10, FALSE otherwise.
#' @param main See \code{\link{plot.default}}. A sensible default exists.
#' @param ylim See \code{\link{plot.default}}. A sensible default exists.
#' @param xlab See \code{\link{plot.default}}. A sensible default exists.
#' @param ylab See \code{\link{plot.default}}. A sensible default exists.
#' @param use_ggplot2 Logical, whether plot_style is "ggplot2" or not. Outdated, use the argument plot_style instead.
#' @param ... Additional parameters to be passed to \code{\link{env_basic_plot}}, \code{\link{dotplot}}
#' (if dotplot=TRUE) or \code{\link{env_ggplot}} (if plot_style="ggplot2").
#'
#' @method plot envelope_test
#' @export
#' @seealso \code{\link{rank_envelope}}, \code{\link{st_envelope}}, \code{\link{qdir_envelope}}
plot.envelope_test <- function(x, plot_style="basic", base_size=15, dotplot=length(x$r)<10,
        main, ylim, xlab, ylab, use_ggplot2, ...) {
    if(!missing(use_ggplot2) && is.logical(use_ggplot2) && use_ggplot2) plot_style <- "ggplot2"
    else use_ggplot2 <- FALSE

    if(missing('main')) main <- env_main_default(x)
    if(missing('ylim')) ylim <- env_ylim_default(x, use_ggplot2)
    if(missing('xlab')) xlab <- expression(italic(r))
    if(missing('ylab')) ylab <- expression(italic(T(r)))

    plot_style <- spatstat::pickoption("ptype", plot_style, c(basic = "basic",
                                                            b = "basic",
                                                            fv = "fv",
                                                            f = "fv",
                                                            ggplot2 = "ggplot2",
                                                            ggplot = "ggplot2",
                                                            g = "ggplot2"))

    switch(plot_style,
           basic = {
             if(dotplot) {
               env_dotplot(x, main, ylim, xlab, ylab, ...)
             }
             else {
               env_basic_plot(x, main, ylim, xlab, ylab, ...)
             }
           },
           fv = {
             spatstat::plot.fv(x, main=main, ylim=ylim, ...)
           },
           ggplot2 = {
             env_ggplot(x, base_size, main, ylim, xlab, ylab, ...)
           })
}

#' Studentised envelope test
#'
#' The studentised envelope test, which takes into account the unequal
#' variances of the test function T(r) for different distances r.
#'
#'
#' @references
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2015). Deviation test construction and power comparison for marked spatial point patterns. Spatial Statistics 11: 19-34. doi: 10.1016/j.spasta.2014.11.004
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' @inheritParams rank_envelope
#' @param savedevs Logical. Should the deviation values u_i, i=1,...,nsim+1 be returned? Default: FALSE.
#' @return An object of class "envelope_test", "envelope" and "fv" (see \code{\link[spatstat]{fv.object}}),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item r = the vector of values of the argument r at which the test was made
#' \item obs = values of the test function for the data point pattern
#' \item lo = the lower envelope based on the simulated functions
#' \item hi = the upper envelope based on the simulated functions
#' \item central = If the curve_set (or envelope object) contains a component 'theo',
#'       then this function is used as the central curve and returned in this component.
#'       Otherwise, the central curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'       Used for visualization only.
#' }
#' Additionally, the return value has attributes
#' \itemize{
#'   \item method = The name of the envelope test ("Studentised envelope test" for the studentised envelope test)
#'   \item alternative = "two-sided
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#'   \item u_alpha = The value of u corresponding to the 100(1-alpha)\% global envelope.
#'   \item u = Deviation values (u[1] is the value for the data pattern). Returned only if savedevs = TRUE.
#'   \item call = The call of the function.
#' }
#' and a punch of attributes for the "fv" object type.
#' @export
#' @importFrom stats sd
#' @examples
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' require(spatstat)
#' pp <- spruces
#' ## Test for complete spatial randomness (CSR)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
#' # The studentised envelope test
#' res <- st_envelope(env)
#' plot(res)
#' # or (requires R library ggplot2)
#' plot(res, plot_style="ggplot2")
#'
#' ## Advanced use:
#' # Create a curve set, choosing the interval of distances [r_min, r_max]
#' curve_set <- crop_curves(env, r_min = 1, r_max = 8)
#' # For better visualisation, take the L(r)-r function
#' curve_set <- residual(curve_set, use_theo = TRUE)
#' # The studentised envelope test
#' res <- st_envelope(curve_set); plot(res, plot_style="ggplot2")
#'
#' ## Random labeling test
#' #----------------------
#' # requires library 'marksummary'
#' mpp <- spruces
#' # Use the test function T(r) = \hat{L}_m(r), an estimator of the L_m(r) function
#' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=2499, r_min=1.5, r_max=9.5)
#' res <- st_envelope(curve_set)
#' plot(res, plot_style="ggplot2", ylab=expression(italic(L[m](r)-L(r))))
st_envelope <- function(curve_set, alpha=0.05, savedevs=FALSE) {

    picked_attr <- pick_attributes(curve_set, alternative="two.sided")
    curve_set <- convert_envelope(curve_set)

    if(alpha < 0 | alpha > 1) stop("Unreasonable value of alpha.")
    if(!is.logical(savedevs)) cat("savedevs should be logical. Using the default FALSE.")

    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    nr <- dim(sim_curves)[2]

    # Define T_0 and residual curve_set
    T_0 <- get_T_0(curve_set)
    curve_set <- residual(curve_set, use_theo = TRUE)

    sdX <- as.vector(apply(curve_set[['sim_m']], MARGIN=1, FUN=stats::sd))

    # Calculate deviation measures
    distance <- array(0, Nsim+1);
    scaled_curve_set <- weigh_curves(curve_set, divisor_to_coeff(sdX))
    #devs <- deviation(scaled_curve_set, measure = 'max', scaling='qdir')
    # u_1
    distance[1] <- max(abs(scaled_curve_set$obs))
    # u_2, ..., u_{s+1}
    distance[2:(Nsim+1)] <- apply(abs(scaled_curve_set[['sim_m']]), 2, max)

    #-- calculate the p-value
    p <- estimate_p_value(x=distance[1], sim_vec=distance[-1])

    #-- calculate the 100(1-alpha)% global envelope
    distancesorted <- sort(distance);
    talpha <- distancesorted[floor((1-alpha)*(Nsim+1))];
    LB <- T_0 - talpha*sdX;
    UB <- T_0 + talpha*sdX;

    res <- structure(data.frame(r=curve_set[['r']], obs=data_curve, central=T_0, lo=LB, hi=UB),
                     class = c("envelope_test", "envelope", "fv", "data.frame"))
    attr(res, "method") <- "Studentised envelope test"
    attr(res, "alternative") <- "two.sided"
    attr(res, "p") <- p
    attr(res, "u_alpha") <- talpha
    if(savedevs) attr(res, "u") <- distance
    # for fv
    attr(res, "fname") <- picked_attr$fname
    attr(res, "argu") <- "r"
    attr(res, "valu") <- "obs"
    attr(res, "ylab") <- picked_attr$ylab
    attr(res, "fmla") <- ". ~ r"
    attr(res, "alim") <- c(min(curve_set[['r']]), max(curve_set[['r']]))
    attr(res, "labl") <- picked_attr$labl
    attr(res, "desc") <- picked_attr$desc
    attr(res, "shade") <- c("lo", "hi")
    attr(res, "call") <- match.call()
    res
}

#' Directional quantile envelope test
#'
#' The directional quantile envelope test, which takes into account the unequal 
#' variances of the test function T(r) for different distances r and is also 
#' protected against asymmetry of T(r).
#'
#' @references
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2015). Deviation test construction and power comparison for marked spatial point patterns. Spatial Statistics 11: 19-34. doi: 10.1016/j.spasta.2014.11.004
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' @inheritParams st_envelope
#' @param probs A two-element vector containing the lower and upper
#'   quantiles for the envelope, in that order and on the interval [0, 1].
#'   The default values are 0.025 and 0.975.
#' @return An object of class "envelope_test", "envelope" and "fv" (see \code{\link[spatstat]{fv.object}}),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item r = the vector of values of the argument r at which the test was made
#' \item obs = values of the test function for the data point pattern
#' \item lo = the lower envelope based on the simulated functions
#' \item hi = the upper envelope based on the simulated functions
#' \item central = If the curve_set (or envelope object) contains a component 'theo',
#'       then this function is used as the central curve and returned in this component.
#'       Otherwise, the central curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'       Used for visualization only.
#' }
#' Additionally, the return value has attributes
#' \itemize{
#'   \item method = The name of the envelope test ("Directional quantile envelope test" for the directional quantile envelope test)
#'   \item alternative = "two-sided
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#'   \item u_alpha = The value of u corresponding to the 100(1-alpha)\% global envelope.
#'   \item u = Deviation values (u[1] is the value for the data pattern). Returned only if savedevs = TRUE.
#'   \item call = The call of the function.
#' }
#' and a punch of attributes for the "fv" object type.
#' @export
#' @importFrom stats quantile
#' @examples
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' require(spatstat)
#' pp <- spruces
#' ## Test for complete spatial randomness (CSR)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
#' # The directional quantile envelope test
#' res <- qdir_envelope(env)
#' plot(res)
#' # or (requires R library ggplot2)
#' plot(res, plot_style="ggplot2")
#'
#' ## Advanced use:
#' # Create a curve set, choosing the interval of distances [r_min, r_max]
#' curve_set <- crop_curves(env, r_min = 1, r_max = 8)
#' # For better visualisation, take the L(r)-r function
#' curve_set <- residual(curve_set, use_theo = TRUE)
#' # The directional quantile envelope test
#' res <- qdir_envelope(curve_set); plot(res, plot_style="ggplot2")
#'
#' ## Random labeling test
#' #----------------------
#' # requires library 'marksummary'
#' mpp <- spruces
#' # Use the test function T(r) = \hat{L}_m(r), an estimator of the L_m(r) function
#' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=2499, r_min=1.5, r_max=9.5)
#' res <- qdir_envelope(curve_set)
#' plot(res, plot_style="ggplot2", ylab=expression(italic(L[m](r)-L(r))))
qdir_envelope <- function(curve_set, alpha=0.05, savedevs=FALSE, probs = c(0.025, 0.975)) {

    picked_attr <- pick_attributes(curve_set, alternative="two.sided")
    curve_set <- convert_envelope(curve_set)
    check_probs(probs)

    if(alpha < 0 | alpha > 1) stop("Unreasonable value of alpha.")
    if(!is.logical(savedevs)) cat("savedevs should be logical. Using the default FALSE.")

    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    nr <- dim(sim_curves)[2]

    # Define T_0 and residual curve_set
    T_0 <- get_T_0(curve_set)
    curve_set <- residual(curve_set, use_theo = TRUE)

    # calculate quantiles for residual curve_set (i.e. for sim_curves - T_0)
    quant_m <- apply(curve_set[['sim_m']], 1, stats::quantile, probs = probs)
    abs_coeff <- divisor_to_coeff(abs(quant_m))
    lower_coeff <- abs_coeff[1, , drop = TRUE]
    upper_coeff <- abs_coeff[2, , drop = TRUE]

    # Calculate deviation measures
    distance <- array(0, Nsim+1);
    # u_1
    scaled_residuals <- weigh_both_sides(curve_set[['obs']], upper_coeff, lower_coeff)
    distance[1] <- max(abs(scaled_residuals))
    # u_2, ..., u_{s+1}
    sim_scaled_residuals <- weigh_both_sides(curve_set[['sim_m']], upper_coeff, lower_coeff)
    distance[2:(Nsim+1)] <- apply(abs(sim_scaled_residuals), 2, max)

    #-- calculate the p-value
    p <- estimate_p_value(x=distance[1], sim_vec=distance[-1])

    #-- calculate the 100(1-alpha)% global envelope
    distancesorted <- sort(distance)
    talpha <- distancesorted[floor((1-alpha)*(Nsim+1))]
    LB <- T_0 - talpha*abs(quant_m[1,])
    UB <- T_0 + talpha*abs(quant_m[2,])

    res <- structure(data.frame(r=curve_set[['r']], obs=data_curve, central=T_0, lo=LB, hi=UB),
                     class = c("envelope_test", "envelope", "fv", "data.frame"))
    attr(res, "method") <- "Directional quantile envelope test"
    attr(res, "alternative") <- "two.sided"
    attr(res, "p") <- p
    attr(res, "u_alpha") <- talpha
    if(savedevs) attr(res, "u") <- distance
    # for fv
    attr(res, "fname") <- picked_attr$fname
    attr(res, "argu") <- "r"
    attr(res, "valu") <- "obs"
    attr(res, "ylab") <- picked_attr$ylab
    attr(res, "fmla") <- ". ~ r"
    attr(res, "alim") <- c(min(curve_set[['r']]), max(curve_set[['r']]))
    attr(res, "labl") <- picked_attr$labl
    attr(res, "desc") <- picked_attr$desc
    attr(res, "shade") <- c("lo", "hi")
    attr(res, "call") <- match.call()
    res
}

#' Unscaled envelope test
#'
#' The unscaled envelope test, which leads to envelopes with constant width
#' over the distances r. It corresponds to the classical maximum deviation test
#' without scaling.
#'
#'
#' This test suffers from unequal variance of T(r) over the distances r and from
#' the asymmetry of distribution of T(r). We recommend to use the rank_envelope
#' (if number of simulations close to 5000 can be afforded) or st_envelope/qdir_envelope
#' (if large number of simulations cannot be afforded) instead.
#'
#' @references
#' Ripley, B.D. (1981). Spatial statistics. Wiley, New Jersey.
#'
#' @inheritParams st_envelope
#' @return An object of class "envelope_test", "envelope" and "fv" (see \code{\link[spatstat]{fv.object}}),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item r = the vector of values of the argument r at which the test was made
#' \item obs = values of the test function for the data point pattern
#' \item lo = the lower envelope based on the simulated functions
#' \item hi = the upper envelope based on the simulated functions
#' \item central = If the curve_set (or envelope object) contains a component 'theo',
#'       then this function is used as the central curve and returned in this component.
#'       Otherwise, the central curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'       Used for visualization only.
#' }
#' Additionally, the return value has attributes
#' \itemize{
#'   \item method = The name of the envelope test ("Studentised envelope test" for the studentised envelope test)
#'   \item alternative = "two-sided
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#'   \item u_alpha = The value of u corresponding to the 100(1-alpha)\% global envelope.
#'   \item u = Deviation values (u[1] is the value for the data pattern). Returned only if savedevs = TRUE.
#'   \item call = The call of the function.
#' }
#' and a punch of attributes for the "fv" object type.
#' @export
#' @examples
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' require(spatstat)
#' pp <- spruces
#' ## Test for complete spatial randomness (CSR)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
#' # The studentised envelope test
#' res <- unscaled_envelope(env)
#' plot(res)
#' # or (requires R library ggplot2)
#' plot(res, plot_style="ggplot2")
#'
#' ## Advanced use:
#' # Create a curve set, choosing the interval of distances [r_min, r_max]
#' curve_set <- crop_curves(env, r_min = 1, r_max = 8)
#' # For better visualisation, take the L(r)-r function
#' curve_set <- residual(curve_set, use_theo = TRUE)
#' # The studentised envelope test
#' res <- unscaled_envelope(curve_set); plot(res, plot_style="ggplot2")
#'
#' ## Random labeling test
#' #----------------------
#' # requires library 'marksummary'
#' mpp <- spruces
#' # Use the test function T(r) = \hat{L}_m(r), an estimator of the L_m(r) function
#' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=2499, r_min=1.5, r_max=9.5)
#' res <- unscaled_envelope(curve_set)
#' plot(res, plot_style="ggplot2", ylab=expression(italic(L[m](r)-L(r))))
unscaled_envelope <- function(curve_set, alpha=0.05, savedevs=FALSE) {

    picked_attr <- pick_attributes(curve_set, alternative="two.sided")
    curve_set <- convert_envelope(curve_set)

    if(alpha < 0 | alpha > 1) stop("Unreasonable value of alpha.")
    if(!is.logical(savedevs)) cat("savedevs should be logical. Using the default FALSE.")

    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    nr <- dim(sim_curves)[2]

    # Define T_0 and residual curve_set
    T_0 <- get_T_0(curve_set)
    curve_set <- residual(curve_set, use_theo = TRUE)

    # Calculate deviation measures
    distance <- array(0, Nsim+1);
    # u_1
    distance[1] <- max(abs(curve_set$obs))
    # u_2, ..., u_{s+1}
    distance[2:(Nsim+1)] <- apply(abs(curve_set[['sim_m']]), 2, max)

    #-- calculate the p-value
    p <- estimate_p_value(x=distance[1], sim_vec=distance[-1])

    #-- calculate the 100(1-alpha)% global envelope
    distancesorted <- sort(distance);
    talpha <- distancesorted[floor((1-alpha)*(Nsim+1))];
    LB <- T_0 - talpha;
    UB <- T_0 + talpha;

    res <- structure(data.frame(r=curve_set[['r']], obs=data_curve, central=T_0, lo=LB, hi=UB),
                     class = c("envelope_test", "envelope", "fv", "data.frame"))
    attr(res, "method") <- "Unscaled envelope test"
    attr(res, "alternative") <- "two.sided"
    attr(res, "p") <- p
    attr(res, "u_alpha") <- talpha
    if(savedevs) attr(res, "u") <- distance
    # for fv
    attr(res, "fname") <- picked_attr$fname
    attr(res, "argu") <- "r"
    attr(res, "valu") <- "obs"
    attr(res, "ylab") <- picked_attr$ylab
    attr(res, "fmla") <- ". ~ r"
    attr(res, "alim") <- c(min(curve_set[['r']]), max(curve_set[['r']]))
    attr(res, "labl") <- picked_attr$labl
    attr(res, "desc") <- picked_attr$desc
    attr(res, "shade") <- c("lo", "hi")
    attr(res, "call") <- match.call()
    res
}
