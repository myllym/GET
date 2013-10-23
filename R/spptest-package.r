#' Spatial point process testing
#'
#' This spptest library provides envelope and deviation tests for spatial point processes.
#' The main motivation for this package are the scalings for deviation tests and global
#' envelope tests.
#'
#'
#' A test is performed by three steps as follows. Examples are provided in the test functions.
#'
#' 1) First, make nsim simulations from the null model and estimate your chosen
#' test function T(r) for your data (T_1(r)) and for each simulation (T_2(r),
#' ..., T_{nsim+1}(r)).
#'
#' For testing simple hypothesis, this can be done as follows
#' \itemize{
#'    \item Complete spatial randomness (CSR):
#'          \itemize{
#'            \item Use \code{\link[spatstat]{envelope}}.
#'                  Important: use the option 'savefuns=TRUE'.
#'                  See the help documentation in \code{\link[spatstat]{spatstat}} for possible test functions.
#'          }
#'
#'    \item Random labeling: use \code{\link{random_labelling}} (requires R library marksummary).
#' }
#'
#' For testing the goodness-of-fit of parametric models,
#' (i) you can utilize \code{\link[spatstat]{spatstat}},
#' \itemize{
#'    \item Fit the model to your data by means of the function
#'          \code{\link[spatstat]{ppm}} or \code{\link[spatstat]{kppm}}.
#'          See the help documentation for possible models.
#'    \item Use \code{\link[spatstat]{envelope}} to create nsim simulations
#'          from the fitted model and to calculate T_1(r), T_2(r), ..., T_{nsim+1}(r).
#'          Important: use the option 'savefuns=TRUE'.
#'          See the help documentation in \code{\link[spatstat]{spatstat}} for possible test functions.
#' }
#' or (ii) use your own programs
#' \itemize{
#'    \item Fit the model and create nsim simulations from the fitted model.
#'    \item Calculate T_1(r), T_2(r), ..., T_{nsim+1}(r).
#'    \item Use \code{\link{create_curve_set}} to create a curve_set object
#'          from the estimated functions T_i(r), i=1,...,s+1.
#' }
#' See \code{\link{rank_envelope}} for examples.
#'
#' 2) Optional: modify the curve set T_1(r), T_2(r), ..., T_{nsim+1}(r) for the test.
#'
#' (i) Choose the interval of distances [r_min, r_max] by \code{\link{crop_curves}}.
#'
#' (ii) For better visualisation, take T(r)-T_0(r) by \code{\link{residual}}.
#' T_0(r) is the expectation of T(r) under the null hypothesis.
#'
#' The function \code{\link{envelope_to_curve_set}} can be used to create a curve_set object
#' from the object returned by \code{\link[spatstat]{envelope}}. The "envelope" object can also
#' directly be given to the functions \code{\link{crop_curves}} and \code{\link{residual}}.
#'
#' 3) Perform the test. The alternatives provided in this library are
#' the following:
#'
#' Envelope tests:
#' \itemize{
#'   \item \code{\link{rank_envelope}}, the completely non-parametric rank envelope test
#'   \item \code{\link{st_envelope}}, the studentised envelope test, protected against unequal variance of T(r) for different distances r
#'   \item \code{\link{qdir_envelope}}, the directional quantile envelope test, protected against unequal variance and asymmetry of T(r) for different distances r
#'   \item \code{\link{normal_envelope}}, a parametric envelope test, which uses a normal approximation of T(r)
#' }
#' Note that the recommended minimum number of simulations for the rank
#' envelope test is nsim=4999, while, for the studentised and directional
#' quantile envelope tests, it is nsim=999. (For the normal test, see its documentation.)
#'
#' Deviation tests (no graphical interpretation):
#' \itemize{
#'   \item \code{\link{deviation_test}}
#' }
#'
#' @author
#' Mari Myllymäki (mari.myllymaki@@aalto.fi, mari.j.myllymaki@@gmail.com),
#' Henri Seijo (henri.seijo@@aalto.fi),
#' Tomáš Mrkvička (mrkvicka.toma@@gmail.com),
#' Pavel Grabarnik (gpya@@rambler.ru)
#'
#' @references
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2013). Deviation test construction and power comparison for marked spatial point patterns. arXiv:1306.1028
#'
#' Myllymäki, M., Mrkvička, T., Seijo, H. and Grabarnik, P. (2013). Global envelope tests for spatial point patterns. arXiv:1307.0239 [stat.ME]
#'
#' @name spptest
#' @docType package
#' @aliases spptest-package spptest
NULL
