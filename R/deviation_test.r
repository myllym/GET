#' Deviation test
#'
#' Crop the curve set to the interval of distances [r_min, r_max],
#' calculate residuals, scale the residuals and perform a deviation test
#' with a chosen deviation measure.
#'
#'
#' The deviation test is based on a test function T(r) and it works as follows:
#'
#' 1) The test function estimated for the data, T_1(r), and for nsim simulations
#' from the null model, T_2(r), ...., T_{nsim+1}(r), must be saved in 'curve_set'
#' and given to the deviation_test function.
#'
#' 2) The deviation_test function then
#'\itemize{
#'   \item Crops the functions to the chosen range of distances [r_min, r_max].
#'   \item If the curve_set does not consist of residuals, i.e. curve_set$is_residual
#'     is FALSE (or does not exists), then the residuals d_i(r) = T_i(r) - T_0(r) are
#'     calculated, where T_0(r) is the expectation of T(r) under the null hypothesis.
#'     If use_theo = TRUE, the theoretical value given in the curve_set$theo is used for
#'     as T_0(r), if it is given. Otherwise, T_0(r) is estimated by the mean of T_j(r),
#'     j=2,...,nsim+1.
#'   \item Scales the residuals. Options are
#'         \itemize{
#'           \item 'none' No scaling. Nothing done.
#'           \item 'q' Quantile scaling.
#'           \item 'qdir' Directional quantile scaling.
#'           \item 'st' Studentised scaling.
#'         }
#'         See for details Myllymäki et al. (2013).
#'   \item Calculates the global deviation measure u_i, i=1,...,nsim+1, see options
#'         for 'measure'.
#'         \itemize{
#'           \item 'max' is the maximum deviation measure
#'
#'           \deqn{u_i = \max_{r \in [r_{\text{min}}, r_{\text{max}}]} | w(r)(T_i(r) - T_0(r))|}{%
#'                 u_i = max_(r in [r_min, r_max]) | w(r)(T_i(r) - T_0(r)) |}
#'
#'           \item 'int2' is the integral deviation measure
#'
#'           \deqn{u_i = \int_{r_{\text{min}}}^{r_{\text{max}}} ( w(r)(T_i(r) - T_0(r)) )^2 dr}{%
#'                 u_i = int_([r_min, r_max]) ( w(r)(T_i(r) - T_0(r)) )^2 dr}
#'
#'           \item 'int' is the 'absolute' integral deviation measure
#'
#'           \deqn{u_i = \int_{r_{\text{min}}}^{r_{\text{max}}} |w(r)(T_i(r) - T_0(r))| dr}{%
#'                 u_i = int_([r_min, r_max]) | w(r)(T_i(r) - T_0(r)) | dr}
#'
#'         }
#'   \item Calculates the p-value.
#'}
#'
#' Currently, there is no special way to take care of the same values of T_i(r)
#' occuring possibly for small distances. Thus, it is preferable to exclude from
#' the test the very small distances r for which ties occur.
#'
#'
#' @inheritParams crop_curves
#' @inheritParams residual
#' @inheritParams deviation
#' @param scaling The name of the scaling to use. Options include 'none',
#'   'q', 'qdir' and 'st'. 'qdir' is default.
#' @param savedevs Logical. Should the global rank values k_i, i=1,...,nsim+1 be returned? Default: FALSE.
#' @return If 'savedevs=FALSE' (default), the p-value is returned.
#' If 'savedevs=TRUE', then a list containing the p-value and calculated deviation measures
#' u_i, i=1,...,nsim+1 (where u_1 corresponds to the data pattern) is returned.
#' @references
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2013). Deviation test construction and power comparison for marked spatial point patterns. arXiv:1306.1028
#' @export
#' @examples
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' require(spatstat)
#' pp <- unmark(spruces)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
#' # The deviation test using the integral deviation measure
#' res <- deviation_test(env, measure = 'int')
#' res
#' # or
#' res <- deviation_test(env, r_min=0, r_max=7, measure='int2')
#'
#' ## Random labeling test
#' #----------------------
#' mpp <- spruces
#' # T(r) = \hat{L}_m(r), an estimator of the L_m(r) function
#' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=999, r_min=1.5, r_max=9.5)
#' res <- deviation_test(curve_set, measure='int2')
#' res
deviation_test <- function(curve_set, r_min = NULL, r_max = NULL,
        use_theo = TRUE, scaling = 'qdir',
        measure = 'max', savedevs=FALSE, ...) {
    curve_set <- crop_curves(curve_set, r_min = r_min, r_max = r_max)
    curve_set <- residual(curve_set, use_theo = use_theo)
    curve_set <- scale_curves(curve_set, scaling = scaling)
    devs <- deviation(curve_set, measure = measure)
    p <- estimate_p_value(devs, ...)
    if(savedevs) res <- list(p=p, devs=devs, call=match.call())
    else res <- list(p=p, call=match.call())
    class(res) <- 'deviation_test'
    res
}

#' Print method for the class 'deviation_test'
#' @usage \method{print}{deviation_test}(x, ...)
#'
#' @param x an 'deviation_test' object
#' @param ... Ignored.
#'
#' @method print deviation_test
#' @export
print.deviation_test <- function(x, ...) {
    with(x, cat("p-value of the test: ", p, "\n", sep=""))
}
