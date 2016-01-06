#' Saplings data set
#'
#' Saplings data set
#'
#'
#' A pattern of small trees (height <= 15 m) originating from an uneven aged multi-species
#' broadleaf nonmanaged forest in Kaluzhskie Zaseki, Russia.
#'
#' The pattern is a sample part of data collected over 10 ha plot as a part of a research
#' program headed by project leader Prof. O.V. Smirnova.
#'
#' @format An object of class \code{\link[spatstat]{ppp.object}} representing the point
#' pattern of tree locations.
#'
#' @usage data(saplings)
#' @references
#' Grabarnik, P. and Chiu, S. N. (2002) Goodness-of-fit test for complete spatial randomness against
#' mixtures of regular and clustered spatial point processes. \emph{Biometrika}, \bold{89}, 411–421.
#'
#' van Lieshout, M.-C. (2010) Spatial point process theory. In Handbook of Spatial Statistics (eds. A. E.
#' Gelfand, P. J. Diggle, M. Fuentes and P. Guttorp), Handbooks of Modern Statistical Methods. Boca
#' Raton: CRC Press.
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2015). Global envelope tests for
#' spatial point patterns. arXiv:1307.0239v4 [stat.ME]
#'
#' @keywords datasets
#' @keywords spatial
#' @name saplings
#' @docType data
#' @examples
#' # This is an example analysis of the saplings data set, see Myllymäki et al. (2015).
#' require(spatstat)
#' data(saplings)
#'
#' # First choose the r-distances for L (r) and J (rJ) functions, respectively.
#' nr <- 500
#' rmin <- 0.3; rminJ <- 0.3
#' rmax <- 10; rmaxJ <- 6;
#' rstep <- (rmax-rmin)/nr; rstepJ <- (rmaxJ-rminJ)/nr;
#' r <- c(0, seq(rmin, rmax, by=rstep))
#' rJ <- c(0, seq(rminJ, rmaxJ, by=rstepJ))
#'
#' #-- CSR test --# (a simple hypothesis)
#' #--------------#
#' # First, a CSR test using the L(r)-r function:
#' # Note: CSR is simulated by fixing the number of points and generating nsim simulations
#' # from the binomial process, i.e. we deal with a simple hypothesis.
#' env <- envelope(saplings, nsim=2499,
#'  simulate=expression(runifpoint(saplings$n, win=saplings$window)), # Simulate CSR
#'  fun="Lest", correction="translate", # T(r) = estimator of L with translational edge correction
#'  transform = expression(.-r),        # Take the L(r)-r function instead of L(r)
#'  r=r,                                # Specify the distance vector
#'  savefuns=TRUE)                      # Save the estimated functions
#' # Crop the curves to the interval of distances [rmin, rmax]
#' # (at the same time create a curve_set from 'env')
#' curve_set <- crop_curves(env, r_min = rmin, r_max = rmax)
#' # Perform the rank envelope test
#' res <- rank_envelope(curve_set)
#' # Plot the result.
#' plot(res, use_ggplot2=TRUE, ylab=expression(italic(hat(L)(r)-r)))
#'
#' # -> The CSR hypothesis is clearly rejected and the rank envelope indicates clear
#' # clustering of saplings. Next we explore the Matern cluster process as a null model.
#'
#' \dontrun{
#' #-- Testing the Matern cluster process --# (a composite hypothesis)
#' #----------------------------------------#
#' # Fit the Matern cluster process to the pattern (using minimum contrast estimation with the pair
#' # correction function)
#' fitted_model <- kppm(saplings, clusters = "MatClust", statistic="pcf")
#' summary(fitted_model)
#'
#' nsim <- 499 # 19 for experimenting; 499 ok for test = 'qdir'
#'
#' # Make the adjusted directional quantile global envelope test using the L(r)-r function
#' # (For the rank envelope test, choose test = "rank" instead and increase nsim.)
#' system.time( # timing; takes a lot of time,
#'              # if nsim is reasonably large (for 'qdir' & nsim=499, about 1,6 h)
#'   adjenvL <- dg.global_envelope(X = fitted_model,
#'                                fun="Lest", correction="translate",
#'                                transform = expression(.-r), r=r,
#'                                test = "qdir", nsim = nsim, nsimsub = nsim,
#'                                r_min=rmin, r_max=rmax,
#'                                save.cons.envelope=TRUE)
#' )
#' # Plot the test result
#' # As save.cons.envelope was set above to TRUE, the unadjusted global envelope can be plotted
#' # together with the adjusted envelope.
#' plot(adjenvL, use_ggplot2=TRUE, ylab=expression(italic(L(r)-r)), plot_unadjusted=TRUE)
#' # Setting plot_unadjusted=FALSE only the adjusted envelope is plotted.
#'
#' # From the test with the L(r)-r function, it appears that the Matern cluster model would be
#' # a reasonable model for the saplings pattern.
#' # To further explore the goodness-of-fit of the Matérn cluster process, test the
#' # model with the J function:
#' system.time( # timing;takes a lot of time if nsim is reasonably large
#'   adjenvJ <- dg.global_envelope(X = fitted_model,
#'                                 fun="Jest", correction="none", r=rJ,
#'                                 test = test, nsim = nsim, nsimsub = nsim,
#'                                 r_min=rminJ, r_max=rmaxJ,
#'                                 save.cons.envelope=TRUE)
#' )
#' # Plot the test result
#' plot(adjenvJ, ylab=expression(italic(J(r))), plot_unadjusted=FALSE)
#' # -> the Matérn cluster process not adequate for the saplings data
#' }
NULL
