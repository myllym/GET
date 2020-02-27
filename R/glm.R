# Check that formula.reduced is nested within formula.full and includes something extra.
check_isnested <- function(formula.full, formula.reduced) {
  # Check that the reduced model is nested within the full model
  if(!all(labels(terms(formula.reduced)) %in% labels(terms(formula.full)))) stop("The reduced model includes some extra variables, not in the full model.\n")
  if(attr(terms(formula.full), "intercept") < attr(terms(formula.reduced), "intercept")) stop("The reduced model includes intercept, but the full model does not.\n")
  # Check that the full model includes something in addition to the reduced model
  if(all(labels(terms(formula.full)) %in% labels(terms(formula.reduced)))) stop("The full model should not be equal to the reduced model.\n")
}

# Preliminary checks for the graph.flm and frank.flm
#' @importFrom stats terms
flm.checks <- function(nsim, formula.full, formula.reduced, curve_sets, factors = NULL, fast = TRUE) {
  # Preliminary checks
  vars <- all.vars(formula.full)
  vars.reduced <- all.vars(formula.reduced)
  # Check that the reduced model is nested within the full model
  if(!all(labels(terms(formula.reduced)) %in% labels(terms(formula.full)))) stop("The reduced model includes some extra variables, not in the full model.\n")
  if(attr(terms(formula.full), "intercept") < attr(terms(formula.reduced), "intercept")) stop("The reduced model includes intercept, but the full model does not.\n")
  # Check that the full model includes something in addition to the reduced model
  if(all(labels(terms(formula.full)) %in% labels(terms(formula.reduced)))) stop("The full model should not be equal to the reduced model.\n")
  if(nsim < 1) stop("Not a reasonable value of nsim.\n")
  if(vars[1] != "Y") stop("The formula should be off the form Y ~ .... where Y is the response.\n")
  if(class(curve_sets)[1] != "list") {
    curve_sets <- list(Y=curve_sets)
  }
  available <- unique(c(names(curve_sets), names(factors)))
  if(!all(vars %in% available)) stop("The variables in the formula not found in the given data (curve_sets and factors).\n")
  if(!all(sapply(curve_sets, function(x) inherits(x, c("curve_set", "fdata"))))) stop("The components of curve_sets do not have a valid class.\n")
  if(inherits(curve_sets[['Y']], "fdata")) {
    einfo <- curve_sets[['Y']]$names
  } else einfo <- NULL
  curve_sets <- lapply(curve_sets, convert_fdata)
  if(!all(sapply(curve_sets, function(x) is.matrix(x[['obs']])))) stop("The curve_set must include data functions (sim_m ignored).\n")
  curve_sets <- check_curve_set_dimensions(curve_sets)
  # Put Y and factors into data.l
  data.l <- list()
  data.l[['Y']] <- t(curve_sets[['Y']][['obs']]) # -> each row corresponds to a data function
  # The argument values
  r <- curve_sets[['Y']][['r']]
  Nfunc <- nrow(data.l[['Y']]) # Number of functions
  nr <- ncol(data.l[['Y']]) # Number of argument values
  vars.csets <- vars[vars %in% names(curve_sets)]
  factors_in_curvesets <- !fast
  if(length(curve_sets) > 1 & length(vars.csets) > 1) { # Factors provided in the curve_sets
    for(i in 2:length(vars.csets)) data.l[[vars.csets[i]]] <- t(curve_sets[[vars.csets[i]]][['obs']])
    factors_in_curvesets <- TRUE
  }
  vars.factors <- vars[vars %in% names(factors)]
  if(!is.null(factors) & length(vars.factors) > 0) {
    if(class(factors)[1] != "data.frame") stop("Invalid factors argument.\n")
    if(nrow(factors) != Nfunc) stop("The dimensions of Y and factors do not match.\n")
    # Expand the factors to each argument value
    for(i in 1:length(vars.factors)) {
      data.l[[vars.factors[i]]] <- if(factors_in_curvesets) {
        matrix(factors[,vars.factors[i]], nrow=nrow(factors), ncol=nr)
      } else { factors[,vars.factors[i]] }
    }
  }
  dfs <- if(factors_in_curvesets) {
    lapply(1:nr, FUN = function(i) { as.data.frame(lapply(data.l, FUN = function(x) x[,i])) })
  } else {
    list(as.data.frame(data.l[-1]))
  }

  list(Y=data.l[['Y']], dfs=dfs, r=r, Nfunc=Nfunc, nr=nr, einfo=einfo)
}

# Regress the given data (true or permuted) against the full model and
# get an effect of interest at all r values in a matrix.
# @param Y True or permuted values of Y inserted to dfs.
# @param dfs A list containing data (Y and factors), all variables in formula.full
# nameinteresting = names of the interesting coefficients
#' @importFrom stats lm
#' @importFrom stats dummy.coef
genCoefmeans.m <- function(Y, dfs, formula.full, nameinteresting, ...) {
  # If covariates are constant with respect to location
  if(length(dfs) == 1) return(genCoefmeans.mlm(Y, dfs[[1]], formula.full, nameinteresting, ...))
  nr <- ncol(Y)
  effect.a <- sapply(1:nr, function(i) {
    df <- dfs[[i]]
    df$Y <- Y[,i]
    M.full <- stats::lm(formula.full, data=df, ...)
    allcoef <- stats::dummy.coef(M.full)
    unlist(allcoef[nameinteresting])
  })
  # Special case for when there is only 1 coefficient of interest
  if(is.vector(effect.a)) {
    name <- names(effect.a)[1]
    dim(effect.a) <- c(length(effect.a), 1)
    dimnames(effect.a) = list(NULL, name)
    effect.a
  } else t(effect.a)
}

# If covariates are constant with respect to location,
# use a multiple linear model (instead of multiple linear models)
#' @importFrom stats lm
#' @importFrom stats coef
#' @importFrom stats dummy.coef
genCoefmeans.mlm <- function(Y, df, formula.full, nameinteresting, ...) {
  df$Y <- Y
  M.full <- stats::lm(formula.full, data=df, ...)
  if(!length(M.full$xlevels)) { # Only continuous factors in the model
    effect.a <- t(coef(M.full)[nameinteresting,, drop=FALSE])
  } else {
    allcoef <- stats::dummy.coef(M.full)
    effect.a <- do.call(cbind, allcoef[nameinteresting])
    factornames <- names(unlist(lapply(allcoef[nameinteresting], function(x) if(is.matrix(x)) x[1,] else x[1])))
    dimnames(effect.a) <- list(NULL, factornames)
  }
  effect.a
}

# The constrasts of the interesting effects in a matrix
#' @importFrom stats lm
#' @importFrom stats dummy.coef
genCoefcontrasts.m <- function(Y, dfs, formula.full, nameinteresting, ...) {
  nr <- ncol(Y)
  effect.a <- genCoefmeans.m(Y=Y, dfs=dfs, formula.full=formula.full, nameinteresting=nameinteresting, ...)
  nameinteresting <- colnames(effect.a) # including levels of factors
  k <- length(nameinteresting)
  if(k == 1) stop("The option \'contrasts\' only valid for discrete factors with at least two levels.\n")
  # contrasts
  cont <- matrix(0, nrow=nr, ncol=k*(k-1)/2)
  cont.names <- vector(length=k*(k-1)/2)
  counter <- 1
  for(i in 1:(k-1)) for(j in (i+1):k)  {
    cont[,counter] <- effect.a[,i] - effect.a[,j] # coef_i - coef_j
    cont.names[counter] <- paste(nameinteresting[i], nameinteresting[j], sep="-")
    counter <- counter + 1
  }
  colnames(cont) <- cont.names
  cont
}


# General F-values from lm-model using lm (slow but ... arguments allowed)
# @param Y True or permuted values of Y inserted to dfs.
# @param dfs The list of (factors) data at each argument value.
#' @importFrom stats lm
#' @importFrom stats anova
genFvaluesLM <- function(Y, dfs, formula.full, formula.reduced, ...) {
  nr <- length(dfs)
  if(nr == 1) return(genFvaluesMLM(Y, dfs[[1]], formula.full, formula.reduced, ...))
  Fvalues <- vector(length=nr)
  for(i in 1:nr) {
    df <- dfs[[i]]
    df$Y <- Y[,i]
    M.full <- stats::lm(formula = formula.full, data = df, ...)
    M.reduced <- stats::lm(formula = formula.reduced, data = df, ...)
    Anova.res <- stats::anova(M.reduced, M.full)
    Fvalues[i] <- Anova.res$F[2]
  }
  Fvalues
}

# If covariates are constant with respect to location,
# use a multiple linear model (instead of multiple linear models)
#' @importFrom stats lm
#' @importFrom stats deviance
#' @importFrom stats df.residual
genFvaluesMLM <- function(Y, df, formula.full, formula.reduced, ...) {
  df$Y <- Y
  M.full <- lm(formula = formula.full, data = df, ...)
  M.reduced <- lm(formula = formula.reduced, data = df, ...)
  # F-statistic
  (deviance(M.full)-deviance(M.reduced))/(df.residual(M.full)-df.residual(M.reduced))/(deviance(M.full)/df.residual(M.full))
}

# Parameter estimate b for lm
bcoef <- function(Y, X) {
  solve(a = t(X)%*%X, b = t(X)%*%Y)
}

# Y = observed data
# X1 = design matrix of the full model
# X2 = design matrix of the reduced model
# This F is the same as obtained by (but faster):
# M.full <- stats::lm(formula = formula.full, data = df, ...)
# M.reduced <- stats::lm(formula = formula.reduced, data = df, ...)
# Anova.res <- stats::anova(M.reduced, M.full); Anova.res$F[2]
Fvalue <- function(Y, X1, X2) {
  # Parameter estimates b
  b1 <- bcoef(Y, X1)
  b2 <- bcoef(Y, X2)
  # Errors
  e1 <- c(Y - X1%*%b1)
  e2 <- c(Y - X2%*%b2)
  # F-statistic
  (length(Y)-length(b1))*(sum(e2^2)-sum(e1^2))/((length(b1)-length(b2))*sum(e1^2))
}

# General F-values from lm-model together with the design matrices
#' @importFrom stats lm
#' @importFrom stats anova
genFvaluesObs <- function(dfs, formula.full, formula.reduced) {
  nr <- length(dfs)
  Fvalues <- vector(length=nr)
  x.full <- x.reduced <- list()
  for(i in 1:nr) {
    df <- dfs[[i]]
    # Call lm to obtain the design matrices (x)
    x.full[[i]] <- stats::lm(formula = formula.full, data = df, x = TRUE)$x
    x.reduced[[i]] <- stats::lm(formula = formula.reduced, data = df, x = TRUE)$x
    Fvalues[i] <- Fvalue(Y = df$Y, X1 = x.full[[i]], X2 = x.reduced[[i]])
  }
  list(Fvalues = Fvalues, full.X = x.full, reduced.X = x.reduced)
}

# General F-value from lm-model
# Y = data, the dependent
# designX.full = design matrix of the full model
# designX.reduced = design matrix of the reduced model
genFvaluesSim <- function(Y, designX.full, designX.reduced) {
  nr <- ncol(Y)
  Fvalues <- vector(length=nr)
  for(i in 1:nr) {
    Fvalues[i] <- Fvalue(Y = Y[,i], X1 = designX.full[[i]], X2 = designX.reduced[[i]])
  }
  Fvalues
}

#' Graphical functional GLM
#'
#' Non-parametric graphical tests of significance in functional general linear model (GLM)
#'
#'
#' The function \code{graph.flm} performs the graphical functional GLM of Mrkvička et al. (2019).
#' This is a nonparametric graphical test of significance of a covariate in functional GLM.
#' The test is able to find not only if the factor of interest is significant, but also which
#' functional domain is responsible for the potential rejection.
#' In the case of functional multi-way main effect ANOVA or functional main effect ANCOVA models,
#' the test is able to find which groups differ (and where they differ).
#' In the case of functional factorial ANOVA or functional factorial ANCOVA models,
#' the test is able to find which combination of levels (which interactions) differ (and where they differ).
#' The described tests are global envelope tests applied in the context of GLMs.
#' The Freedman-Lane algorithm (Freedman and Lane, 1983) is applied to permute the functions
#' (to obtain the simulations under the null hypothesis of "no effects");
#' consequently, the test approximately achieves the desired significance level.
#'
#' The specification of the full and reduced formulas is important. The reduced model should be
#' nested within the reduced model. The full model should include in addition to the reduced
#' model the interesting factors whose effects are under investigation.
#' The implementation to find the coefficients of the interesting factors is based on dummy.coef and
#' the restrictions there apply.
#'
#' There are different versions of the implementation depending on the application.
#' Given that the argument \code{fast} is TRUE, then
#' \itemize{
#' \item If all the covariates are constant across the functions, i.e. they can be provided in the
#' argument \code{factors}, then a linear model is fitted separately by least-squares estimation to
#' the data at each argument value of the functions fitting a multiple linear model by \code{\link[stats]{lm}}.
#' The possible extra arguments passed in \code{...} to \code{\link[stats]{lm}} must be of the form that
#' \code{\link[stats]{lm}} accepts for fitting a multiple linear model. In the basic case, no extra arguments are
#' needed.
#' \item If some of the covariates vary across the space and there are user specified extra arguments given in
#' \code{...}, then the implementation fits a linear model at each argument value of the functions using
#' \code{\link[stats]{lm}}, which can be rather slow. The arguments \code{...} are passed to \code{\link[stats]{lm}}
#' for fitting each linear model.
#' }
#' By setting \code{fast = FALSE}, it is possible to use the slow version for any case. Usually this is not desired.
#'
#' @inheritParams graph.fanova
#' @inheritParams GET.composite
#' @param formula.full The formula specifying the general linear model,
#' see \code{formula} in \code{\link[stats]{lm}}.
#' @param formula.reduced The formula of the reduced model with nuisance factors only. This model
#' should be nested within the full model.
#' @param curve_sets A named list of sets of curves giving the dependent variable (Y), and
#' possibly additionally all the factors. The dimensions of the elements should
#' match with each other, i.e. the factor values should be given for each argument value
#' and each function. If factors are given in the argument \code{factors}, then can also be just
#' the curve set representing Y. Also \code{\link[fda.usc]{fdata}} objects allowed.
#' @param factors A data frame of factors. An alternative way to specify factors when they
#' are constant for all argument values. The number of rows of the data frame should be equal
#' to the number of curves. Each column should specify the values of a factor.
#' @param ... Additional arguments to be passed to \code{\link[stats]{lm}}. See details.
#' @param GET.args A named list of additional arguments to be passed to \code{\link{global_envelope_test}}.
#' @param mc.args A named list of additional arguments to be passed to \code{\link{mclapply}}.
#' Only relevant if \code{mc.cores} is more than 1.
#' @param cl Allows parallelization through the use of \code{\link{parLapply}} (works also
#' in Windows), see the argument \code{cl} there, and examples.
#' @param fast Logical. See details.
#' @return A \code{global_envelope} or \code{combined_global_envelope} object,
#' which can be printed and plotted directly.
#' @export
#' @references
#' Mrkvička, T., Roskovec, T. and Rost, M. (2019) A nonparametric graphical tests of significance in functional GLM. arXiv:1902.04926 [stat.ME]
#'
#' Freedman, D., & Lane, D. (1983) A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics, 1(4), 292-298. doi:10.2307/1391660
#' @importFrom stats lm
#' @importFrom stats predict.lm
#' @importFrom stats dummy.coef
#' @importFrom parallel mclapply
#' @importFrom parallel parLapply
#' @importFrom parallel clusterEvalQ
#' @examples
#' data(rimov)
#' \donttest{
#' res <- graph.flm(nsim=19, # Increase the number of simulations for serious analysis!
#'                  formula.full = Y~Year,
#'                  formula.reduced = Y~1,
#'                  curve_sets = list(Y=rimov), factors = data.frame(Year = 1979:2014))
#' }
#' \dontshow{
#' res <- graph.flm(nsim = 3,
#'                  formula.full = Y~Year,
#'                  formula.reduced = Y~1,
#'                  curve_sets = list(Y=rimov), factors = data.frame(Year = 1979:2014),
#'                  GET.args = list(alpha=0.25))
#' }
#' plot(res)
#'
#' # Test if there is a change in the slope in 1994,
#' # i.e. the full model is T = a + b*year + c*year:group,
#' res <- graph.flm(nsim = 19, # Increase the number of simulations for serious analysis!
#'                  formula.full = Y ~ Year + Year:Group,
#'                  formula.reduced = Y ~ Year,
#'                  curve_sets = list(Y=rimov),
#'                  factors = data.frame(Year = 1979:2014,
#'                                      Group = factor(c(rep(1,times=24), rep(2,times=12)),
#'                                                     levels=1:2)),
#'                  contrasts = FALSE)
#' plot(res)
#'
#' \donttest{
#' data(GDPtax)
#' factors.df <- data.frame(Group = GDPtax$Group, Tax = GDPtax$Profittax)
#' res.tax_within_group <- graph.flm(nsim = 999,
#'                                   formula.full = Y~Group+Tax+Group:Tax,
#'                                   formula.reduced = Y~Group+Tax,
#'                                   curve_sets = list(Y=GDPtax$GDP),
#'                                   factors = factors.df)
#' plot(res.tax_within_group)
#' }
graph.flm <- function(nsim, formula.full, formula.reduced, curve_sets, factors = NULL,
                      contrasts = FALSE, savefuns = FALSE, ..., GET.args = NULL,
                      mc.cores = 1L, mc.args = NULL, cl = NULL,
                      fast = TRUE) {
  # Set up the contrasts
  op <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(op))
  if(!is.null(cl)) {
    clusterEvalQ(cl, { op <- options(contrasts = c("contr.sum", "contr.poly")) })
    on.exit({ clusterEvalQ(cl, { options(op) }) }, add=TRUE)
  }
  # Preliminary checks and formulation of the data to suitable form for further processing
  X <- flm.checks(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
                   curve_sets=curve_sets, factors=factors, fast=fast)

  nameinteresting <- setdiff(labels(terms(formula.full)), labels(terms(formula.reduced)))

  # setting that 'fun' is a function
  if(!contrasts) genCoef <- genCoefmeans.m
  else genCoef <- genCoefcontrasts.m

  # Fit the full model to the data and obtain the coefficients
  obs <- genCoef(X$Y, X$dfs, formula.full, nameinteresting, ...)

  #-- Freedman-Lane procedure
  # Fit the reduced model at each argument value to get the fitted values and residuals
  loopfun1 <- function(i, ...) {
    if(length(X$dfs) == 1) {
      df <- X$dfs[[1]]
      df$Y <- X$Y[,i]
    } else df <- X$dfs[[i]]
    mod.red <- lm(formula.reduced, data=df, ...)
    list(fitted.m = mod.red$fitted.values, res.m = mod.red$residuals)
  }
  if(is.null(cl)) mclapply_res <- do.call(mclapply, c(list(X=1:X$nr, FUN=loopfun1, mc.cores=mc.cores), mc.args, ...))
  else mclapply_res <- parLapply(cl, 1:X$nr, loopfun1, ...)
  fitted.m <- sapply(mclapply_res, function(x) x$fitted.m)
  res.m <- sapply(mclapply_res, function(x) x$res.m)

  # Simulations by permuting the residuals + calculate the coefficients
  Yperm <- function() { # Permutation
    permutation <- sample(1:nrow(res.m), size=nrow(res.m), replace=FALSE)
    # Permute the residuals (rows in res.m) and create new 'y'
    fitted.m + res.m[permutation, ]
  }
  loopfun2 <- function(i, ...) {
    # Regress the permuted data against the full model and get a new effect of interest
    genCoef(Yperm(), X$dfs, formula.full, nameinteresting, ...)
  }
  if(is.null(cl)) sim <- do.call(mclapply, c(list(X=1:nsim, FUN=loopfun2, mc.cores=mc.cores), mc.args, ...))
  else sim <- parLapply(cl, 1:nsim, loopfun2, ...)
  sim <- simplify2array(sim)
  complabels <- colnames(obs)

  csets <- NULL
  for(i in 1:ncol(obs)) {
    csets[[complabels[i]]] <- create_curve_set(list(r = X$r,
                                                    obs = obs[,i],
                                                    sim_m = sim[,i,]))
  }
  res <- do.call(global_envelope_test,
                 c(list(curve_sets=csets, alternative="two.sided", nstep=1), GET.args))
  attr(res, "method") <- "Graphical functional GLM" # Change method name
  attr(res, "labels") <- complabels
  # Take the xlab and ylab from a fdata object, if such is given:
  if(!is.null(X$einfo)) {
    if(!is.null(X$einfo$xlab)) {
      if(inherits(res, "global_envelope")) {
        attr(res, "xlab") <- attr(res, "xexp") <- X$einfo$xlab
        if(inherits(X$einfo$xlab, "expression")) attr(res, "xlab") <- deparse(X$einfo$xlab)
      }
      if(inherits(res, "combined_global_envelope")) {
        attr(attr(res, "level2_ge"), "xlab") <- attr(attr(res, "level2_ge"), "xexp") <- X$einfo$xlab
        if(inherits(X$einfo$xlab, "expression")) attr(attr(res, "level2_ge"), "xlab") <- deparse(X$einfo$xlab)
      }
    }
    if(!is.null(X$einfo$ylab)) {
      if(inherits(res, "global_envelope")) {
        attr(res, "ylab") <- attr(res, "yexp") <- X$einfo$ylab
        if(inherits(attr(res, "ylab"), "expression")) attr(res, "ylab") <- deparse(attr(res, "ylab"))
      }
      if(inherits(res, "combined_global_envelope")) {
        attr(attr(res, "level2_ge"), "ylab") <- attr(attr(res, "level2_ge"), "yexp") <- X$einfo$ylab
        if(inherits(X$einfo$ylab, "expression")) attr(attr(res, "level2_ge"), "ylab") <- deparse(X$einfo$ylab)
      }
    }
  }
  attr(res, "call") <- match.call()
  if(savefuns) attr(res, "simfuns") <- csets
  res
}

#' F rank functional GLM
#'
#' Multiple testing in permutation inference for the general linear model (GLM)
#'
#'
#' The function \code{frank.flm} performs
#' a nonparametric test of significance of a covariate in the functional GLM.
#' Similarly as in the graphical functional GLM (\code{\link{graph.flm}}),
#' the Freedman-Lane algorithm (Freedman and Lane, 1983) is applied to permute the functions
#' (to obtain the simulations under the null hypothesis of "no effects");
#' consequently, the test approximately achieves the desired significance level.
#' In contrast to the graphical functional GLM, the F rank functional GLM is based on the
#' F-statistics that are calculated at each argument value of the functions.
#' The global envelope test is applied to the observed and simulated F-statistics.
#' The test is able to find if the factor of interest is significant and also which
#' argument values of the functional domain are responsible for the potential rejection.
#'
#' The specification of the full and reduced formulas is important. The reduced model should be
#' nested within the reduced model. The full model should include in addition to the reduced
#' model the interesting factors whose effects are under investigation. Please avoid use of
#' '*' when specifying interactions, e.g. factor1*factor2; instead explicitly specify all
#' components of the model.
#'
#' There are different versions of the implementation depending on the application.
#' Given that the argument \code{fast} is TRUE, then
#' \itemize{
#' \item If all the covariates are constant across the functions, i.e. they can be provided in the
#' argument \code{factors}, then a linear model is fitted separately by least-squares estimation to
#' the data at each argument value of the functions fitting a multiple linear model by \code{\link[stats]{lm}}.
#' The possible extra arguments passed in \code{...} to \code{\link[stats]{lm}} must be of the form that
#' \code{\link[stats]{lm}} accepts for fitting a multiple linear model. In the basic case, no extra arguments are
#' needed.
#' \item If some of the covariates vary across the space, i.e. they are provided in the list of curve sets in
#' the argument \code{curve_sets} together with the dependent functions, but there are no extra arguments given
#' by the user in \code{...}, there is a rather fast implementation of the F-value calculation (which does not
#' use \code{\link[stats]{lm}}).
#' \item If some of the covariates vary across the space and there are user specified extra arguments given in
#' \code{...}, then the implementation fits a linear model at each argument value of the functions using
#' \code{\link[stats]{lm}}, which can be rather slow. The arguments \code{...} are passed to \code{\link[stats]{lm}}
#' for fitting each linear model.
#' }
#' By setting \code{fast = FALSE}, the latter version is used even in a case where faster implementation would be
#' available. Usually this is not desired.
#'
#' @inheritParams graph.flm
#' @return A \code{global_envelope} object, which can be printed and plotted directly.
#' @export
#' @references
#' Mrkvička, T., Myllymäki, M. and Narisetty, N. N. (2019) New methods for multiple testing in permutation inference for the general linear model. arXiv:1906.09004 [stat.ME]
#'
#' Freedman, D., & Lane, D. (1983) A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics, 1(4), 292-298. doi:10.2307/1391660
#' @examples
#' data(GDPtax)
#' factors.df <- data.frame(Group = GDPtax$Group, Tax = GDPtax$Profittax)
#' \donttest{
#' res.tax_within_group <- frank.flm(nsim = 999,
#'                                   formula.full = Y~Group+Tax+Group:Tax,
#'                                   formula.reduced = Y~Group+Tax,
#'                                   curve_sets = list(Y=GDPtax$GDP),
#'                                   factors = factors.df)
#' }
#' \dontshow{
#' res.tax_within_group <- frank.flm(nsim = 4,
#'                                   formula.full = Y~Group+Tax+Group:Tax,
#'                                   formula.reduced = Y~Group+Tax,
#'                                   curve_sets = list(Y=GDPtax$GDP),
#'                                   factors = factors.df,
#'                                   GET.args = list(alpha=0.20))
#' }
#' plot(res.tax_within_group)
# Freedman-Lane procedure (Freedman and Lane, 1983, p. 385)
#' @importFrom parallel mclapply
#' @importFrom parallel parLapply
#' @importFrom stats lm
frank.flm <- function(nsim, formula.full, formula.reduced, curve_sets, factors = NULL,
                      savefuns = TRUE, ..., GET.args = NULL,
                      mc.cores = 1, mc.args = NULL, cl = NULL,
                      fast = TRUE) {
  # Preliminary checks and formulation of the data to suitable form for further processing
  X <- flm.checks(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
                   curve_sets=curve_sets, factors=factors, fast=fast)

  extraargs <- list(...)
  if(length(extraargs) > 0) fast <- FALSE
  # The fast version is meant for the case where there are covariates that vary across the function.
  # Then length(X$dfs) > 1.
  if(length(X$dfs) == 1) fast <- FALSE
  # Calculate the F-statistic for the data, and, if fast, obtain also the design matrices
  if(fast) obs <- genFvaluesObs(X$dfs, formula.full, formula.reduced)
  else obs <- genFvaluesLM(X$Y, X$dfs, formula.full, formula.reduced, ...)
  # Freedman-Lane procedure
  # Fit the reduced model at each argument value to get fitted values and residuals
  if(fast) {
    loopfun1 <- function(i, ...) { # ... ignored
      b <- bcoef(Y = X$dfs[[i]]$Y, X = obs$reduced.X[[i]])
      fit <- obs$reduced.X[[i]]%*%b
      list(fitted.m = fit, res.m = X$dfs[[i]]$Y - fit)
    }
  }
  else {
    loopfun1 <- function(i, ...) {
      if(length(X$dfs) == 1) {
        df <- X$dfs[[1]]
        df$Y <- X$Y[,i]
      } else df <- X$dfs[[i]]
      mod.red <- lm(formula.reduced, data=df, ...)
      list(fitted.m = mod.red$fitted.values, res.m = mod.red$residuals)
    }
  }
  if(is.null(cl)) mclapply_res <- do.call(mclapply, c(list(X=1:X$nr, FUN=loopfun1, mc.cores=mc.cores), mc.args, ...))
  else mclapply_res <- parLapply(cl, 1:X$nr, loopfun1, ...)
  fitted.m <- sapply(mclapply_res, function(x) x$fitted.m)
  res.m <- sapply(mclapply_res, function(x) x$res.m)
  # Simulations by permuting the residuals + F-values for each permutation
  Yperm <- function() { # Permutation
    permutation <- sample(1:nrow(res.m), size=nrow(res.m), replace=FALSE)
    # Permute the residuals (rows in res.m) and create new 'y'
    fitted.m + res.m[permutation, ]
  }
  # Regress the permuted data against the full model and get a new effect of interest
  if(fast) {
    loopfun2 <- function(i, ...) {
      genFvaluesSim(Yperm(), obs$full.X, obs$reduced.X)
    }
  }
  else {
    loopfun2 <- function(i, ...) {
      genFvaluesLM(Yperm(), X$dfs, formula.full, formula.reduced, ...)
    }
  }
  if(is.null(cl)) sim <- do.call(mclapply, c(list(X=1:nsim, FUN=loopfun2, mc.cores=mc.cores), mc.args, ...))
  else sim <- parLapply(cl, 1:nsim, loopfun2)
  sim <- simplify2array(sim)
  if(fast) obs <- obs$Fvalues

  cset <- create_curve_set(list(r = X$r, obs = obs, sim_m = sim))
  res <- do.call(global_envelope_test, c(list(curve_sets=cset, alternative="greater"), GET.args))
  attr(res, "ylab") <- "F(r)"
  attr(res, "yexp") <- quote(F(r))
  attr(res, "fname") <- "F"
  if(savefuns) attr(res, "simfuns") <- cset
  res
}

#' Graphical functional GLM for images
#'
#' Non-parametric graphical tests of significance in functional general linear model (GLM)
#' for images
#'
#' @inheritParams graph.flm
#' @param image_sets A named list of sets of images giving the dependent variable (Y), and
#' possibly additionally all the factors. The dimensions of the elements should
#' match with each other, i.e. the factor values should be given for each argument value
#' and each function.
#' @param ... Additional parameters to be passed to \code{\link{graph.flm}}.
#' The possibly saved simulations are currently only provided in a vector format.
#' @seealso \code{\link{graph.flm}}, \code{\link{frank.flm2d}}
#' @return A \code{global_envelope2d} or \code{combined_global_envelope2d} object,
#' which can be printed and plotted directly.
#' @export
#' @references
#' Mrkvička, T., Roskovec, T. and Rost, M. (2019) A nonparametric graphical tests of significance in functional GLM. arXiv:1902.04926 [stat.ME]
#'
#' Freedman, D., & Lane, D. (1983) A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics, 1(4), 292-298. doi:10.2307/1391660
#' @examples
#' \donttest{
#' data("imageset2")
#' # Testing discrete factor group
#' res.g <- graph.flm2d(nsim = 19, # Increase nsim for serious analysis!
#'                      formula.full = Y ~ group + z,
#'                      formula.reduced = Y ~ z,
#'                      image_sets = list(Y = imageset2$image_set),
#'                      factors = data.frame(group = imageset2$Group,
#'                                           z = imageset2$z))
#' plot(res.g)
#' # Testing discrete factor group with contrasts
#' res.gc <- graph.flm2d(nsim = 19, # Increase nsim for serious analysis!
#'                       formula.full = Y ~ group + z,
#'                       formula.reduced = Y ~ z,
#'                       image_sets = list(Y = imageset2$image_set),
#'                       factors = data.frame(group = imageset2$Group,
#'                                            z = imageset2$z),
#'                       contrasts = TRUE)
#' plot(res.gc)
#'
#' # Testing continuous factor z
#' res.z <- graph.flm2d(nsim = 19, # Increase nsim for serious analysis!
#'                      formula.full = Y ~ group + z,
#'                      formula.reduced = Y ~ group,
#'                      image_sets = list(Y = imageset2$image_set),
#'                      factors = data.frame(group = imageset2$Group,
#'                                           z = imageset2$z))
#' plot(res.z)
#' }
graph.flm2d <- function(nsim, formula.full, formula.reduced, image_sets, factors = NULL, ...) {
  if(class(image_sets)[1] == "image_set") image_sets <- list(image_sets)
  obs_d <- lapply(image_sets, function(x) { dim(x$obs) })
  sim_d <- lapply(image_sets, function(x) { dim(x$sim_m) })
  # Check that dimensions of different image sets are the same
  if(!all(sapply(obs_d, FUN=identical, y=obs_d[[1]]))) stop("Dimensions of image sets (obs) do not match.\n")
  if(!all(sapply(sim_d, FUN=identical, y=sim_d[[1]]))) stop("Dimensions of image sets (sim_m) do not match.\n")
  # Check dimensions of each image set
  image_sets <- lapply(image_sets, check_image_set_dimensions)
  # Check equalities of the r values
  if(!all(sapply(image_sets, FUN = function(x) { identical(x$r, y=image_sets[[1]]$r) }))) stop("The r values of image sets do not match.\n")
  # Create curve sets transforming the 2d functions (matrices) to vectors
  curve_sets_v <- lapply(image_sets, image_set_to_curve_set)

  res_v <- graph.flm(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
                      curve_sets=curve_sets_v, factors=factors, ...)
  # Transform the results to 2d
  res <- curve_set_results_to_image_results(res_v, image_sets)
  attr(res, "call") <- match.call()
  if(!is.null(attr(res_v, "simfuns"))) attr(res, "simfuns") <- attr(res_v, "simfuns")
  res
}

#' F rank functional GLM for images
#'
#' Multiple testing in permutation inference for the general linear model (GLM)
#'
#' @inheritParams frank.flm
#' @inheritParams graph.flm2d
#' @param ... Additional parameters to be passed to \code{\link{frank.flm}}.
#' The possibly saved simulations are currently only provided in a vector format.
#' @seealso \code{\link{frank.flm}}, \code{\link{graph.flm2d}}
#' @export
#' @references
#' Mrkvička, T., Myllymäki, M. and Narisetty, N. N. (2019) New methods for multiple testing in permutation inference for the general linear model. arXiv:1906.09004 [stat.ME]
#' @return A \code{global_envelope2d} object, which can be printed and plotted
#' directly.
#' @examples
#' \donttest{
#' data("imageset2")
#' # Testing discrete factor group
#' res.g <- frank.flm2d(nsim = 19, # Increase nsim for serious analysis!
#'                      formula.full = Y ~ group + z,
#'                      formula.reduced = Y ~ z,
#'                      image_sets = list(Y = imageset2$image_set),
#'                      factors = data.frame(group = imageset2$Group,
#'                                           z = imageset2$z))
#' plot(res.g)
#'
#' # Testing continuous factor z
#' res.z <- frank.flm2d(nsim = 19, # Increase nsim for serious analysis!
#'                      formula.full = Y ~ group + z,
#'                      formula.reduced = Y ~ group,
#'                      image_sets = list(Y = imageset2$image_set),
#'                      factors = data.frame(group = imageset2$Group,
#'                                           z = imageset2$z))
#' plot(res.z)
#' }
frank.flm2d <- function(nsim, formula.full, formula.reduced, image_sets, factors = NULL, ...) {
  if(class(image_sets)[1] == "image_set") image_sets <- list(image_sets)
  obs_d <- lapply(image_sets, function(x) { dim(x$obs) })
  sim_d <- lapply(image_sets, function(x) { dim(x$sim_m) })
  # Check that dimensions of different image sets are the same
  if(!all(sapply(obs_d, FUN=identical, y=obs_d[[1]]))) stop("Dimensions of image sets (obs) do not match.\n")
  if(!all(sapply(sim_d, FUN=identical, y=sim_d[[1]]))) stop("Dimensions of image sets (sim_m) do not match.\n")
  # Check dimensions of each image set
  image_sets <- lapply(image_sets, check_image_set_dimensions)
  # Check equalities of the r values
  if(!all(sapply(image_sets, FUN = function(x) {identical(x$r, y=image_sets[[1]]$r)}))) stop("The r values of image sets do not match.\n")
  # Create curve sets transforming the 2d functions (matrices) to vectors
  curve_sets_v <- lapply(image_sets, image_set_to_curve_set)

  res_v <- frank.flm(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
                      curve_sets=curve_sets_v, factors=factors, ...)
  # Transform the results to 2d
  res <- curve_set_results_to_image_results(res_v, image_sets)
  attr(res, "call") <- match.call()
  if(!is.null(attr(res_v, "simfuns"))) attr(res, "simfuns") <- attr(res_v, "simfuns")
  res
}
