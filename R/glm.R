# Preliminary checks for the graph.fglm and frank.fglm
#' @importFrom stats terms
fglm.checks <- function(nsim, formula.full, formula.reduced, curve_sets, factors = NULL) {
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
  if(!all(lapply(curve_sets[['obs']], is.matrix))) stop("The curve_set must include data functions (sim_m ignored).\n")
  curve_sets <- check_curve_set_dimensions(curve_sets)
  # Put the functions and factors into data.l
  data.l <- list()
  data.l[['Y']] <- t(curve_sets[[vars[1]]][['obs']]) # -> each row corresponds to a data function
  # The argument values
  r <- curve_sets[[vars[1]]][['r']]
  Nfunc <- nrow(data.l[[1]]) # Number of functions
  nr <- ncol(data.l[[1]]) # Number of argument values
  vars.csets <- vars[vars %in% names(curve_sets)]
  if(length(curve_sets) > 1 & length(vars.csets) > 1) { # Factors provided in the curve_sets
    for(i in 2:length(vars.csets)) data.l[[vars.csets[i]]] <- t(curve_sets[[vars.csets[i]]][['obs']])
  }
  vars.factors <- vars[vars %in% names(factors)]
  if(!is.null(factors) & length(vars.factors) > 0) {
    if(class(factors)[1] != "data.frame") stop("Invalid factors argument.\n")
    if(nrow(factors) != nrow(data.l[[1]])) stop("The dimensions of Y and factors do not match.\n")
    # Expand the factors to each argument value
    for(i in 1:length(vars.factors)) data.l[[vars.factors[i]]] <- matrix(factors[,vars.factors[i]], nrow=nrow(factors), ncol=nr)
  }
  list(data.l=data.l, r=r, Nfunc=Nfunc, nr=nr, einfo=einfo)
}

# M1 = full model; M2 = reduced model
# Return the names of the coefficients that are in M1, but not in M2
factorname_diff <- function(M1, M2) {
  # Names of the coefficients in the full model
  namecoef.full <- names(unlist(stats::dummy.coef(M1)))
  # Names of the coefficients in the reduced model
  namecoef.red <- names(unlist(stats::dummy.coef(M2)))
  # Change the names of the coefficients in the reduced model,
  # if the full model includes discrete factors, but the reduced model not
  if(length(M2$xlevels) == 0 & length(M1$xlevels) > 0) {
    namecoef.red <- paste(namecoef.red, ".", namecoef.red, sep="")
  }
  # The interesting coefficients
  setdiff(namecoef.full, namecoef.red)
}

# Regress the given data (true or permuted) against the full model and
# get an effect of interest at all r values in a matrix.
# data.l = a list containing data (Y and factors), all variables in formula.full
# nameinteresting = names of the interesting coefficients
#' @importFrom stats lm
#' @importFrom stats dummy.coef
genCoefmeans.m <- function(data.l, formula.full, nameinteresting, ...) {
  nr <- ncol(data.l[[1]])
  effect.a <- matrix(0, nrow=nr, ncol=length(nameinteresting))
  dimnames(effect.a) <- list(NULL, nameinteresting)
  for(i in 1:nr) {
    df <- as.data.frame(lapply(data.l, FUN = function(x) x[,i]))
    M.full <- stats::lm(formula.full, data=df, ...)
    allcoef <- unlist(stats::dummy.coef(M.full))
    effect.a[i,] <- allcoef[nameinteresting]
  }
  effect.a
}

# The constrasts of the interesting effects in a matrix
#' @importFrom stats lm
#' @importFrom stats dummy.coef
genCoefcontrasts.m <- function(data.l, formula.full, nameinteresting, ...) {
  nr <- ncol(data.l[[1]])
  k <- length(nameinteresting)
  effect.a <- genCoefmeans.m(data.l=data.l, formula.full=formula.full, nameinteresting=nameinteresting, ...)
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
#' @importFrom stats lm
#' @importFrom stats anova
genFvaluesLM <- function(data.l, formula.full, formula.reduced, ...) {
  nr <- ncol(data.l[[1]])
  Fvalues <- vector(length=nr)
  for(i in 1:nr) {
    df <- as.data.frame(lapply(data.l, FUN = function(x) x[,i]))
    M.full <- stats::lm(formula = formula.full, data = df, ...)
    M.reduced <- stats::lm(formula = formula.reduced, data = df, ...)
    Anova.res <- stats::anova(M.reduced, M.full)
    Fvalues[i] <- Anova.res$F[2]
  }
  Fvalues
}

# Parameter estimate b for lm
bcoef <- function(Y, X) {
  solve(a = t(X)%*%X, b = t(X) %*% Y)
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
genFvaluesObs <- function(data.l, formula.full, formula.reduced) {
  nr <- ncol(data.l[[1]])
  Fvalues <- vector(length=nr)
  # Case i = 1
  df <- as.data.frame(lapply(data.l, FUN = function(x) x[,1]))
  # Call first lm to obtain the design matrices
  M.full <- stats::lm(formula = formula.full, data = df, x = TRUE)
  M.reduced <- stats::lm(formula = formula.reduced, data = df, x = TRUE)
  # Design matrices are M.full$x, M.reduced$x (the same for each i)
  Fvalues[1] <- Fvalue(Y = df$Y, X1 = M.full$x, X2 = M.reduced$x)
  for(i in 2:nr) {
    Fvalues[i] <- Fvalue(Y = data.l$Y[,i], X1 = M.full$x, X2 = M.reduced$x)
  }
  list(Fvalues = Fvalues, full.X = M.full$x, reduced.X = M.reduced$x)
}

# General F-value from lm-model
# Y = data, the dependent
# designX.full = design matrix of the full model
# designX.reduced = design matrix of the reduced model
genFvaluesSim <- function(Y, designX.full, designX.reduced) {
  nr <- ncol(Y)
  Fvalues <- vector(length=nr)
  for(i in 1:nr) {
    Fvalues[i] <- Fvalue(Y = Y[,i], X1 = designX.full, X2 = designX.reduced)
  }
  Fvalues
}

#' Graphical functional GLM
#'
#' Non-parametric graphical tests of significance in functional general linear model (GLM)
#'
#'
#' The function \code{graph.fglm} performs the graphical functional GLM of Mrkvička et al. (2019).
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
#' @inheritParams graph.fanova
#' @inheritParams GET.composite
#' @param formula.full The formula specifying the general linear model,
#' see \code{formula} in \code{\link[stats]{lm}}.
#' @param formula.reduced The formula of the reduced model with nuisance factors only.
#' @param curve_sets A named list of sets of curves giving the dependent variable (Y), and
#' possibly additionally all the factors. The dimensions of the elements should
#' match with each other, i.e. the factor values should be given for each argument value
#' and each function. If factors are given in the argument \code{factors}, then can also be just
#' the curve set representing Y. Also \code{\link[fda.usc]{fdata}} objects allowed.
#' @param factors A data frame of factors. An alternative way to specify factors when they
#' are constant for all argument values. The number of rows of the data frame should be equal
#' to the number of curves. Each column should specify the values of a factor.
#' @param ... Additional arguments to be passed to \code{\link[stats]{lm}}.
#' @param GET.args A named list of additional arguments to be passed to \code{\link{global_envelope_test}}.
#' @param mc.args A named list of additional arguments to be passed to \code{\link{mclapply}}.
#' Only relevant if \code{mc.cores} is more than 1.
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
#' @examples
#' data(rimov)
#' \donttest{
#' res <- graph.fglm(nsim=19, # Increase the number of simulations for serious analysis!
#'                   formula.full = Y~Year,
#'                   formula.reduced = Y~1,
#'                   curve_sets = list(Y=rimov), factors = data.frame(Year = 1979:2014))
#' }
#' \dontshow{
#' res <- graph.fglm(nsim=3,
#'                   formula.full = Y~Year,
#'                   formula.reduced = Y~1,
#'                   curve_sets = list(Y=rimov), factors = data.frame(Year = 1979:2014),
#'                   GET.args = list(alpha=0.25))
#' }
#' plot(res)
#'
#' \donttest{
#' data(GDPtax)
#' factors.df <- data.frame(Group = GDPtax$Group, Tax = GDPtax$Profittax)
#' res.tax_within_group <- graph.fglm(nsim = 999,
#'                                   formula.full = Y~Group+Tax+Group:Tax,
#'                                   formula.reduced = Y~Group+Tax,
#'                                   curve_sets = list(Y=GDPtax$GDP),
#'                                   factors = factors.df)
#' plot(res.tax_within_group)
#' }
graph.fglm <- function(nsim, formula.full, formula.reduced, curve_sets, factors = NULL,
                       summaryfun = c("means", "contrasts"),
                       savefuns = FALSE, ..., GET.args = NULL, mc.cores = 1L, mc.args = NULL) {
  # Set up the contrasts
  op <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(op))
  # Preliminary checks and formulation of the data to suitable form for further processing
  X <- fglm.checks(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
                   curve_sets=curve_sets, factors=factors)
  summaryfun <- match.arg(summaryfun)
  # setting that 'fun' is a function
  switch(summaryfun, 
         means = {genCoef = genCoefmeans.m},
         contrasts = {genCoef = genCoefcontrasts.m}
  )

  # Fit the reduced model at each argument value
  loopfun1 <- function(i, ...) {
    df <- as.data.frame(lapply(X$data.l, FUN = function(x) x[,i])) # create the data.frame at the ith argument value
    mod.red <- lm(formula.reduced, data=df, ...)
    list(fitted.m = mod.red$fitted.values, res.m = mod.red$residuals)
  }
  # Simulations by permuting the residuals + calculate the coefficients
  loopfun2 <- function(i, ...) {
    permutation <- sample(1:nrow(res.m), size=nrow(res.m), replace=FALSE)
    # Permute the residuals (rows in res.m) and create new 'y'
    X$data.l[['Y']] <- fitted.m + res.m[permutation, ]
    # Regress the permuted data against the full model and get a new effect of interest
    genCoef(X$data.l, formula.full, nameinteresting, ...)
  }

  # Fit the full model at the first argument value to get the names of interesting coefficients
  df <- as.data.frame(lapply(X$data.l, FUN = function(x) x[,1]))
  mod.full <- lm(formula = formula.full, data=df, ...)
  mod.red <- lm(formula = formula.reduced, data=df, ...)
  nameinteresting <- factorname_diff(mod.full, mod.red)
  message("The inspected coefficients are: \n", paste(nameinteresting, " "), "\n")
  # Fit the full model to the data and obtain the coefficients
  obs <- genCoef(X$data.l, formula.full, nameinteresting, ...)

  #-- Freedman-Lane procedure
  # Fit the reduced model at each argument value to get the fitted values and residuals
  mclapply_res <- do.call(mclapply, c(list(X=1:X$nr, FUN=loopfun1, mc.cores=mc.cores), mc.args, ...))
  fitted.m <- sapply(mclapply_res, function(x) x$fitted.m)
  res.m <- sapply(mclapply_res, function(x) x$res.m)

  # Simulations by permuting the residuals + coefficients for each permutation
  sim <- do.call(parallel::mclapply, c(list(X=1:nsim, FUN=loopfun2, mc.cores=mc.cores), mc.args, ...))
  sim <- sapply(sim, function(x) x, simplify="array")
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
#' The function \code{frank.fglm} performs
#' a nonparametric test of significance of a covariate in the functional GLM.
#' Similarly as in the graphical functional GLM (\code{\link{graph.fglm}}),
#' the Freedman-Lane algorithm (Freedman and Lane, 1983) is applied to permute the functions
#' (to obtain the simulations under the null hypothesis of "no effects");
#' consequently, the test approximately achieves the desired significance level.
#' In contrast to the graphical functional GLM, the F rank functional GLM is based on the F-statistics,
#' that are calculated at each argument value of the functions.
#' The global envelope test is applied to the observed and simulated F-statistics.
#' The test is able to find if the factor of interest is significant and also which
#' argument values of the functional domain are responsible for the potential rejection.
#'
#' @inheritParams graph.fglm
#' @param fast Logical. If TRUE and no additional parameters are passed to \code{\link[stats]{lm}}
#' in \code{...}, then a faster implementation to calculate the test statistics is used.
#' If FALSE, then \code{\link[stats]{lm}} is utilized for the test statistic calculation (slow).
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
#' res.tax_within_group <- frank.fglm(nsim = 999,
#'                                    formula.full = Y~Group+Tax+Group:Tax,
#'                                    formula.reduced = Y~Group+Tax,
#'                                    curve_sets = list(Y=GDPtax$GDP),
#'                                    factors = factors.df)
#' }
#' \dontshow{
#' res.tax_within_group <- frank.fglm(nsim = 4,
#'                                    formula.full = Y~Group+Tax+Group:Tax,
#'                                    formula.reduced = Y~Group+Tax,
#'                                    curve_sets = list(Y=GDPtax$GDP),
#'                                    factors = factors.df,
#'                                    GET.args = list(alpha=0.20))
#' }
#' plot(res.tax_within_group)
# Freedman-Lane procedure (Freedman and Lane, 1983, p. 385)
frank.fglm <- function(nsim, formula.full, formula.reduced, curve_sets, factors = NULL,
                       savefuns = TRUE, ..., GET.args = NULL, mc.cores = 1, mc.args = NULL,
                       fast = TRUE) {
  # Preliminary checks and formulation of the data to suitable form for further processing
  X <- fglm.checks(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
                   curve_sets=curve_sets, factors=factors)

  extraargs <- list(...)
  if(length(extraargs) < 1 & fast) {
    # Fit the reduced model at each argument value
    loopfun1 <- function(i, ...) { # ... ignored
      df <- as.data.frame(lapply(X$data.l, FUN = function(x) x[,i])) # create the data.frame at the ith argument value
      b <- bcoef(Y = df$Y, X = obs$reduced.X)
      fit <- obs$reduced.X%*%b
      list(fitted.m = fit, res.m = df$Y - fit)
    }
    # Simulations by permuting the residuals + calculate F-values
    loopfun2 <- function(i, ...) { # ... ignored
      permutation <- sample(1:nrow(res.m), size=nrow(res.m), replace=FALSE)
      # Permute the residuals (rows in res.m) and create new 'y'
      Yperm <- fitted.m + res.m[permutation, ]
      # Regress the permuted data against the full model and get a new effect of interest
      genFvaluesSim(Yperm, obs$full.X, obs$reduced.X)
    }
  }
  else {
    fast <- FALSE
    # Fit the reduced model at each argument value
    loopfun1 <- function(i, ...) {
      df <- as.data.frame(lapply(X$data.l, FUN = function(x) x[,i])) # create the data.frame at the ith argument value
      mod.red <- lm(formula.reduced, data=df, ...)
      # Save predictions and residuals
      list(fitted.m = mod.red$fitted.values, res.m = mod.red$residuals)
    }
    # Simulations by permuting the residuals + calculate F-values
    loopfun2 <- function(i, ...) {
      permutation <- sample(1:nrow(res.m), size=nrow(res.m), replace=FALSE)
      # Permute the residuals (rows in res.m) and create new 'y'
      X$data.l[['Y']] <- fitted.m + res.m[permutation, ]
      # Regress the permuted data against the full model and get a new effect of interest
      genFvaluesLM(X$data.l, formula.full, formula.reduced, ...)
    }
  }
  # Calculate the F-statistic for the data, and, if fast, obtain the design matrices (obs$full.X and obs$reduced.X)
  if(fast) obs <- genFvaluesObs(X$data.l, formula.full, formula.reduced)
  else obs <- genFvaluesLM(X$data.l, formula.full, formula.reduced, ...)
  # Freedman-Lane procedure
  # Fit the reduced model at each argument value to get fitted values and residuals
  mclapply_res <- do.call(mclapply, c(list(X=1:X$nr, FUN=loopfun1, mc.cores=mc.cores), mc.args, ...))
  fitted.m <- sapply(mclapply_res, function(x) x$fitted.m)
  res.m <- sapply(mclapply_res, function(x) x$res.m)
  # Simulations by permuting the residuals + F-values for each permutation
  sim <- do.call(parallel::mclapply, c(list(X=1:nsim, FUN=loopfun2, mc.cores=mc.cores), mc.args, ...))
  sim <- sapply(sim, function(x) x, simplify="array")
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
#' @inheritParams graph.fglm
#' @param image_sets A named list of sets of images giving the dependent variable (Y), and
#' possibly additionally all the factors. The dimensions of the elements should
#' match with each other, i.e. the factor values should be given for each argument value
#' and each function.
#' @param ... Additional parameters to be passed to \code{\link{graph.fglm}}.
#' The possibly saved simulations are currently only provided in a vector format.
#' @seealso \code{\link{graph.fglm}}, \code{\link{frank.fglm2d}}
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
#' res.g <- graph.fglm2d(nsim = 19, # Increase nsim for serious analysis!
#'                       formula.full = Y ~ group + z,
#'                       formula.reduced = Y ~ z,
#'                       image_sets = list(Y = imageset2$image_set),
#'                       factors = data.frame(group = imageset2$Group,
#'                                            z = imageset2$z))
#' plot(res.g)
#' # Testing discrete factor group with contrasts
#' res.gc <- graph.fglm2d(nsim = 19, # Increase nsim for serious analysis!
#'                        formula.full = Y ~ group + z,
#'                        formula.reduced = Y ~ z,
#'                        image_sets = list(Y = imageset2$image_set),
#'                        factors = data.frame(group = imageset2$Group,
#'                                             z = imageset2$z),
#'                        summaryfun = "contrasts")
#' plot(res.gc)
#'
#' # Testing continuous factor z
#' res.z <- graph.fglm2d(nsim = 19, # Increase nsim for serious analysis!
#'                       formula.full = Y ~ group + z,
#'                       formula.reduced = Y ~ group,
#'                       image_sets = list(Y = imageset2$image_set),
#'                       factors = data.frame(group = imageset2$Group,
#'                                            z = imageset2$z))
#' plot(res.z)
#' }
graph.fglm2d <- function(nsim, formula.full, formula.reduced, image_sets, factors = NULL, ...) {
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

  res_v <- graph.fglm(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
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
#' @inheritParams frank.fglm
#' @inheritParams graph.fglm2d
#' @param ... Additional parameters to be passed to \code{\link{frank.fglm}}.
#' The possibly saved simulations are currently only provided in a vector format.
#' @seealso \code{\link{frank.fglm}}, \code{\link{graph.fglm2d}}
#' @export
#' @references
#' Mrkvička, T., Myllymäki, M. and Narisetty, N. N. (2019) New methods for multiple testing in permutation inference for the general linear model. arXiv:1906.09004 [stat.ME]
#' @return A \code{global_envelope2d} object, which can be printed and plotted
#' directly.
#' @examples
#' \donttest{
#' data("imageset2")
#' # Testing discrete factor group
#' res.g <- frank.fglm2d(nsim = 19, # Increase nsim for serious analysis!
#'                        formula.full = Y ~ group + z,
#'                        formula.reduced = Y ~ z,
#'                        image_sets = list(Y = imageset2$image_set),
#'                        factors = data.frame(group = imageset2$Group,
#'                                             z = imageset2$z))
#' plot(res.g)
#'
#' # Testing continuous factor z
#' res.z <- frank.fglm2d(nsim = 19, # Increase nsim for serious analysis!
#'                       formula.full = Y ~ group + z,
#'                       formula.reduced = Y ~ group,
#'                       image_sets = list(Y = imageset2$image_set),
#'                       factors = data.frame(group = imageset2$Group,
#'                                            z = imageset2$z))
#' plot(res.z)
#' }
frank.fglm2d <- function(nsim, formula.full, formula.reduced, image_sets, factors = NULL, ...) {
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

  res_v <- frank.fglm(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
                      curve_sets=curve_sets_v, factors=factors, ...)
  # Transform the results to 2d
  res <- curve_set_results_to_image_results(res_v, image_sets)
  attr(res, "call") <- match.call()
  if(!is.null(attr(res_v, "simfuns"))) attr(res, "simfuns") <- attr(res_v, "simfuns")
  res
}
