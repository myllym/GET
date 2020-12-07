# Preliminary checks for the graph.flm and frank.flm
#' @importFrom stats terms
flm.checks <- function(nsim, formula.full, formula.reduced, curve_sets, factors = NULL, fast = TRUE) {
  # Preliminary checks
  vars <- all.vars(formula.full)
  vars.reduced <- all.vars(formula.reduced)
  # Check that formula.reduced is nested within formula.full
  check_isnested(formula.full, formula.reduced)
  if(nsim < 1) stop("Not a reasonable value of nsim.")
  if(vars[1] != "Y") stop("The formula should be off the form Y ~ .... where Y is the response.")
  if(!(length(class(curve_sets)) == 1 && class(curve_sets) == "list")) {
    curve_sets <- list(Y=curve_sets)
  }
  available <- unique(c(names(curve_sets), names(factors)))
  if(!all(vars %in% available)) stop("The variables in the formula not found in the given data (curve_sets and factors).")
  if(!all(sapply(curve_sets, function(x) inherits(x, c("curve_set", "fdata"))))) stop("The components of curve_sets do not have a valid class.")
  curve_sets <- lapply(curve_sets, convert_fdata)
  if(any(sapply(curve_sets, function(x) curve_set_is1obs(x)))) stop("All (data) functions of the curve_set must be equal.")
  curve_sets <- check_curve_set_dimensions(curve_sets)
  # Put Y and factors into data.l
  data.l <- list()
  data.l[['Y']] <- data_and_sim_curves(curve_sets[['Y']]) # -> each row corresponds to a data function
  # The argument values
  r <- curve_sets[['Y']][['r']]
  Nfunc <- nrow(data.l[['Y']]) # Number of functions
  nr <- ncol(data.l[['Y']]) # Number of argument values
  vars.csets <- vars[vars %in% names(curve_sets)]
  factors_in_curvesets <- !fast
  if(length(curve_sets) > 1 & length(vars.csets) > 1) { # Factors provided in the curve_sets
    for(i in 2:length(vars.csets)) data.l[[vars.csets[i]]] <- data_and_sim_curves(curve_sets[[vars.csets[i]]])
    factors_in_curvesets <- TRUE
  }
  vars.factors <- vars[vars %in% names(factors)]
  if(!is.null(factors) & length(vars.factors) > 0) {
    if(!inherits(factors, "data.frame")) stop("Invalid factors argument.")
    if(nrow(factors) != Nfunc) stop("The dimensions of Y and factors do not match.")
    # Expand the factors to each argument value
    for(i in 1:length(vars.factors)) {
      data.l[[vars.factors[i]]] <- if(factors_in_curvesets) {
        A <- rep(factors[,vars.factors[i]], times=nr)
        dim(A) <- c(nrow(factors), nr)
        A
      } else { factors[,vars.factors[i]] }
    }
  }
  dfs <- if(factors_in_curvesets) {
    lapply(1:nr, FUN = function(i) { as.data.frame(lapply(data.l, FUN = function(x) x[,i])) })
  } else {
    list(as.data.frame(data.l[-1]))
  }

  list(Y=data.l[['Y']], dfs=dfs, r=r, Nfunc=Nfunc, nr=nr)
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
  if(k == 1) stop("The option \'contrasts\' only valid for discrete factors with at least two levels.")
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

# Compute large multiple linear models in parts to save memory
#' @importFrom stats lm deviance df.residual
gendevianceMLM <- function(Y, df, formula, partsize=200, ...) {
  nr <- ncol(Y)
  nparts <- ceiling(nr/partsize)
  i0 <- round(seq(1, nr+1, length=nparts+1))
  dev <- numeric(nr)
  for(i in 1:nparts) {
    ii <- i0[i]:(i0[i+1]-1)
    df$Y <- Y[,ii]
    M.full <- lm(formula = formula, data = df, model = FALSE, ...)
    dev[ii] <- deviance(M.full)
    dof <- df.residual(M.full)
  }
  list(dev=dev, df=dof)
}

# If covariates are constant with respect to location,
# use a multiple linear model (instead of multiple linear models)
genFvaluesMLM <- function(Y, df, formula.full, formula.reduced, ...) {
  full <- gendevianceMLM(Y, df, formula.full, ...)
  red <- gendevianceMLM(Y, df, formula.reduced, ...)
  # F-statistic
  (full$dev - red$dev)/(full$df - red$df)/(full$dev/full$df)
}

# Compute F-statistics for several pixels at a time in the simplest case of no extra arguments
# and no curve sets in factors.
#' @importFrom stats lm df.residual
genFvaluesSimple <- function(Y, df, formula.full, formula.reduced, partsize=200, ...) {
  df$Y <- Y[,1]
  fx <- function(fit) {
    X <- fit$x
    list(X1 = X %*% solve(t(X)%*%X), X2 = t(X), dof=df.residual(fit))
  }
  x.full <- fx(lm(formula.full, df, model=FALSE, x=TRUE))
  x.red <- fx(lm(formula.reduced, df, model=FALSE, x=TRUE))

  nr <- ncol(Y)
  Fstat <- numeric(nr)

  nparts <- ceiling(nr/partsize)
  i0 <- round(seq(1, nr+1, length=nparts+1))
  for(i in 1:nparts) {
    ii <- i0[i]:(i0[i+1]-1)
    y <- Y[,ii]
    rssf <- colSums((y - x.full[[1]] %*% (x.full[[2]] %*% y))^2)
    rssr <- colSums((y - x.red[[1]] %*% (x.red[[2]] %*% y))^2)

    Fstat[ii] <- (rssf - rssr)/(x.full[[3]] - x.red[[3]])/(rssf/x.full[[3]])
  }
  Fstat
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
#' The function \code{graph.flm} performs the graphical functional GLM of Mrkvička et al. (2019),
#' described also in Section 3.6 of Myllymäki and Mrkvička (2020) (type \code{vignette("GET")} in R).
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
#' nested within the full model. The full model should include in addition to the reduced
#' model the interesting factors whose effects are under investigation.
#' The implementation to find the coefficients of the interesting factors is based on
#' \code{\link{dummy.coef}} and the restrictions there apply.
#'
#' The regression coefficients serve as test functions in the graphical functional GLM.
#' For a continuous interesting factor, the test function is its regression coefficient across the
#' functional domain. For a discrete factor, there are two possibilities that are controlled by
#' the arguments \code{contrasts}. If \code{contrasts = FALSE}, then the test statistic is
#' the function/long vector where the coefficients related to all levels of the factor are joined
#' together. If \code{contrasts = TRUE}, then the differences between the same coefficients are
#' considered instead. Given the coefficients in a specific order that is obtained through the use
#' of \code{\link{lm}} and \code{\link{dummy.coef}}, the differences are taken for couples i and j
#' where i < j and reducing j from i (e.g. for three groups 1,2,3, the constrasts are 1-2, 1-3, 2-3).
#'
#' There are different versions of the implementation depending on the application.
#' Given that the argument \code{fast} is TRUE, then
#' \itemize{
#' \item If all the covariates are constant across the functions, i.e. they can be provided in the
#' argument \code{factors}, then a linear model is fitted separately by least-squares estimation to
#' the data at each argument value of the functions fitting a multiple linear model by \code{\link[stats]{lm}}.
#' The possible extra arguments passed in \code{lm.args} to \code{\link[stats]{lm}} must be of the form that
#' \code{\link[stats]{lm}} accepts for fitting a multiple linear model. In the basic case, no extra arguments are
#' needed.
#' \item If some of the covariates vary across the space and there are user specified extra arguments given in
#' \code{lm.args}, then the implementation fits a linear model at each argument value of the functions using
#' \code{\link[stats]{lm}}, which can be rather slow. The arguments \code{lm.args} are passed to \code{\link[stats]{lm}}
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
#' possibly additionally factors whose values vary across the argument values of the functions.
#' The dimensions of the elements should match with each other.
#' Note that factors that are fixed across the functions can be given in the argument \code{factors}.
#' Also \code{\link[fda.usc]{fdata}} objects allowed.
#' @param factors A data frame of factors. An alternative way to specify factors when they
#' are constant for all argument values of the functions. The number of rows of the data frame should be equal
#' to the number of curves. Each column should specify the values of a factor.
#' @param lm.args A named list of additional arguments to be passed to \code{\link[stats]{lm}}. See details.
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
#' Mrkvička, T., Roskovec, T. and Rost, M. (2019) A nonparametric graphical tests of significance in functional GLM. Methodology and Computing in Applied Probability. doi: 10.1007/s11009-019-09756-y
#'
#' Myllymäki, M and Mrkvička, T. (2020). GET: Global envelopes in R. arXiv:1911.06583 [stat.ME]
#'
#' Freedman, D., & Lane, D. (1983) A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics, 1(4), 292-298. doi:10.2307/1391660
#' @importFrom stats lm
#' @importFrom stats predict.lm
#' @importFrom stats dummy.coef
#' @importFrom parallel mclapply
#' @importFrom parallel parLapply
#' @importFrom parallel clusterEvalQ
#' @examples
#' data("rimov")
#' res <- graph.flm(nsim=19, # Increase the number of simulations for serious analysis!
#'                  formula.full = Y~Year,
#'                  formula.reduced = Y~1,
#'                  curve_sets = list(Y=rimov), factors = data.frame(Year = 1979:2014))
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
#' # An example of testing the joint effect of a discrete and a continuous variable
#' \dontshow{nsim <- 19}
#' \donttest{nsim <- 999}
#' data("GDPtax")
#' factors.df <- data.frame(Group = GDPtax$Group, Tax = GDPtax$Profittax)
#' res.tax_within_group <- graph.flm(nsim = nsim,
#'                                   formula.full = Y~Group+Tax+Group:Tax,
#'                                   formula.reduced = Y~Group+Tax,
#'                                   curve_sets = list(Y=GDPtax$GDP),
#'                                   factors = factors.df)
#' plot(res.tax_within_group)
#'
#' # Image data examples
#'
#' data("abide_9002_23")
#' iset <- abide_9002_23$curve_set
#' \dontshow{
#' # Cut the data to reduce time
#' iset$r <- iset$r[1:29,]
#' iset$funcs <- iset$funcs[1:29, ]
#' }
#'
#' # Figure of an image in the group 1 and group 2
#' plot(iset, idx=c(1, 27))
#'
#' # Testing the discrete factor 'group' with contrasts
#' # (Use contrasts = FALSE for 'means')
#' res <- graph.flm(nsim = 19, # Increase nsim for serious analysis!
#'                  formula.full = Y ~ Group + Sex + Age,
#'                  formula.reduced = Y ~ Sex + Age,
#'                  curve_sets = list(Y = iset),
#'                  factors = abide_9002_23[['factors']],
#'                  contrasts = TRUE,
#'                  GET.args = list(type = "area"))
#' plot(res)
#'
#' # Testing the continuous factor 'age'
#' res.a <- graph.flm(nsim = 19, # Increase nsim for serious analysis!
#'                  formula.full = Y ~ Group + Sex + Age,
#'                  formula.reduced = Y ~ Group + Sex,
#'                  curve_sets = list(Y = iset),
#'                  factors = abide_9002_23[['factors']],
#'                  GET.args = list(type = "area"))
#' plot(res.a)
graph.flm <- function(nsim, formula.full, formula.reduced, curve_sets, factors = NULL,
                      contrasts = FALSE, savefuns = FALSE, lm.args = NULL, GET.args = NULL,
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

  nameinteresting <- labels_interesting(formula.full, formula.reduced)

  # setting that 'fun' is a function
  if(!contrasts) {
    genCoef <- genCoefmeans.m
    ylab <- expression(italic(hat(beta)[i](r)))
  }
  else {
    genCoef <- genCoefcontrasts.m
    ylab <- expression(italic(hat(beta)[i](r)-hat(beta)[j](r)))
  }

  # Fit the full model to the data and obtain the coefficients
  obs <- do.call(genCoef, c(list(X$Y, X$dfs, formula.full, nameinteresting), lm.args))

  #-- Freedman-Lane procedure
  # Fit the reduced model at each argument value to get the fitted values and residuals
  if(length(X$dfs) == 1) {
    df <- X$dfs[[1]]
    df$Y <- X$Y
    fit <- do.call(lm, c(list(formula.reduced, data=df, model=FALSE), lm.args))
    fitted.m <- fit$fitted.values
    res.m <- fit$residuals
    fit <- NULL
  } else {
    loopfun1 <- function(i) {
      mod.red <- do.call(lm, c(list(formula.reduced, data=X$dfs[[i]]), lm.args))
      list(fitted.m = mod.red$fitted.values, res.m = mod.red$residuals)
    }
    if(is.null(cl)) mclapply_res <- do.call(mclapply, c(list(X=1:X$nr, FUN=loopfun1, mc.cores=mc.cores), mc.args))
    else mclapply_res <- parLapply(cl, 1:X$nr, loopfun1)
    fitted.m <- sapply(mclapply_res, function(x) x$fitted.m)
    res.m <- sapply(mclapply_res, function(x) x$res.m)
  }
  # Simulations by permuting the residuals + calculate the coefficients
  Yperm <- function() { # Permutation
    permutation <- sample(1:nrow(res.m), size=nrow(res.m), replace=FALSE)
    # Permute the residuals (rows in res.m) and create new 'y'
    fitted.m + res.m[permutation, ]
  }
  loopfun2 <- function(i) {
    # Regress the permuted data against the full model and get a new effect of interest
    do.call(genCoef, c(list(Yperm(), X$dfs, formula.full, nameinteresting), lm.args))
  }
  if(is.null(cl)) sim <- do.call(mclapply, c(list(X=1:nsim, FUN=loopfun2, mc.cores=mc.cores), mc.args))
  else sim <- parLapply(cl, 1:nsim, loopfun2)
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
  # Re-define the default ylab
  res <- envelope_set_labs(res, ylab=ylab)
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
#' consequently, the test achieves the desired significance level only approximately.
#' If the reduced model contains only a constant, then the algorithm corresponds to
#' simple permutation of raw data.
#' In contrast to the graphical functional GLM, the F rank functional GLM is based on the
#' F-statistics that are calculated at each argument value of the functions.
#' The global envelope test is applied to the observed and simulated F-statistics.
#' The test is able to find if the factor of interest is significant and also which
#' argument values of the functional domain are responsible for the potential rejection.
#'
#' The specification of the full and reduced formulas is important. The reduced model should be
#' nested within the full model. The full model should include in addition to the reduced
#' model the interesting factors whose effects are under investigation.
#'
#' There are different versions of the implementation depending on the application.
#' \itemize{
#' \item If all the covariates are constant across the functions, i.e. they can be provided in the
#' argument \code{factors}, and there are no extra arguments given by the user in \code{lm.args}, then a
#' fast implementation is used to directly compute the F-statistics.
#' \item If all the covariates are constant across the functions, but there are some extra arguments,
#' then a linear model is fitted separately by least-squares estimation to
#' the data at each argument value of the functions fitting a multiple linear model by \code{\link[stats]{lm}}.
#' The possible extra arguments passed in \code{lm.args} to \code{\link[stats]{lm}} must be of the form that
#' \code{\link[stats]{lm}} accepts for fitting a multiple linear model. In the basic case, no extra arguments are
#' needed.
#' \item If some of the covariates vary across the space, i.e. they are provided in the list of curve sets in
#' the argument \code{curve_sets} together with the dependent functions, but there are no extra arguments given
#' by the user in \code{lm.args}, there is a rather fast implementation of the F-value calculation (which does not
#' use \code{\link[stats]{lm}}).
#' \item If some of the covariates vary across the space and there are user specified extra arguments given in
#' \code{lm.args}, then the implementation fits a linear model at each argument value of the functions using
#' \code{\link[stats]{lm}}, which can be rather slow. The arguments \code{lm.args} are passed to \code{\link[stats]{lm}}
#' for fitting each linear model.
#' }
#' By default the fastest applicable method is used. This can be changed by setting \code{method} argument.
#' The cases above correspond to "simple", "mlm", "complex" and "lm". Changing the default can be useful for
#' checking the validity of the implementation.
#'
#' @inheritParams graph.flm
#' @param savefuns Logical or "return". If TRUE, then the functions from permutations are saved to the attribute simfuns.
#' If "return", then the function returns the permutations in a curve_set, instead of the result of the envelope test on those;
#' this can be used by \code{\link{partial_forder}}.
#' @param method For advanced use.
#' @return A \code{global_envelope} object, which can be printed and plotted directly.
#' @export
#' @references
#' Mrkvička, T., Myllymäki, M., Kuronen, M. and Narisetty, N. N. (2020) New methods for multiple testing in permutation inference for the general linear model. arXiv:1906.09004 [stat.ME]
#'
#' Freedman, D., & Lane, D. (1983) A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics, 1(4), 292-298. doi:10.2307/1391660
#' @examples
#' data("GDPtax")
#' factors.df <- data.frame(Group = GDPtax$Group, Tax = GDPtax$Profittax)
#' \dontshow{nsim <- 19}
#' \donttest{nsim <- 999}
#' res.tax_within_group <- frank.flm(nsim = nsim,
#'                                   formula.full = Y~Group+Tax+Group:Tax,
#'                                   formula.reduced = Y~Group+Tax,
#'                                   curve_sets = list(Y=GDPtax$GDP),
#'                                   factors = factors.df)
#' plot(res.tax_within_group)
#'
#' # Image set examples
#' data("abide_9002_23")
#' iset <- abide_9002_23$curve_set
#' \dontshow{
#' # Cut the data to reduce time
#' iset$r <- iset$r[1:29,]
#' iset$funcs <- iset$funcs[1:29, ]
#' }
#' res.F <- frank.flm(nsim = 19, # Increase nsim for serious analysis!
#'                    formula.full = Y ~ Group + Age + Sex,
#'                    formula.reduced = Y ~ Age + Sex,
#'                    curve_sets = list(Y = iset),
#'                    factors = abide_9002_23[['factors']],
#'                    GET.args = list(type = "area"))
#' plot(res.F)
# Freedman-Lane procedure (Freedman and Lane, 1983, p. 385)
#' @importFrom parallel mclapply
#' @importFrom parallel parLapply
#' @importFrom stats lm
frank.flm <- function(nsim, formula.full, formula.reduced, curve_sets, factors = NULL,
                      savefuns = TRUE, lm.args = NULL, GET.args = NULL,
                      mc.cores = 1, mc.args = NULL, cl = NULL,
                      method = c("best", "simple", "mlm", "complex", "lm")) {
  # Preliminary checks and formulation of the data to suitable form for further processing
  method <- match.arg(method)
  X <- flm.checks(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
                   curve_sets=curve_sets, factors=factors, fast=method %in% c("mlm", "simple", "best"))

  if(method == "best") {
    method <- if(length(lm.args) > 0) {
      if(length(X$dfs) > 1) "lm"
      else "mlm"
    } else {
      if(length(X$dfs) > 1) "complex"
      else "simple"
    }
  }
  if(length(lm.args) > 0) {
    # This is only a warning because the user might know that the extra args in lm.args are ok.
    if(method %in% c("simple", "complex")) warning("Method ", method, " doesn't work with extra arguments (lm.args).")
  }
  if(length(X$dfs) > 1 && method %in% c("mlm", "simple")) {
    stop("Curvesets in factors not allowed with method='", method, "'")
  }
  # Calculate the F-statistic for the data, and, if method="complex", obtain also the design matrices
  obs <- switch(method,
                complex=genFvaluesObs(X$dfs, formula.full, formula.reduced),
                simple=genFvaluesSimple(X$Y, X$dfs[[1]], formula.full, formula.reduced),
                lm=,
                mlm=do.call(genFvaluesLM, c(list(X$Y, X$dfs, formula.full, formula.reduced), lm.args)))
  # Freedman-Lane procedure
  # Fit the reduced model at each argument value to get fitted values and residuals
  if(method %in% c("mlm", "simple")) {
    df <- X$dfs[[1]]
    df$Y <- X$Y
    fit <- do.call(lm, c(list(formula.reduced, data=df, model=FALSE), lm.args))
    fitted.m <- fit$fitted.values
    res.m <- fit$residuals
    fit <- NULL
  } else {
    if(method == "complex") {
      loopfun1 <- function(i) {
        b <- bcoef(Y = X$dfs[[i]]$Y, X = obs$reduced.X[[i]])
        fit <- obs$reduced.X[[i]]%*%b
        list(fitted.m = fit, res.m = X$dfs[[i]]$Y - fit)
      }
    } else { # method == "lm"
      loopfun1 <- function(i) {
        mod.red <- do.call(lm, c(list(formula.reduced, data=X$dfs[[i]]), lm.args))
        list(fitted.m = mod.red$fitted.values, res.m = mod.red$residuals)
      }
    }
    if(is.null(cl)) mclapply_res <- do.call(mclapply, c(list(X=1:X$nr, FUN=loopfun1, mc.cores=mc.cores), mc.args))
    else mclapply_res <- parLapply(cl, 1:X$nr, loopfun1)
    fitted.m <- sapply(mclapply_res, function(x) x$fitted.m)
    res.m <- sapply(mclapply_res, function(x) x$res.m)
  }
  # Original data is not needed anymore
  X$Y <- NULL
  # Simulations by permuting the residuals + F-values for each permutation
  Yperm <- function() { # Permutation
    permutation <- sample(1:nrow(res.m), size=nrow(res.m), replace=FALSE)
    # Permute the residuals (rows in res.m) and create new 'y'
    fitted.m + res.m[permutation, ]
  }
  # Regress the permuted data against the full model and get a new effect of interest
  loopfun2 <- switch(method,
                     lm=,
                     mlm=function(i) do.call(genFvaluesLM, c(list(Yperm(), X$dfs, formula.full, formula.reduced), lm.args)),
                     complex=function(i) genFvaluesSim(Yperm(), obs$full.X, obs$reduced.X),
                     simple=function(i) genFvaluesSimple(Yperm(), X$dfs[[1]], formula.full, formula.reduced))

  if(is.null(cl)) sim <- do.call(mclapply, c(list(X=1:nsim, FUN=loopfun2, mc.cores=mc.cores), mc.args))
  else sim <- parLapply(cl, 1:nsim, loopfun2)
  sim <- simplify2array(sim)
  if(method == "complex") obs <- obs$Fvalues

  cset <- create_curve_set(list(r = X$r, obs = obs, sim_m = sim))
  if(savefuns=="return") return(cset)
  res <- do.call(global_envelope_test, c(list(curve_sets=cset, alternative="greater"), GET.args))
  res <- envelope_set_labs(res, ylab = expression(italic(F(r))))
  if(savefuns) attr(res, "simfuns") <- cset
  res
}
