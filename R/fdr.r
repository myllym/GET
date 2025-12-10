# Returns the last component of a vector
vlast <- function(x) {
  x[length(x)]
}

# Helper functions for FDR envelopes

# Calculate E(R0) and R.obs for the rank envelopes (gamma = 1,2,3,4,...)
# E(R0) is calculated using the theoretical expectation.
fdr_rejections_rank <- function(curve_set, alternative) {
  curve_set <- convert_envelope(curve_set)

  obs_curve <- curve_set_1obs(curve_set)
  sim_curves <- data_or_sim_curves(curve_set) # Exclude the data function
  nsim <- curve_set_nfunc(curve_set) - 1
  nr <- curve_set_narg(curve_set) # m in the fdr paper

  # Count the numbers of rejections R0_j for data and for all simulations j,
  # j=1,...,nsim2 and for all envelopes l=1,2,3,...,floor(nsim/2)
  if(alternative=="two.sided") {
    max_l <- floor(nsim/2)
    mult <- 2
  }
  else {
    max_l <- floor(3/4*nsim) # In the one-sided case, go up to 25% envelope).
    mult <- 1
  }

  # The case of just one set of curves
  # Use the theoretical expectation for ER0:
  all_curves <- data_and_sim_curves(curve_set) # Include the data function

  # Calculate pointwise ranks for each argument value (r)
  calc_pointwiserank <- find_calc_pointwiserank("rank", alternative)
  dataranks <- vector(length=nr)
  for(i in 1:nr) {
    dataranks[i] <- calc_pointwiserank(all_curves[,i])[1]
  }
  R.obs <- sapply(1:max_l, FUN = function(i) sum(dataranks <= i))

  # The expected numbers of rejections R0 for all envelopes l=1,2,3...,max_l
  R0 <- nr * mult * (1:max_l)/nsim

  list(R0=R0, R.obs=R.obs, dataranks=dataranks)
}

# Storey & Tibshirani (2001)
fdr_step13 <- function(R.obs, R0) {
  R0 / pmax(R.obs, 1)
}

# Find the index from rejs (R.obs, R0) for which the fdr is below the given limit
find_ind <- function(R.obs, R0, limit, leftright="right") {
  switch(leftright,
    right = { # The last one that is still within the limit
      max(0, vlast(which(fdr_step13(R.obs, R0) <= limit)))
    },
    left = { # The last one before fdr first goes above the limit
      max(0, (which(fdr_step13(R.obs, R0) > limit))[1]) - 1
    }
  )
}

# Calculate the lower and upper bound of the envelope
fdr_envelope_rank <- function(ind, curve_set) {
  # ind corresponds to the l of the lth minimal and maximal pointwise values
  sim_curves <- data_or_sim_curves(curve_set)
  nr <- curve_set_narg(curve_set)
  nsim <- dim(sim_curves)[1]
  LB <- UB <- array(0, nr)
  for(i in 1:nr){
    Hod <- sort(sim_curves[,i])
    LB[i] <- Hod[ind]
    UB[i] <- Hod[nsim-ind+1]
  }
  list(LB=LB, UB=UB)
}

# Given the index of the 'l'ind'th rank envelope, and the number of data rejections for that ind_Robs (to handle a special case),
# find the gamma between suitable the indices [ind, ind+1].
# Note that anyway only the one for which 'ind' is given is useful and the function using this function
# should take care of that.
# Calculates the rejections for the FDRest choice "ST" used in ATSE and IATSE algorithms. No other options.
# llimit, ulimit = Lower and upper limits between which the functions lie,
# e.g. local correlation must be between -1 and 1. Affects the case ind = 0.
# NULL means -Inf and Inf for the lower and upper limits, respectively.
fdr_rejections_between_two_rank <- function(ind, ind_Robs, limit, curve_set,
                                            alternative, llimit=NULL, ulimit=NULL) {
  failed <- FALSE
  obs_curve <- curve_set_1obs(curve_set)
  nr <- curve_set_narg(curve_set)
  nsim <- curve_set_nfunc(curve_set) - 1
  noutside <- function(curve, lower_bound, upper_bound, alternative) {
    switch(alternative,
      two.sided = { sum(curve < lower_bound | curve > upper_bound) },
      less = { sum(curve < lower_bound) },
      greater = { sum(curve > upper_bound) }
    )
  }
  if(alternative == "two.sided") mult <- 2 else mult <- 1
  if(ind == 0) { # The value of gamma is between 0 and 1
    # Find the envelope of rank 1; THE envelope is wider than this
    e2 <- fdr_envelope_rank(1, curve_set)
    # Define the gamma values
    # Smallest gamma needs to be smaller than nsim * limit / ( nr * mult )
    d <- max(0.0001, min(0.0075, nsim * limit / ( nr * mult )))
    gammas <- seq(0+d/2, 1-d/2, by=d)
    # Calculate rejections for each gamma
    e <- vector("list", length(gammas))
    R.obs_new <- R0_new <- vector(length=length(gammas))
    if(is.null(llimit) & is.null(ulimit)) {
      for(i in seq_along(gammas)) {
        if(alternative != "greater") LB <- e2$LB + log(gammas[i]) * (e2$UB - e2$LB) else LB <- -Inf
        if(alternative != "less") UB <- e2$UB - log(gammas[i]) * (e2$UB - e2$LB) else UB <- Inf
        R.obs_new[i] <- noutside(obs_curve, LB, UB, alternative)
        R0_new[i] <- nr*mult*gammas[i]/nsim
        e[[i]] <- list(LB=LB, UB=UB)
      }
    }
    else { # llimit or ulimit
      for(i in seq_along(gammas)) {
        if(alternative != "greater") {
          if(is.null(llimit)) LB <- e2$LB + log(gammas[i]) * (e2$UB - e2$LB)
          else LB <- llimit + (gammas[i])*(e2$LB-llimit)
        }
        else LB <- -Inf
        if(alternative != "less") {
          if(is.null(ulimit)) UB <- e2$UB - log(gammas[i]) * (e2$UB - e2$LB)
          else UB <- ulimit - (gammas[i])*(ulimit-e2$UB)
        }
        else UB <- Inf
        R.obs_new[i] <- noutside(obs_curve, LB, UB, alternative)
        R0_new[i] <- nr*mult*gammas[i]/nsim
        e[[i]] <- list(LB=LB, UB=UB)
      }
    }
  }
  else { # The value of gamma lies between ind and ind + 1
    # Calculate the envelopes between which THE envelope lies
    sim_curves <- data_or_sim_curves(curve_set)
    e1l <- e1u <- e2l <- e2u <- array(0, nr)
    for(i in 1:nr) {
      Hod <- sort(sim_curves[,i])
      # larger envelope
      e1l[i] <- Hod[ind]
      e1u[i] <- Hod[nsim-ind+1]
      # smaller envelope
      e2l[i] <- Hod[ind+1]
      e2u[i] <- Hod[nsim-ind]
    }
    switch(alternative,
           "two.sided" = {},
           "less" = { e1u <- e2u <- Inf },
           "greater" = { e1l <- e2l <- -Inf })
    # Define gamma values
    gammas <- seq(ind, ind+1, length=18)
    gammas <- gammas[-c(1, length(gammas))]
    # Calculate rejections for each gamma
    e <- vector("list", length(gammas))
    R.obs_new <- R0_new <- vector(length=length(gammas))
    for(i in seq_along(gammas)) {
      if(alternative != "greater") LB <- e1l + (gammas[i] - ind)*(e2l-e1l) else LB <- -Inf
      if(alternative != "less") UB <- e1u - (gammas[i] - ind)*(e1u-e2u) else UB <- Inf
      R.obs_new[i] <- noutside(obs_curve, LB, UB, alternative)
      R0_new[i] <- nr*mult*gammas[i]/nsim
      e[[i]] <- list(LB=LB, UB=UB)
    }
  }
  #-- Find the gamma
  find_gamma <- function(ind_new) {
    if(ind_new == 0) {
      if(ind != 0) {
        gamma <- ind
        e <- list(LB=e1l, UB=e1u)
        r1 <- ind_Robs # Robs[ind]
      }
      else {
        warning("Envelope was not found. Returning the largest envelope.")
        gamma <- gammas[1]
        e <- e[[1]]
        r1 <- R.obs_new[1]
        failed <- TRUE
      }
    }
    else {
      gamma <- gammas[ind_new]
      e <- e[[ind_new]]
      r1 <- R.obs_new[ind_new]
    }
    list(gamma=gamma, e=e, r1=r1, failed=failed)
  }
  # Indices for the three methods
  ind_new_ST <- find_ind(R.obs_new, R0_new, limit)

  list(ST=find_gamma(ind_new_ST))
}

# Specifications for the three steps of the ATSE and IATSE algorithms

limit1 <- function(algorithm) {
  switch(algorithm,
         ATSE = {
           function(alpha) alpha/(1+alpha)
         },
         IATSE = {
           function(alpha) alpha
         })
}

limit2 <- function(algorithm) {
  switch(algorithm,
         ATSE = {
           function(alpha, pi0) alpha/pi0
         },
         IATSE = {
           function(alpha, pi0) alpha/vlast(pi0)
         })
}

get_pi0est <- function(algorithm, alternative, nsim, nr) {
  switch(algorithm,
         ATSE = {
           function(r1, gamma=NULL) { # r1 (gamma ignored) from step1
             (nr-r1)/nr
           }
         },
         IATSE = {
           function(r1, gamma=NULL) { # r1 & gamma from step1
             if(alternative == "two.sided") gammastar <- 2*gamma/nsim
             else gammastar <- gamma/nsim
             min(1, ((nr-r1)/nr) / (1-gammastar))
           }
         })
}

find_ind_local <- function(FDRest = "ST") {
  switch(FDRest,
         ST = {
           function(rejs, limit) find_ind(rejs$R.obs, rejs$R0, limit)
         }
  )
}

general_FDRenvelope_algorithm <- function(algorithm = c("ATSE", "IATSE"), FDRest = "ST",
                                          alpha, alternative, curve_set,
                                          llimit=NULL, ulimit=NULL,
                                          rejs1=NULL) {
  pi0 <- 1
  gammastar <- step2 <- NULL
  nr <- curve_set_narg(curve_set)
  nsim <- curve_set_nfunc(curve_set) - 1

  limit1 <- limit1(algorithm) # Step 1. limit
  limit2 <- limit2(algorithm) # Step 3. limit
  pi0est <- get_pi0est(algorithm, alternative, nsim, nr) # Step 2. pi0 estimation
  find_ind_local <- find_ind_local(FDRest) # Find the index from rejs (R.obs, R0) for which the fdr is below the given limit

  # Calculate E(R0), R.obs, E(Q)
  if(is.null(rejs1)) rejs1 <- fdr_rejections_rank(curve_set, alternative)
  # Find the largest rejection region G* for which E(Q) <= limit1(alpha) (alpha/(1+alpha) for ATSE; alpha for IATSE)
  ind1 <- find_ind_local(rejs1, limit1(alpha))
  step1 <- fdr_rejections_between_two_rank(ind1, rejs1$R.obs[ind1], limit1(alpha),
                                           curve_set, alternative,
                                           llimit, ulimit) # pi0=1
  # Number of rejected hypotheses
  r1 <- step1[[FDRest]]$r1
  if(!(r1==0 | r1==nr) | (algorithm == "IATSE" & r1==nr)) {
    # Estimate \pi_0
    pi0 <- pi0est(r1, step1[[FDRest]]$gamma)
    # Find the largest rejection region G for which E(R0) / max(R.obs, 1) <= alpha/pi0
    ind2 <- find_ind_local(rejs1, limit2(alpha, pi0))
    step2 <- fdr_rejections_between_two_rank(ind2, rejs1$R.obs[ind2], limit2(alpha, pi0),
                                             curve_set, alternative,
                                             llimit, ulimit)
  }
  else # Take the first step envelope as the envelope
    step2 <- step1

  list(step2=step2[[FDRest]], pi0=pi0, rejs=rejs1, step1=step1[[FDRest]]) # what else? Note: rejs includes dataranks
}


# Calculates the fdr envelope for one curve_set
individual_fdr_envelope <- function(curve_set, alpha = 0.05,
                                    alternative = c("two.sided", "less", "greater"),
                                    algorithm=c("IATSE", "ATSE"),
                                    savefdr=FALSE, lower=NULL, upper=NULL) {
  if(!is.numeric(alpha) || (alpha < 0 | alpha >= 1)) stop("Unreasonable value of alpha.")
  alternative <- match.arg(alternative)
  algorithm <- match.arg(algorithm)
  isenvelope <- inherits(curve_set, "envelope")
  picked_attr <- pick_attributes(curve_set, alternative=alternative) # Saving for attributes / plotting purposes

  curve_set <- convert_envelope(curve_set)
  if(!curve_set_is1obs(curve_set)) stop("The curve_set must contain an observed function.")
  obs_curve <- curve_set_1obs(curve_set)

  nr <- curve_set_narg(curve_set)
  if(!is.null(lower)) {
    if(!is.numeric(lower) & is.finite(lower)) stop("Invalid lower.")
    if(length(lower) == 1) lower <- rep(lower, times=nr)
    else if(length(lower) != nr) stop("Invalid lower.")
  }
  if(!is.null(upper)) {
    if(!is.numeric(upper) & is.finite(upper)) stop("Invalid upper.")
    if(length(upper) == 1) upper <- rep(upper, times=nr)
    else if(length(upper) != nr) stop("Invalid upper.")
  }

  # Go through the algorithm to find the 'indices' for constructing envelopes
  inds <- general_FDRenvelope_algorithm(algorithm=algorithm, FDRest="ST",
                                        alpha, alternative, curve_set,
                                        llimit=lower, ulimit=upper)

  res <- make_envelope_object(type=algorithm, curve_set, LB=inds$step2$e$LB, UB=inds$step2$e$UB,
                              T_0=rowMeans(curve_set[['funcs']][,-1]),
                              picked_attr=picked_attr, isenvelope=isenvelope,
                              Malpha=NULL, alpha=alpha, distance=NULL)
  # Change the "method" attribute
  attr(res, "method") <- "FDR envelope"

  class(res) <- c("fdr_envelope", class(res))
  attr(res, "rej") <- obs_curve < inds$step2$e$LB | obs_curve > inds$step2$e$UB
  attr(res, "gamma") <- inds$step2$gamma
  attr(res, "pi0") <- inds$pi0
  if(savefdr) {
    attr(res, "R0") <- inds$rejs$R0
    attr(res, "R") <- inds$rejs$R.obs
  }
  res
}

#' The FDR envelope
#'
#' Calculate the FDR envelope based on the ATSE or IATSE algorithm
#' of Mrkvička and Myllymäki (2023).
#'
#'
#' Typical use of this function is through other functions.
#' \code{fdr_envelope(cset)} is the same as \code{global_envelope_test(cset, typeone = "fdr")}.
#' Functions such as \code{\link{graph.fanova}}, \code{\link{graph.flm}}, \code{\link{frank.flm}}
#' allow to use the FDR control by specifying \code{typeone = "fdr"} appropriately
#' (passing this to \code{global_envelope_test}).
#'
#' @inheritParams global_envelope_test
#' @inheritParams forder
#' @param algorithm The algorithm for the computation of the FDR envelope.
#' Either "IATSE" or "ATSE" standing for the iteratively adaptive two-stage
#' envelope and the adaptive two-stage envelope, respectively, see Mrkvička and Myllymäki (2023).
#' @param lower A single number (or a vector of suitable length) giving a lower bound
#' for the functions. Used only for the extension of the FDR envelope.
#' @param upper A single number (or a vector of suitable length) giving an upper bound
#' for the functions. Used only for the extension of the FDR envelope.
#' @export
#' @references
#' Mrkvička and Myllymäki (2023). False discovery rate envelopes. Statistics and Computing 33, 109. https://doi.org/10.1007/s11222-023-10275-7
#' @examples
#' # A GLM example
#' data(rimov)
#' \donttest{nsim <- 1000 # Number of simulations}
#' \dontshow{nsim <- 20}
#' res <- graph.flm(nsim=nsim,
#'                  formula.full = Y~Year,
#'                  formula.reduced = Y~1,
#'                  curve_sets = list(Y=rimov),
#'                  factors = data.frame(Year = 1979:2014),
#'                  GET.args = list(typeone = "fdr"))
#' plot(res)
#'
fdr_envelope <- function(curve_sets, alpha = 0.05,
                         alternative = c("two.sided", "less", "greater"),
                         algorithm = c("IATSE", "ATSE"),
                         lower = NULL, upper = NULL) {
  # Consider the case of several curves
  if(!is_a_single_curveset(curve_sets)) {
    if(length(curve_sets) > 1) { # Combine
      nr <- lapply(curve_sets, curve_set_narg)
      nfuns <- length(nr)
      is2d <- curve_set_is2d(curve_sets[[1]])
      csetnames <- names(curve_sets)
      curve_sets <- combine_curve_sets(curve_sets, equalr=FALSE) # Allow also different numbers of hypotheses
      # Test based on the combined sets
      res <- individual_fdr_envelope(curve_sets, alpha=alpha,
                                     alternative=alternative,
                                     algorithm=algorithm, savefdr=FALSE,
                                     lower=lower, upper=upper)
      # Transform the envelope to a combined envelope
      idx <- rep(1:nfuns, times=nr)
      # Split the envelopes to the original groups
      res_ls <- split(res, f=idx)
      # Redefine attributes method and rej
      rejs <- split(attr(res, "rej"), f=idx)
      for(i in 1:nfuns) {
        attr(res_ls[[i]], "method") <- paste0("1/", nfuns, "th of a combined FDR envelope")
        attr(res_ls[[i]], "rej") <- rejs[[i]]
      }
      # Set unreasonable attributes of individuals sets of curves to NULL
      anames <- c("alpha", "R0", "R", "pi0")
      anames <- anames[anames %in% names(attributes(res_ls[[1]]))]
      for(name in anames) {
        for(i in 1:nfuns) attr(res_ls[[i]], name) <- NULL
      }
      mostattributes(res_ls) <- attributes(res)
      attr(res_ls, "row.names") <- NULL
      if(!is.null(names(curve_sets))) names(res_ls) <- csetnames
      attr(res_ls, "method") <- "Combined FDR envelope"
      attr(res_ls, "nstep") <- 1
      class(res_ls) <- c("combined_fdr_envelope", "combined_global_envelope", "list")
      if(is2d)
        class(res_ls) <- c("combined_global_envelope2d", class(res_ls))
      return(res_ls)
    }
    else if(length(curve_sets) == 1) {
      curve_sets <- curve_sets[[1]] # unlist, and call individual_...
    }
    else
      stop("The given list of curve_sets is empty.")
  }
  return(individual_fdr_envelope(curve_sets, alpha=alpha,
                                 alternative=alternative,
                                 algorithm=algorithm, savefdr=FALSE,
                                 lower=lower, upper=upper))
}


#' Print method for the class 'fdr_envelope'
#'
#' @param x An 'fdr_envelope' object
#' @param ... Ignored.
#' @export
print.fdr_envelope <- function(x, ...) {
  cat(attr(x, "method"), "based on", attr(x, "einfo")$nsim, "simulations:\n",
      "Number of rejected hypotheses:", sum(attr(x, "rej")), "\n",
      "Number of accepted hypotheses:", length(attr(x, "rej")) - sum(attr(x, "rej")), "\n",
      "Total number of hypotheses:", length(attr(x, "rej")), "\n")
}
