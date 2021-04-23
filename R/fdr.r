# Returns the last component of a vector
vlast <- function(x) {
  x[length(x)]
}

#' Adaptive two stage method for FDR
#' @inheritParams fdr_envelope
#' @export
#' @examples
#' if(require("spatstat", quietly=TRUE)) {
#'   # CSR test
#'   #=========
#'   X <- unmark(spruces)
#'   # Number of simulations should be larger than
#'   alpha <- 0.10
#'   1 / (alpha/(1+alpha) / 513) # 513 is the number of r-values
#'   \donttest{nsim <- 6000 # Number of simulations}
#'   \dontshow{nsim <- 20}
#'   env <- envelope(X, fun="Lest", nsim=nsim,
#'                   savefuns=TRUE, # save the functions
#'                   correction="translate", # edge correction for L
#'                   transform = expression(.-r), # centering
#'                   simulate=expression(runifpoint(ex=X))) # Simulate CSR
#'   res1 <- fdr_ATS(env, alpha=alpha)
#'   plot(res1$pvalues, col=res1$rej*1+1) # red = rejected
#' }
fdr_ATS <- function(curve_set, alpha = 0.05,
                    alternative = c("two.sided", "less", "greater")) {
  isenvelope <- inherits(curve_set, "envelope")
  alternative <- match.arg(alternative)
  if(!is.numeric(alpha) || (alpha < 0 | alpha > 1)) stop("Unreasonable value of alpha.")
  curve_set <- convert_envelope(curve_set)
  if(!curve_set_is1obs(curve_set)) stop("The curve_set must contain an observed function.")
  data_and_sim_curves <- data_and_sim_curves(curve_set) # Include the data function
  nr <- curve_set_narg(curve_set) # m in the fdr paper

  # Calculate pointwise p-values for each argument value (r)
  calc_pointwiserank <- find_calc_pointwiserank("rank", alternative)
  pps <- vector(length=nr)
  for(i in 1:nr) {
    pranks <- calc_pointwiserank(data_and_sim_curves[,i])
    pps[i] <- estimate_p_value(x=-pranks[1], sim_vec=-pranks[-1], ties = 'conservative')
  }
  # Linear step-up procedure at level alpha/(1+alpha)
  pps.sorted <- sort(pps)
  rej.sorted <- which(pps.sorted <= (1:nr)*((alpha/(1+alpha))/nr))
  r1 <- vlast(rej.sorted)
  pi0 <- NULL
  if(length(rej.sorted) < 1) rej <- rep(FALSE, times=nr)
  else if(r1 == nr) rej <- rep(TRUE, times=nr)
  else {
    pi0 <- (nr-r1)/nr
    # Linear step-up procedure at level alpha/pi0
    rej.sorted <- which(pps.sorted <= (1:nr)*((alpha/pi0)/nr))
    r1 <- vlast(rej.sorted)
    rej <- pps <= pps.sorted[r1]
  }
  list(pvalues=pps, rej=rej, pi0=pi0)
}


# Helper functions for FDR envelopes

# Calculate E(R0) and R.obs for the rank envelopes
fdr_rejections_rank <- function(curve_set, alternative, curve_set2=NULL, fast=TRUE) {
  curve_set <- convert_envelope(curve_set)

  obs_curve <- curve_set_1obs(curve_set)
  sim_curves <- data_or_sim_curves(curve_set) # Exclude the data function
  #all_curves <- data_and_sim_curves(curve_set) # Include the data function
  nsim <- curve_set_nfunc(curve_set) - 1
  nr <- curve_set_narg(curve_set) # m in the fdr paper

  if(!is.null(curve_set2)) {
    curve_set2 <- convert_envelope(curve_set2)

    # From the second set of simulations, calculate R0 (the numbers of rejections under the null)
    sim_curves2 <- data_or_sim_curves(curve_set2) # Exclude the data function (if it exists)
    # FIXME, there should be check that curve_set and curve_set2 correspond, and obs_curve is the same
    all_curves2 <- data_and_sim_curves(curve_set2) # The data function and the second set of simulated functions
    nsim2 <- curve_set_nfunc(curve_set2) - 1

    if(fast) {
      all_curves <- rbind(obs_curve, sim_curves2, sim_curves)

      # Calculate pointwise ranks for each argument value (r)
      calc_pointwiserank <- find_calc_pointwiserank("rank", alternative)
      for(i in 1:nr) {
        sim_curves[,i] <- calc_pointwiserank(sim_curves[,i]) # overwriting curves by their ranks
        all_curves[,i] <- calc_pointwiserank(all_curves[,i])
        all_curves2[,i] <- calc_pointwiserank(all_curves2[,i])
      }
      # Ranks of the simulated curves from which the envelopes are constructed
      simranks <- sim_curves
      # Ranks of the data curve and all simulated curves from curve_set2 among all the curves
      datasim2ranks <- all_curves[1:(1+nsim2),]
      # Ranks of the data curve and all simulated curves from curve_set2 among curve_set curves
      datasim2ranks <- datasim2ranks - (all_curves2 - 1)

      # Count the numbers of rejections R0_j for data and for all simulations j, j=1,...,nsim2 and for all envelopes l=1,2,3,...,floor(nsim/2)
      max_l <- floor(nsim/2)
      rejs <- matrix(nrow=1+nsim2, ncol=max_l)
      for(i in 1:max_l) {
        rejected <- datasim2ranks <= i
        rejs[,i] <- rowSums(rejected)
      }
      R0 <- apply(rejs[-1,], MARGIN = 2, FUN=mean)
      R.obs <- rejs[1,]
    }
    else {
      # Calculate pointwise ranks for each argument value (r) for the functions in the first set of curves
      calc_pointwiserank <- find_calc_pointwiserank("rank", alternative)
      simranks <- sim_curves
      for(i in 1:nr) {
        simranks[,i] <- calc_pointwiserank(sim_curves[,i]) # NOT overwriting curves by their ranks
      }

      # Calculate ranks for data curve
      all_curves <- rbind(obs_curve, sim_curves)
      for(i in 1:nr) {
        all_curves[,i] <- calc_pointwiserank(all_curves[,i]) # overwriting curves by their ranks. FIXME(for computational reasons)? Only the rank of data needed here.
      }
      dataranks <- all_curves[1,]

      # and all curves from the second set of simulations
      sim2_ranks <- sim_curves2
      for(j in 1:nsim2) {
        all_curves2 <- rbind(sim_curves2[j,], sim_curves)
        for(i in 1:nr) {
          all_curves2[,i] <- calc_pointwiserank(all_curves2[,i]) # overwriting curves by their ranks. FIXME(for computational reasons)? Only the rank of data needed here.
        }
        sim2_ranks[j,] <- all_curves2[1,]
      }

      # Count the numbers of rejections
      max_l <- floor(nsim/2)
      rejs <- matrix(nrow=nsim2, ncol=max_l)
      for(i in 1:max_l) {
        rejected <- sim2_ranks <= i
        rejs[,i] <- rowSums(rejected)
      }
      R0 <- apply(rejs, MARGIN = 2, FUN=mean)
      R.obs <- sapply(1:floor(nsim/2), FUN = function(i) sum(dataranks <= i))
    }
  }
  else { # The case of just one set of curves (curve_set2 = NULL)
    all_curves <- data_and_sim_curves(curve_set) # Include the data function

    # Calculate pointwise ranks for each argument value (r)
    calc_pointwiserank <- find_calc_pointwiserank("erl", alternative)
    for(i in 1:nr) {
      sim_curves[,i] <- calc_pointwiserank(sim_curves[,i]) # overwriting curves by their ranks
      all_curves[,i] <- calc_pointwiserank(all_curves[,i]) # overwriting curves by their ranks. FIXME(for computational reasons)? Only the rank of data needed here.
    }
    simranks <- sim_curves
    dataranks <- all_curves[1,]

    # Count the numbers of rejections R0_i for all simulations i and for all envelopes l=1,2,3...,floor(nsim/2)
    max_l <- floor(nsim/2)
    rejs <- matrix(nrow=nsim, ncol=max_l)
    for(i in 1:max_l) {
      rejected <- simranks < i
      rejs[,i] <- rowSums(rejected)
    }
    R0 <- apply(rejs, MARGIN = 2, FUN=mean)
    R.obs <- sapply(1:floor(nsim/2), FUN = function(i) sum(dataranks <= i))
  }

  list(R0=R0, R.obs=R.obs)
}

fdr_step13 <- function(R.obs, R0) {
  R0 / pmax(R.obs, 1)
}

# Find the index from rejs for which the fdr is below the given limit
find_ind <- function(R.obs, R0, limit, leftright) {
  switch(leftright,
    right = { # The last one that is still within the limit
      max(0, vlast(which(fdr_step13(R.obs, R0) <= limit)))
    },
    left = { # The last one before fdr first goes above the limit
      max(0, (which(fdr_step13(R.obs, R0) > limit))[1]) - 1
    }
  )
}
# Find first the index, i.e. integer value
# Then find the gamma between suitable indices
find_gamma_and_newrankenvelope <- function(R.obs, R0, limit, curve_set, curve_set2, leftright="right") { # FIXME: Only two.sided implemented!
  ind <- find_ind(R.obs, R0, limit, leftright)
  obs_curve <- curve_set_1obs(curve_set)
  sim_curves2 <- data_or_sim_curves(curve_set2) # The second set of simulations; calculate the rejections for these
  nsim2 <- dim(sim_curves2)[1]
  if(ind == 0) { # The value of gamma is between 0 and 1
    # Find the envelope of rank 1; THE envelope is wider than this
    e1 <- fdr_envelope_rank(1, curve_set)
    # Define the gamma values
    gammas <- seq(0, 1, length=130)
    gammas <- gammas[-c(1, length(gammas))]
    # Calculate rejections for each gamma
    e <- list()
    R.obs_new <- R0_new <- vector(length=128)
    for(i in 1:length(gammas)) {
      LB <- e1$LB + log(gammas[i]) * (e1$UB - e1$LB)
      UB <- e1$UB - log(gammas[i]) * (e1$UB - e1$LB)
      R0_new[i] <- mean(sapply(1:nsim2, FUN = function(j) { sum(sim_curves2[j,] < LB | sim_curves2[j,] > UB) }))
      R.obs_new[i] <- sum(obs_curve < LB | obs_curve > UB)
      e[[i]] <- list(LB=LB, UB=UB)
    }
  }
  else { # The value of gamma lies between ind and ind + 1
    # Calculate the envelopes between which THE envelope lies
    sim_curves <- data_or_sim_curves(curve_set)
    nr <- curve_set_narg(curve_set)
    nsim <- dim(sim_curves)[1]
    e1l <- e1u <- e2l <- e2u <- array(0, nr)
    for(i in 1:nr){
      Hod <- sort(sim_curves[,i])
      e1l[i] <- Hod[ind]
      e1u[i] <- Hod[nsim-ind+1]
      e2l[i] <- Hod[ind+1]
      e2u[i] <- Hod[nsim-ind]
    }
    # Define gamma values
    gammas <- seq(ind, ind+1, length=18)
    gammas <- gammas[-c(1, length(gammas))]
    # Calculate rejections for each gamma
    e <- list()
    R.obs_new <- R0_new <- vector(length=length(gammas))
    for(i in 1:length(gammas)) {
      LB <- e1l + (gammas[i] - ind)*(e2l-e1l)
      UB <- e1u - (gammas[i] - ind)*(e1u-e2u)
      R0_new[i] <- mean(sapply(1:nsim2, FUN = function(j) { sum(sim_curves2[j,] < LB | sim_curves2[j,] > UB) }))
      R.obs_new[i] <- sum(obs_curve < LB | obs_curve > UB)
      e[[i]] <- list(LB=LB, UB=UB)
    }
  }
  # Find the gamma
  ind_new <- find_ind(R.obs_new, R0_new, limit, leftright)
  if(ind_new == 0) {
    gamma <- ind
    e <- list(LB=e1l, UB=e1u)
    r1 <- R0[ind]
  }
  else {
    gamma <- gammas[ind_new]
    e <- e[[ind_new]]
    r1 <- R.obs_new[ind_new]
  }

  #list(ind, R.obs_new=R.obs_new, R0_new=R0_new, gammas=gammas, ind_new=ind_new, e=e, r1=r1)
  list(gamma=gamma, e=e, r1=r1)
}

atsge_algorithm <- function(rejs, alpha, curve_set, curve_set2, leftright="right") {
  pi0 <- step2 <- NULL
  nr <- curve_set_narg(curve_set)
  # Find the largest rejection region G* for which E(R0) / max(R.obs, 1) <= alpha/(1+alpha)
  step1 <- find_gamma_and_newrankenvelope(rejs$R.obs, rejs$R0, alpha/(1+alpha), curve_set, curve_set2, leftright)
  # Number of rejected hypotheses
  r1 <- step1$r1
  if(!(r1==0 | r1==nr)) {
    # Estimate \pi_0
    pi0 <- (nr-r1)/nr
    # Find the largest rejection region G for which E(R0) / max(R.obs, 1) <= alpha/pi0
    step2 <- find_gamma_and_newrankenvelope(rejs$R.obs, rejs$R0, alpha/pi0, curve_set, curve_set2, leftright)
  }
  else # Take the first step envelope as the envelope
    step2 <- step1
  list(step2=step2, pi0=pi0)
}

iatsge_algorithm <- function(rejs, alpha, curve_set, curve_set2, alternative, leftright="right") {
  pi0 <- gammastar <- step2 <- NULL
  nr <- curve_set_narg(curve_set)
  nsim <- curve_set_nfunc(curve_set) - 1
  # Find the largest rejection region G* for which E(R0) / max(R.obs, 1) <= alpha
  if(alternative != "two.sided") stop("One-sided alternative not implemented yet.")
  step1 <- find_gamma_and_newrankenvelope(rejs$R.obs, rejs$R0, alpha, curve_set, curve_set2, leftright)
  # Number of rejected hypotheses
  r1 <- step1$r1
  if(!(r1==0 | r1==nr)) {
    # Estimate \pi_0
    #####pi0 <- min( 1, (nr-r1+(nr-r1+nr*alpha)*alpha)/nr )
    #if(type == "rank")
    if(alternative == "two.sided") gammastar <- 2*step1$gamma/nsim
    else gammastar <- step1$gamma/nsim
    pi0[1] <- min(1, (nr-r1+nr*gammastar)/nr)
    pi0[2] <- min(1, (nr-r1+pi0[1]*nr*gammastar)/nr)
    j <- 2
    while(pi0[j-1]-pi0[j] > 10^-10) {
      j <- j+1
      pi0[j] <- min(1, (nr-r1+pi0[j-1]*nr*gammastar)/nr)
    }
    # Find the largest rejection region G for which E(R0) / max(R.obs, 1) <= alpha/pi0
    step2 <- find_gamma_and_newrankenvelope(rejs$R.obs, rejs$R0, alpha/vlast(pi0), curve_set, curve_set2, leftright)
  }
  else # Take the first step envelope as the envelope
    step2 <- step1
  list(step2=step2, pi0=pi0, gammastar=gammastar)
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

#' The FDR envelope
#'
#' Calculate the FDR envelope based on ATSGE or IATSGE algorithm.
#'
#' @param curve_set A \code{curve_set} (see \code{\link{create_curve_set}}) or an
#'   \code{envelope} object of \pkg{spatstat} containing the observed function
#'   and the functions from which the envelope is to be constructed.
#' @param curve_set2 A \code{curve_set} giving another set of simulations from which
#' the numbers of rejections under the null are calculated.
#' @param type The type of the global envelope with current options for 'rank', 'erl', 'cont', 'area'.
#' See \code{\link{forder}} for description of the measures.
#' @inheritParams global_envelope_test
#' @param method Either "ATSGE" or "IATSGE" standing for adaptive two stage
#'   global envelope, and iteratively adaptive two stage global envelope methods,
#'   respectively.
#' @export
#' @examples
#' # Just for fun
#' if(require("spatstat", quietly=TRUE)) {
#'   # CSR test
#'   #=========
#'   X <- unmark(spruces)
#'   \donttest{nsim <- 2000 # Number of simulations}
#'   \dontshow{nsim <- 20}
#'   env1 <- envelope(X, fun="Lest", nsim=nsim,
#'                   savefuns=TRUE, # save the functions
#'                   correction="translate", # edge correction for L
#'                   transform = expression(.-r), # centering
#'                   simulate=expression(runifpoint(ex=X))) # Simulate CSR
#'   env2 <- envelope(X, fun="Lest", nsim=nsim,
#'                   savefuns=TRUE, # save the functions
#'                   correction="translate", # edge correction for L
#'                   transform = expression(.-r), # centering
#'                   simulate=expression(runifpoint(ex=X))) # Simulate CSR
#'   res <- fdr_envelope(env1, env2)
#'   plot(res)
#' }
#'
#' # A GLM example
#' data(rimov)
#' \donttest{nsim <- 1000 # Number of simulations}
#' \dontshow{nsim <- 20}
#' sims1 <- graph.flm(nsim=nsim,
#'                  formula.full = Y~Year,
#'                  formula.reduced = Y~1,
#'                  curve_sets = list(Y=rimov), factors = data.frame(Year = 1979:2014),
#'                  savefuns = TRUE)
#' sims2 <- graph.flm(nsim=nsim,
#'                  formula.full = Y~Year,
#'                  formula.reduced = Y~1,
#'                  curve_sets = list(Y=rimov), factors = data.frame(Year = 1979:2014),
#'                  savefuns = TRUE)
#' res <- fdr_envelope(attr(sims1, "simfuns")[[1]], attr(sims2, "simfuns")[[1]])
#' plot(res)
#'
fdr_envelope <- function(curve_set, curve_set2, type = "rank", alpha = 0.05,
                         alternative = c("two.sided", "less", "greater"),
                         method = c("IATSGE", "ATSGE")) {
  isenvelope <- inherits(curve_set, "envelope")
  type <- match.arg(type)
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  if(!is.numeric(alpha) || (alpha < 0 | alpha > 1)) stop("Unreasonable value of alpha.")
  picked_attr <- pick_attributes(curve_set, alternative=alternative) # saving for attributes / plotting purposes
  curve_set <- convert_envelope(curve_set)
  curve_set2 <- convert_envelope(curve_set2)
  if(!curve_set_is1obs(curve_set)) stop("The curve_set must contain an observed function.")
  obs_curve <- curve_set_1obs(curve_set)
  nsim <- curve_set_nfunc(curve_set) - 1
  nr <- curve_set_narg(curve_set) # m in the fdr paper
  if(curve_set_is1obs(curve_set2))
    if(!identical(obs_curve, curve_set_1obs(curve_set2)))
      stop("The observed functions in the two sets of curves differ.")
  if(!identical(nr, curve_set_narg(curve_set2)))
    stop("The number of hypotheses in the two sets of curves differ.")
  # FIXME Add further checks for the two curve sets

  # Calculate E(R0) and R.obs
  rejs <- fdr_rejections_rank(curve_set, alternative)

  # Go through the algorithm to find the 'indices' for constructing envelopes
  switch(method,
         ATSGE = {
           inds <- atsge_algorithm(rejs, alpha, curve_set, curve_set2)
         },
         IATSGE = {
           inds <- iatsge_algorithm(rejs, alpha, curve_set, curve_set2, alternative)
         })

  res <- make_envelope_object(type, curve_set, LB=inds$step2$e$LB, UB=inds$step2$e$UB,
                              T_0=rowMeans(curve_set[['funcs']][,-1]),
                              picked_attr=picked_attr, isenvelope=isenvelope,
                              Malpha=NULL, alpha=alpha, distance=NULL)
  # Change the "method" attribute
  attr(res, "method") <- "FDR envelope"

  class(res) <- c("fdr_envelope", class(res))
  attr(res, "R0") <- rejs$R0
  attr(res, "R") <- rejs$R.obs
  attr(res, "pi0") <- inds$pi0
  # FIXME, add further attributes from 'inds'
  attr(res, "rej") <- obs_curve < inds$step2$e$LB | obs_curve > inds$step2$e$UB
  res
}


#' Print method for the class 'fdr_envelope'
#'
#' @param x An 'fdr_envelope' object
#' @param ... Ignored.
#' @export
print.fdr_envelope <- function(x, ...) {
  cat(attr(x, "method"), "based on", attr(x, "einfo")$nsim, "simulations:\n",
      "Number of rejected hypotheses:", sum(attr(x, "rej")), "\n",
      "Number of accepted hypotheses:", length(x[['r']]) - sum(attr(x, "rej")), "\n",
      "Total number of hypotheses:", length(x[['r']]), "\n")
}
