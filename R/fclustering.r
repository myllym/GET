#' Functional clustering
#'
#' Functional clustering based on a specified measure.
#' The options of the measures can be found in \code{\link{central_region}}.
#'
#'
#' Functional clustering joins the list of \code{curve_set} objects in one \code{curve_set} with long functions and
#' applies on the differences of all functions the specified measure. This provides a dissimilarity matrix
#' which is used in partitioning around medoids procedure. The resulting clusters can then be shown by plotting
#' the function respectively for each \code{curve_set}. Thus for each \code{curve_set}, the panel with all the medoids
#' is shown followed by all clusters represented by central region, medoid and all curves belonging to it, when
#' the result object is plotted.
#'
#' If there are less than three curves in some of the groups, then the central region is not plotted.
#' This leads to a warning message from ggplot2.
#'
#' @param curve_sets A \code{curve_set} object or a list of \code{curve_set} objects to which
#' the functional clustering is to be applied. If list of \code{curve_set} objects is provided,
#' then the joined functional clustering is applied, which provides an equal weight combination
#' of \code{curve_set} objects, if the \code{curve_set} objects contain the same numbers of elements
#' (same lengths of vector \eqn{r}{r}).
#' @param k The number of clusters.
#' @param type The measure which is used to compute the dissimilarity matrix. The preferred options
#' are \code{"area"} and \code{"st"}, but \code{"erl"} and \code{"cont"} can be also used with caution.
#' @param ... Additional parameters to be passed to \code{\link{central_region}},
#' which is responsible for calculating the central region (global envelope)
#' on which the functional clustering is based.
#'
#' @return An object having the class \code{fclust}, containing
#' \itemize{
#' \item curve_sets = The set(s) of functions determined for clustering
#' \item k = Number of clusters
#' \item type = Type of clustering method
#' \item triangineq = The proportion of combinations of functions which satisfies the triangular inequality.
#' The triangular inequality must hold to ensure the chosen measure forms a metric. In some weird cases
#' it does not hold for ‘area’ measure, therefore this check is provided to ensure the data forms metric
#' with the ‘area’ measure. The triangineq must be 1 to ensure the inequality holds for all functions.
#' \item dis = The joined dissimilarity matrix
#' \item pam = Results of the partitioning around medoids (pam) method applied on the joined functions
#' with the dissimilarity matrix (dis). See \code{\link{pam}}.
#' }
#' @references
#' Dai, W., Athanasiadis, S., Mrkvička, T. (2021) A new functional clustering method with joined dissimilarity sources and graphical interpretation. Journal of Multivariate Analysis.
#'
#' @export
#' @seealso \code{\link{central_region}}, \code{\link{plot.fclust}}
#' @examples
#' # Read raw data from population growth rdata
#' # with countries over million inhabitants
#' data("popgrowthmillion")
#' \dontshow{
#' popgrowthmillion <- popgrowthmillion[1:10, 1:50]
#' }
#' # Create centred data
#' m <- apply(popgrowthmillion, 2, mean) # Country-wise means
#' cpopgrowthmillion <- popgrowthmillion
#' for(i in 1:dim(popgrowthmillion)[1]) {
#'   cpopgrowthmillion[i,] <- popgrowthmillion[i,] - m
#' }
#'
#' # Create scaled data
#' t2 <- function(v) { sqrt(sum(v^2)) }
#' s <- apply(cpopgrowthmillion, 2, t2)
#' spopgrowthmillion <- popgrowthmillion
#' for(i in 1:dim(popgrowthmillion)[1]) {
#'   spopgrowthmillion[i,] <- cpopgrowthmillion[i,]/s
#' }
#'
#' # Create curve sets
#' r <- 1951:2015
#' \dontshow{
#' r <- r[1:10]
#' }
#' cset1 <- create_curve_set(list(r = r, obs = popgrowthmillion))
#' cset2 <- create_curve_set(list(r = r, obs = spopgrowthmillion))
#' csets <- list(Raw = cset1, Shape = cset2)
#'
#' # Functional clustering with respect to joined "st" difference measure
#' # and "joined" central regions of each group
#' res <- fclustering(csets, k=3, type="area")
#' p <- plot(res, plotstyle = "marginal", coverage = 0.5)
#' p[[1]] # Central functions
#' p[[2]] # Groups: central functions and regions
#' # To collect the two figures into one use, e.g., patchwork:
#' if(require("patchwork", quietly=TRUE)) {
#'   p[[1]] + p[[2]] + plot_layout(widths = c(1, res$k))
#' }
#' # Silhouette plot of pam
#' plot(res$pam)
#' @importFrom stats dist
#' @importFrom cluster pam
#' @importFrom utils combn
fclustering <- function(curve_sets, k, type = c("area", "st", "erl", "cont"), ...) {
  # Check curve_sets. Note: k is checked by pam
  if(!(length(class(curve_sets)) == 1 && class(curve_sets) == "list")) {
    curve_sets <- list(curve_sets) # Make a list of a single curve_set to treat it similarly as several sets of curves
  }
  # Convert (e.g. fdata) to curve_sets and check the content
  curve_sets <- lapply(curve_sets, FUN=convert_to_curveset)
  if(any(sapply(curve_sets, FUN=curve_set_is1obs))) stop("Some or all sets of curves have one observed function.")
  # Check consistency of the given curve sets
  checkequal <- function(f) {
    all(sapply(curve_sets, FUN=function(curve_set) { f(curve_set) == f(curve_sets[[1]]) }))
  }
  if(!checkequal(curve_set_nfunc))
    stop("The numbers of functions in curve sets differ.")
  if(!checkequal(curve_set_narg))
    warning("The numbers of points where the curves are observed (r) differ in the given curve sets.")

  type <- match.arg(type)
  if(type %in% c("erl", "cont"))
    warning("Type is ", type, ". Check if the triangular inequality is equal to 1. \n",
            "If not, the measure does not form a metric. pam should not be used.")
  if(type == "st") {
    measure <- "max"; scaling <- "st"
  }
  else {
    measure <- type; scaling <- NULL
  }
  # Number of curve_sets
  n <- length(curve_sets)

  # Compute the dissimilarity matrix
  # lstats = original statistics, jstats = joined statistics
  lstats <- list()
  jr <- NULL
  jstats <- NULL
  for(j in 1:n) {
    nr <- curve_set_narg(curve_sets[[j]])
    nfunc <- curve_set_nfunc(curve_sets[[j]])
    funcs <- curve_set_funcs(curve_sets[[j]])
    co <- combn(nfunc, 2, simplify = FALSE) # All combinations
    stats <- array(0, c(2*length(co), nr))
    for(i in 1:length(co)) {
      stats[i,] <- funcs[, co[[i]][1]]-funcs[, co[[i]][2]]
    }
    for(i in 1:length(co)) {
      stats[length(co)+i, ] <- funcs[, co[[i]][2]]-funcs[, co[[i]][1]]
    }
    lstats[[j]] <- stats

    # Combined
    jr <- c(jr, curve_sets[[1]][['r']])
    jstats <- cbind(jstats, lstats[[j]]) # j = joint
  }
  # Add the zero curve
  jstats <- rbind(jstats, rep(0, length(jr))) # length(jr) is same as dim(jstats)[2]
  data <- create_curve_set(list(r=jr, obs=t(jstats))) # Non-increasing r-values provided, but not used

  # dist gives distances between the rows of a data matrix. each row = a function
  b <- dist(data_and_sim_curves(curve_sets[[1]])) # Or t(curve_set_funcs(curve_sets[[1]]))
  # Fill it in with other dissimilarity matrix

  if(measure == "max")
    fo <- forder(data, measure="max", scaling=scaling)
  else
    fo <- 1 - forder(data, measure=measure, scaling=scaling)
  b[1:length(b)] <- fo[1:length(b)] # for(i in 1:length(b)) { b[i] <- fo[i] }
  resultpamF <- pam(b, k=k)
  bb <- as.matrix(b)
  nr1 <- curve_set_narg(curve_sets[[1]])
  # Triangular inequality
  r <- NULL # Here r is a logical vector for something
  for(i in 1:(nr1-2)) {
    for(j in (i+1):(nr1-1)) {
      for(kk in (j+1):(nr1)) {
        r <- c(r, bb[i,j] + bb[j,kk] >= bb[i,kk] - 0.00000001)
      }
    }
  }

  res <- list(curve_sets=curve_sets, k=k, type=type,
              pam=resultpamF, dis=b,
              triangineq = sum(r)/length(r))
  class(res) <- "fclust"
  res
}

#' Print method for the class 'fclust'
#'
#' Print method for the 'fclust' objects returned by \code{\link{fclustering}}.
#' @param x A object of class 'fclust'.
#' @param ... Ignored.
#' @export
print.fclust <- function(x, ...) {
  cat("Functional clustering:\n")
  cat(" * Based on the measure \'", x$type, "\' at ", x$k, " clusters \n", sep="")
  cat(" * Triangular inequality is satisfied in", x$triangineq*100, "% of cases.\n")
  cat("The object contains:\n")
  cat("$pam - Results of the partitioning around medoids\n")
  cat("$pam$clustering - The clustering vector\n")
  cat("$dis - The dissimilarity matrix\n")
}

#' Plot method for the class 'fclust'
#'
#' Plot method for the 'fclust' objects returned by \code{\link{fclustering}}.
#'
#'
#' The clusters are shown respectively for each \code{curve_set}. Thus for each \code{curve_set}
#' the panel with all the medoids is shown followed by all clusters represented by central region,
#' medoid and all curves belonging to it.
#' 
#' For all sources, the function plots the deepest curves for all clusters and
#' the deepest curve of each cluster together with the desired central region and
#' all the curves of the group.
#' @param x An 'fclust' object.
#' @param ncol The maximum number of columns for the figures determining the dimensions of the graphical output.
#' Default to a one row figure.
#' @param ncol The number of columns in the graphical output, when there is just one set of curves
#' that has been ordered. If not given, \code{c(1, k+1)} is used, which gives all plots in one row.
#' For more sets of curves, the rows are fixed to correspond to the sets (one row for each set).
#' @param plotstyle The resulting central regions of clusters can be plotted by sorting the
#' appropriate \code{curve_set} only 'marginal' or by sorting the joined list of \code{curve_set}
#' objects 'joined'. If 'joined' is used the shown central regions corresponds to the joined
#' ordering used to cluster the functional data. If 'marginal' is used the shown central regions
#' do not correspond to the joined ordering used to cluster the functional data, but better express
#' the shape of cluster with respect to given \code{curve_set}.
#' @param coverage The coverage of central regions to be used to show the clusters.
#' @param ... Ignored.
#' @param nstep 1 or 2 for how to contruct a combined (joined) global envelope
#' if there are more than one sets of curves. Default to 1, if the numbers of points
#' where the curves are observed (r) are the same in each set, and 2 otherwise.
#'
#' @references
#' Dai, W., Athanasiadis, S., Mrkvička, T. (2021) A new functional clustering method with joined dissimilarity sources and graphical interpretation. Journal of Multivariate Analysis.
#' @export
#' @importFrom ggplot2 ggplot aes_ geom_line guides facet_wrap facet_grid vars
plot.fclust <- function(x, plotstyle = c("marginal", "joined"), coverage = 0.5,
                        nstep, ncol, ...) {
  csets <- x$curve_sets
  n <- length(csets)
  k <- x$k
  type <- x$type
  plotstyle <- match.arg(plotstyle)
  if(coverage < 0 | coverage > 1) stop("Unreasonable coverage.")
  if(missing(ncol)) ncol <- k+1
  if(!is.numeric(ncol) || length(ncol) != 1) stop("Unreasonable ncol.")
  if(type %in% "st") centralf <- min
  else centralf <- max
  if(missing(nstep)) {
    nstep <- 1
    checkequal <- function(f) {
      all(sapply(csets, FUN=function(curve_set) { f(curve_set) == f(csets[[1]]) }))
    }
    if(!checkequal(curve_set_narg))
      nstep <- 2
  }
  if(!(nstep %in% c(1,2))) stop("Invalid value of nstep.")

  if(!all(sapply(csets, FUN = curve_set_is1d)))
    stop("Plots are valid only for 1-dimensional curve_set objects.")

  resultpamF <- x$pam
  resultpamT <- x$pam$clustering

  nis <- as.numeric(table(resultpamT)) # Numbers of functions in the groups i, i = 1,...,k
  #- Calculate central regions (either)
  cr_all <- vector("list", n)
  if(plotstyle == "marginal") {
    for(j in 1:n) { # For all sets of curves
      cr <- vector("list", k)
      ff <- vector("list", k)
      # Calculate the central regions for each cluster
      for(i in 1:k) { # For all groups (within each set of curves)
        # Calculate the central regions for each cluster
        clusteri <- subset(csets[[j]], resultpamT==i)
        if(nis[i] > 3) { # At least 3 functions
          cr[[i]] <- central_region(clusteri, type=type, coverage=coverage)
        }
        else {
          cr[[i]] <- data.frame(r=clusteri[['r']], central=1, lo=NA, hi=NA)
        }
        # Redefine the central curve to be the medoid
        cr[[i]]$central <- curve_set_funcs(csets[[j]])[, resultpamF$id.med[i]]
      }
      cr_all[[j]] <- cr
    }
  }
  else { # "joined"
    cr <- vector("list", k)
    ff <- vector("list", k)
    for(i in 1:k) { # Calculate the central regions for each cluster
      if(nis[i] > 3) {
        clusteri <- lapply(csets, FUN=subset, subset=resultpamT==i) # Curves for the group i
        cr[[i]] <- central_region(clusteri, type=type, coverage=coverage, nstep=nstep) # Joint central region for the group i
        # Redefine the central curve to be the medoid
        if(n > 1) for(z in 1:n) cr[[i]][[z]]$central <- curve_set_funcs(csets[[z]])[, resultpamF$id.med[i]]
        else cr[[i]]$central <- curve_set_funcs(csets[[1]])[, resultpamF$id.med[i]]
      }
    }
    for(j in 1:n) cr_all[[j]] <- lapply(cr, FUN=function(z) { z[[j]] })
  }

  #- Setup for the plotting
  alt <- "two.sided"
  d <- list(xlab="", ylab="")
  labels <- paste("Group", 1:k)
  if(n < 2) {
    n_of_plots <- k + 1
    ncols_of_plots <- min(n_of_plots, ncol)
    nrows_of_plots <- ceiling(n_of_plots / ncols_of_plots)
    scales <- "fixed"
  }
  else scales <- "free_y"
  if(!is.null(names(csets))) csetnames <- names(csets)
  else csetnames <- paste("Curves", 1:n)

  #- df: Central functions and regions
  df <- NULL
  for(i in 1:n) {
    df_i <- combined_df_construction(cr_all[[i]], labels=labels)
    df_i$cset <- csetnames[i]
    df <- rbind(df, df_i)
  }
  df$cset <- factor(df$cset, levels=csetnames)
  #- df: All single functions
  funcs <- lapply(csets, FUN=curve_set_funcs)
  nams <- lapply(funcs, FUN=function(x) { n <- colnames(x); if(!is.null(n)) n else 1:ncol(x) })
  funcs.df <- do.call(rbind, lapply(1:length(funcs), FUN = function(i) {
    data.frame(r = rep(cr_all[[1]][[i]][['r']], times=ncol(funcs[[i]])),
               curves = c(funcs[[i]]),
               id = factor(rep(nams[[i]], each=length(cr_all[[1]][[i]][['r']])), levels = nams[[i]]),
               plotmain = factor(rep(paste("Group", resultpamT), each=length(cr_all[[1]][[i]][['r']])), levels = labels),
               cset = csetnames[i])
  }))
  #- Plot
  size <- 0.3
  df$Group <- df$plotmain
  p1 <- ( ggplot()
          + geom_line(data = df, aes_(x = ~r, y = ~curves, group = ~Group,
                                    linetype = ~Group), size = size, col='black')
          + set_envelope_legend_position()
          + labs(title = "", x = d$xlab, y = d$ylab)
  )
  if(n > 1) p1 <- ( p1 + facet_grid(cset~type, scales="free_y") )
  else p1 <- ( p1 + facet_wrap(~type, scales="free_y") )

  p2 <- ( ggplot()
         + geom_line(data = funcs.df, aes_(x = ~r, y = ~curves, group = ~id), col='grey70')
         + basic_stuff_for_fclustplot(df, d$xlab, d$ylab, main="", size=size)
         + guides(linetype = "none")
  )
  if(n > 1) p2 <- (p2 + facet_grid(cset~plotmain, scales=scales))
  else p2 <- (p2 + facet_wrap(~ plotmain, scales=scales,
                              nrow=nrows_of_plots, ncol=ncols_of_plots))

  #- Collect
  list(p1, p2)
}
