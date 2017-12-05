#' Year temperature curves
#'
#' Year temperature curves
#'
#'
#' The water temperature data sampled at the water level of Rimov
#' reservoir in Czech republic every day for the 36 years between 1979 and 2014.
#'
#'
#' @format A data frame of water temperatures in 365 days of the 36 years (rows).
#' I.e. each row consist of year temperatures in a year.
#'
#' @usage data(rimov)
#' @references
#' Mrkvička, T., Hahn, U. and Myllymäki, M. (2016)
#' A one-way ANOVA test for functional data with graphical interpretation.
#' arXiv:1612.03608 [stat.ME] (http://arxiv.org/abs/1612.03608)
#' @keywords datasets
#' @keywords curves
#' @name rimov
#' @docType data
#' @examples
#' \donttest{
#' # This is an example analysis of the water temperature data set in Mrkvicka et al. (2017).
#' data(rimov)
#' groups <- factor(c(rep(1, times=12), rep(2, times=12), rep(3, times=12)))
#'
#' nsim <- 50000 # try smaller for exploring
#'
#' res <- graph.fanova(nsim=nsim, x=rimov, groups=groups, summaryfun="means")
#' plot(res)
#' res2 <- graph.fanova(nsim=nsim, x=rimov, groups=groups, summaryfun="contrasts")
#' plot(res2)
#' }
NULL
