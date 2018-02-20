#' Year temperature curves
#'
#' Year temperature curves
#'
#'
#' The water temperature data sampled at the water level of Rimov
#' reservoir in Czech republic every day for the 36 years between 1979 and 2014.
#'
#'
#' @format A \code{curve_set} object with water temperatures in 365 days of the 36 years.
#' The component \code{curve_set[['r']]} is a vector of days (from 1 to 365),
#' whereas \code{curve_set[['obs']]} contains the water temperatures such that
#' each column gives year temperatures in a year.
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
#' system.time(res <- graph.fanova(nsim=2499, curve_set=rimov, groups=groups, summaryfun="means", type="erl"))
#' plot(res)
#' system.time(res2 <- graph.fanova(nsim=2499, curve_set=rimov, groups=groups, summaryfun="contrasts", type="erl"))
#' plot(res2)
#' }
NULL
