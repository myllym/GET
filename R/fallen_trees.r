#' Fallen trees
#'
#' Fallen trees
#'
#'
#' The dataset comprised the locations and heights of 232 trees, which fell during two large wind
#' gusts (1967 and 1990) in the west of France (Pontailler et al., 1997). The study area was a
#' biological reserve, which had been preserved for at least four centuries, with little human
#' influence for a long period (Guinier, 1950). Thus, the forest stand followed almost natural
#' dynamics. It was an uneven-aged beech stand with a few old oaks.
#' 
#' The data was analysed in Myllymäki et al. (2017, Supplementary material).
#'
#' @format A \code{\link{list}} of two data frames, where \code{trees} contains the locations (x and y coordinates)
#' and heights (=marks) of 232 trees in a window with polygonal boundary, and \code{window} species the polygonal
#' window (see examples).
#'
#' @usage data("fallen_trees")
#' @references
#' Guinier, P. (1950) Foresterie et protection de la nature. l’exemple de fontainebleau. Rev Forestière Fr., II, 703-717.
#'
#' Pontailler, J.-Y., Faille, A. and Lemée, G. (1997) Storms drive successional dynamics in natural forests: a case study in fontainebleau forest (france). Forest Ecol. Manag., 98, 1-15.
#' 
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#' @keywords datasets
#' @keywords spatial
#' @name fallen_trees
#' @docType data
#' @examples
#' data("fallen_trees")
#' if(require("spatstat.geom", quietly=TRUE)) {
#'   fallen_trees <- as.ppp(fallen_trees$trees, W = owin(poly=fallen_trees$window))
#'   plot(fallen_trees)
#' }
NULL
