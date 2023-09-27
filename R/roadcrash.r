#' Road crashes
#'
#' Road crashes
#'
#'
#' Mrkvička et al. (2023) worked with the database of road crashes
#' reported to the Police in the Czech Republic from 1 January 2016 to 31
#' December 2020. The data available here is a subpattern of this full data set,
#' included here with the permission of the Police in the Czech Republic.
#' The full data is published as open data, see
#' \url{https://policyvisuals.eu/traffic-accidents-data-in-the-czech-republic/}.
#' The subpattern 7700 crashes lying on a linear network with 269
#' vertices and 354 lines. Further average traffic volume (number of vehicles
#' per 24 hours), forest density and building density in the cell are available
#' in the region of the linear network.
#'
#'
#' @format A \code{list} with
#' \itemize{
#' \item x = x-coordinates of road accidents
#' \item y = y-coordinates of road accidents
#' \item xrange = x coordinate limits of enclosing box (-774936.86,-727048.86)
#' \item yrange = y coordinate limits of enclosing box (-1201599.83,-1125679.83)
#' \item Vertices.x = x-coordinates of vertices of the linear network
#' \item Vertices.y = y-coordinates of vertices of the linear network
#' \item Edges = a 2 column matrix giving the ID (index) of the origin and destination vectices
#' (in Vertices.x and Vertices.y)
#' \item Traffic = matrix of traffic volume
#' \item ForestDensity = matrix of forest density
#' \item BuildingDensity = matrix of building density
#' }
#'
#' @usage data("roadcrash")
#' @references
#' Mrkvička, T., Kraft, S., Blažek, V. and Myllymäki, M. (2023) Hotspots detection on a linear network with presence of covariates: a case study on road crash data.
#' @keywords datasets
#' @name roadcrash
#' @docType data
#' @examples
#' if(require("spatstat.geom", quietly = TRUE) & require("spatstat.linnet", quietly = TRUE)) {
#'   data("roadcrash")
#'   win <- owin(xrange = roadcrash$xrange,
#'               yrange = roadcrash$yrange)
#'   X <- ppp(x = roadcrash$x, y = roadcrash$y, window = win)
#'   Vertices.pp <- ppp(x = roadcrash$Vertices.x,
#'                      y = roadcrash$Vertices.y,
#'                      window=win)
#'   L <- linnet(vertices=Vertices.pp,
#'               edges = roadcrash$Edges)
#'   PP <- lpp(X, L)
#'   z1 <- im(roadcrash$Traffic,
#'            xrange = roadcrash$xrange,
#'            yrange = roadcrash$yrange)
#'   z2 <- im(roadcrash$ForestDensity,
#'            xrange = roadcrash$xrange,
#'            yrange = roadcrash$yrange)
#'   z3 <- im(roadcrash$BuildingDensity,
#'            xrange = roadcrash$xrange,
#'            yrange = roadcrash$yrange)
#' }
NULL
