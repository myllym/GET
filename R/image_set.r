# Check image_set dimensions and correct r
check_image_set_dimensions <- function(image_set) {
  # Check dimensions
  obs_d <- dim(image_set$obs)
  sim_d <- dim(image_set$sim_m)
  # Check r
  r <- image_set[['r']]
  if(length(r) > 0L) {
    if(!all(unlist(lapply(r, FUN=is.vector))) || !all(unlist(lapply(r, FUN=is.numeric))) || !all(unlist(lapply(r, FUN=is.finite)))) stop("Error in image_set[[\'r\']].")
    nr <- unlist(lapply(r, FUN=length))
    if(length(r) == 2L) {
      if(!all(names(r) %in% c("x", "y"))) stop("Dimension names should be x and y")
      if(is.null(names(r))) {
        names(r) <- c("x", "y")
        image_set$r <- r
      }
      w <- diff(r$x)
      h <- diff(r$y)
      allequal <- function(x) all(abs(x - x[1]) < 1e-10*x[1])
      if(!allequal(w) || !allequal(h)) stop("Unequal gridsize detected, please specify width and height of cells.")
      if(length(r$x)!=obs_d[1] | length(r$y)!=obs_d[2]) stop("Unsuitable image_set[[\'r\']].")
    } else if(identical(sort(names(r)), c("height", "width", "x", "y"))) {
      if(length(r$width) == 1) image_set$r$width <- rep(r$width, times=length(r$x))
      else if(length(r$width)!=length(r$x)) stop("Unsuitable image_set[[\'r\']]: width should have the same length as x.")
      if(length(r$height) == 1) image_set$r$height <- rep(r$height, times=length(r$y))
      else if(length(r$height)!=length(r$x)) stop("Unsuitable image_set[[\'r\']]: height should have the same length as y.")
      if(length(r$x)!=obs_d[1]) stop("Unsuitable image_set[[\'r\']]: the length of x does not match the dimension of obs matrix.")
      if(length(r$y)!=obs_d[2]) stop("Unsuitable image_set[[\'r\']]: the length of y does not match the dimension of obs matrix.")
    } else if(identical(sort(names(r)), c("xmax", "xmin", "ymax", "ymin"))) {
      if(length(r$xmin)!=obs_d[1]) stop("Unsuitable image_set[[\'r\']]: the length of xmin does not match the dimension of obs matrix.")
      if(length(r$ymin)!=obs_d[2]) stop("Unsuitable image_set[[\'r\']]: the length of ymin does not match the dimension of obs matrix.")
      if(length(r$xmax)!=obs_d[1]) stop("Unsuitable image_set[[\'r\']]: the length of xmax does not match the dimension of obs matrix.")
      if(length(r$ymax)!=obs_d[2]) stop("Unsuitable image_set[[\'r\']]: the length of ymax does not match the dimension of obs matrix.")
    } else {
      stop("Unsuitable image_set[[\'r\']].")
    }
  }
  else {
    image_set$r <- list(x = 1:obs_d[1], y = 1:obs_d[2])
  }
  # If obs_d is 3, then set sim_m and theo to NULL
  if(length(obs_d) == 3) {
    if(!is.null(image_set$sim_m) | !is.null(image_set$theo)) {
      warning("As dim(obs) is three, sim_m (and theo) set to NULL and all the data assumed to be in obs.")
      image_set$sim_m <- NULL
      image_set$theo <- NULL
    }
  }
  image_set
}

expand_image_set_r <- function(r) {
  if(identical(sort(names(r)), c("x", "y"))) {
    xy <- expand.grid(x=r[['x']], y=r[['y']], KEEP.OUT.ATTRS = FALSE)
    xy$width <- min(diff(r[['x']]))
    xy$height <- min(diff(r[['y']]))
    xy
  } else if(identical(sort(names(r)), c("height", "width", "x", "y"))) {
    cbind(expand.grid(x=r[['x']], y=r[['y']], KEEP.OUT.ATTRS = FALSE),
          expand.grid(width=r[['width']], height=r[['height']], KEEP.OUT.ATTRS = FALSE))
  } else if(identical(sort(names(r)), c("xmax", "xmin", "ymax", "ymin"))) {
    cbind(expand.grid(xmin=r[['xmin']], ymin=r[['ymin']], KEEP.OUT.ATTRS = FALSE),
          expand.grid(xmax=r[['xmax']], ymax=r[['ymax']], KEEP.OUT.ATTRS = FALSE))
  } else {
    stop("Invalid image_set r")
  }
}

# Turn an image set object into a curve_set object vectorizing the matrices.
# The r-values are expanded to the format expected by create_curve_set.
# Should be preceeded by a call of check_image_set_dimensions.
image_set_to_curve_set <- function(image_set, ...) {
  obs_d <- dim(image_set$obs)
  sim_d <- dim(image_set$sim_m)
  if(!(length(obs_d) %in% c(2,3))) stop("Error in the dimension of image_set[['r']].")
  # Create curve_set transforming the 2d functions to vectors
  if(length(obs_d) == 3) {
    obs_v <- matrix(nrow=obs_d[1]*obs_d[2], ncol=obs_d[3])
    for(i in 1:obs_d[3]) obs_v[,i] <- as.vector(image_set$obs[,,i])
    curve_set_v <- create_curve_set(list(r=expand_image_set_r(image_set[['r']]),
                                         obs=obs_v))
  }
  else {
    n_theo <- length(image_set[['theo']])
    if(n_theo > 0L) {
      if(all(dim(theo)==obs_d) | (is.numeric(theo) & length(theo) == 1)) {
        if(length(theo) == 1) theo <- array(theo, dim=obs_d)
        theo_v <- as.vector(theo)
      }
      else stop("Unsuitable theo.")
    } else theo_v <- NULL
    if(!all(obs_d == sim_d[1:2])) stop("Something wrong with the dimensions of obs and sim_m.")
    sim_v <- matrix(nrow=sim_d[1]*sim_d[2], ncol=sim_d[3])
    for(i in 1:sim_d[3]) sim_v[,i] <- as.vector(image_set$sim_m[,,i])
    curve_set_v <- create_curve_set(list(r=expand_image_set_r(image_set[['r']]),
                                         obs=as.vector(image_set$obs),
                                         sim_m=sim_v))
    curve_set_v$theo <- theo_v
  }
  # and check them
  check_curve_set_content(curve_set_v, ...)
  curve_set_v
}

#' Create a curve set of images
#'
#' Create a curve set consisting of a set of images, given a list containing
#' the values of the 2d functions in the right form.
#' Only 2d functions in a rectangular windows are supported; the values are provided
#' in matrices (arrays). For more general 2d functions see \code{\link{create_curve_set}}.
#'
#' @param image_set A list containing elements \code{r}, \code{obs}, \code{sim_m} and \code{theo}.
#'   \code{r}, \code{sim_m} and \code{theo} are optional, \code{obs} needs to be provided always.
#'   If provided, \code{r} must be a \code{list} describing the argument values
#'   where the images have been observed (or simulated). It must consist of the following two or
#'   four components:
#'   a) "x" and "y" giving the equally spaced argument values for the x- and y-coordinates
#'   (first and second dimension of the 2d functions) where the data have been observed,
#'   b) "x", "y", "width" and "height", where the width and height give the width and height of the
#'   pixels placed at x and y, or
#'   c) "xmin", "xmax", "ymin" and "ymax" giving the corner coordinates of the pixels
#'   where the data have been observed.
#'   If not given, r is set to be a list of values from 1 to the number of first/second dimension
#'   of 2d functions in \code{obs}.
#'   \code{obs} must be either a 2d matrix (dimensions matching the lengths of r vectors)
#'   or 3d array containing the observed 2d functions (the third dimension matching the number
#'   of functions).
#'   If \code{obs} is a 3d array, then \code{sim_m} is ignored.
#'   If \code{obs} is a 2d array, then \code{sim_m} must be a 3d array containing the simulated
#'   images (2d functions) (the third dimension matching the number of functions).
#'   If included, \code{theo} corresponds to the theoretical function
#'   (e.g., under the null hypothesis) and its dimensions must either match the dimensions
#'   of 2d functions in \code{obs} or it must be a constant.
#' @param ... Do not use. (For internal use only.)
#' @return The given list as a \code{curve_set}.
#' @export
#' @examples
#' a <- create_image_set(list(obs=array(runif(4*5*6), c(4,5,6))))
#' plot(a)
#' plot(a, idx=1:6)
#'
#' a <- create_image_set(list(r=list(x=c(10,20,30,40), y=1:5*0.1),
#'                            obs=array(runif(4*5*6), c(4,5,6))))
#' plot(a)
#'
#' a <- create_image_set(list(r=list(xmin=c(1, 2, 4, 7), xmax=c(2, 4, 7, 11),
#'                                   ymin=c(1,1.1,2,2.1,3), ymax=c(1.1,2,2.1,3,3.1)),
#'                            obs=array(runif(4*5*6), c(4,5,6))))
#' plot(a)
#' plot(a, idx=1:5)
create_image_set <- function(image_set, ...) {
  image_set <- check_image_set_dimensions(image_set) # Check image_set dimensions and assign r if it does not exist
  image_set <- image_set_to_curve_set(image_set) # convert the images to functions
  image_set
}
