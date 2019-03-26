# Check image_set dimensions
check_image_set_dimensions <- function(image_set) {
  # Check dimensions
  obs_d <- dim(image_set$obs)
  sim_d <- dim(image_set$sim_m)
  # Check r
  r <- image_set[['r']]
  if(length(r) > 0L) {
    if(!all(unlist(lapply(r, FUN=is.vector))) || !all(unlist(lapply(r, FUN=is.numeric))) || !all(unlist(lapply(r, FUN=is.finite)))) stop("Error in image_set[[\'r\']].\n")
    nr <- unlist(lapply(r, FUN=length))
    if(length(nr) != 2 || !all(nr == c(obs_d[length(obs_d)-1], obs_d[length(obs_d)]))) stop("Unsuitable image_set[[\'r\']].\n")
  }
  else {
    image_set$r <- list(x = 1:obs_d[length(obs_d)-1], y = 1:obs_d[length(obs_d)])
  }
  # If obs_d is 3, then set sim_m and theo to NULL
  if(length(obs_d) == 3) {
    cat("dim(obs) is three, sim_m (and theo) ignored (set to NULL) and all the data assumed to be in obs.\n")
    image_set$sim_m <- NULL
    image_set$theo <- NULL
  }
  image_set
}

# Turn an image set object into a curve_set object vectorizing the matrices.
# The r-values are taken to be just from 1 to number of values of an observed function.
# Should be preceeded by a call of check_image_set_dimensions.
image_set_to_curve_set <- function(image_set, ...) {
  obs_d <- dim(image_set$obs)
  sim_d <- dim(image_set$sim_m)
  if(!(length(obs_d) %in% c(2,3))) stop("Error in the dimension of image_set[['r']].\n")
  # Create curve_set transforming the 2d functions to vectors
  if(length(obs_d) == 3) {
    obs_v <- matrix(nrow=obs_d[2]*obs_d[3], ncol=obs_d[1])
    for(i in 1:obs_d[1]) obs_v[,i] <- as.vector(image_set$obs[i,,])
    curve_set_v <- create_curve_set(list(r=1:(obs_d[2]*obs_d[3]),
                                         obs=obs_v))
  }
  else {
    n_theo <- length(image_set[['theo']])
    if(n_theo > 0L) {
      if(all(dim(theo)==obs_d) | (is.numeric(theo) & length(theo) == 1)) {
        if(length(theo) == 1) theo <- array(theo, dim=obs_d)
        theo_v <- as.vector(theo)
      }
      else stop("Unsuitable theo\n")
    } else theo_v <- NULL
    if(!all(obs_d == sim_d[2:3])) stop("Something wrong with the dimensions of obs and sim_m.\n")
    sim_v <- matrix(nrow=sim_d[2]*sim_d[3], ncol=sim_d[1])
    for(i in 1:sim_d[1]) sim_v[,i] <- as.vector(image_set$sim_m[i,,])
    curve_set_v <- create_curve_set(list(r=1:(obs_d[1]*obs_d[2]),
                                         obs=as.vector(image_set$obs),
                                         sim_m=sim_v))
    curve_set_v$theo <- theo_v
  }
  # and check them
  check_curve_set_content(curve_set_v, ...)
  curve_set_v
}

# Check the content validity of a potential image_set object.
check_image_set_content <- function(image_set) {
  image_set <- check_image_set_dimensions(image_set)
  # convert the images to functions and check the values
  image_set_to_curve_set(image_set)
  image_set
}

#' Create an image set out of a list in the right form.
#'
#' Create an image set out of a list in the right form containing the values of the 2d functions.
#' Only 2d functions in a rectangular windows are currently supported; the values are provided
#' in matrices (arrays).
#'
#' @param image_set A list containing elements \code{r}, \code{obs}, \code{sim_m} and \code{theo}.
#'   \code{r}, \code{sim_m} and \code{theo} are optional, \code{obs} needs to be provided always.
#'   \code{r} must be a list of two vectors describing the argument values where
#'   images have been observed (or simulated).
#'   r[[1]] should give the argument values for x-coordinate (first dimension of the 2d functions).
#'   r[[2]] should give the argument values for y-coordinate (second dimension of the 2d functions).
#'   If not given, r is set to be a list of values from 1 to the number of first/second dimension
#'   of 2d functions in \code{obs}.
#'   \code{obs} must be either a 2d matrix or 3d array containing the observed 2d functions.
#'   If \code{obs} is a 3d array, then \code{sim_m} is ignored.
#'   If \code{obs} is a 2d array, then \code{sim_m} must be a 3d array containing the simulated
#'   images (2d functions).
#'   If included, \code{theo} corresponds to the theoretical function
#'   (e.g., under the null hypothesis) and its dimensions must either match the dimensions
#'   of 2d functions in \code{obs} or it must be a constant.
#' @param ... Do not use. (For internal use only.)
#' @return The given list adorned with the proper class name.
#' @export
create_image_set <- function(image_set, ...) {
  check_image_set_content(image_set)
  class(image_set) <- 'image_set'
  image_set
}

#' Print method for the class 'image_set'
#' @usage \method{print}{image_set}(x, ...)
#'
#' @param x an 'image_set' object
#' @param ... Ignored.
#'
#' @method print image_set
#' @export
#' @importFrom utils str
print.image_set <- function(x, ...) {
  cat("image_set object containing:\n")
  utils::str(x)
}