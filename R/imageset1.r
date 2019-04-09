#' A simulated set of images
#'
#' A simulated set of images with two simulated covariates
#'
#'
#' We considered a categorical factor \code{Group} obtaining the values 1 or 2
#' according to the group to which the image belongs to (10 images in the first group,
#' 10 images in the second), and a continous factor z that was generated from the uniform
#' distribution on (0,1). The images were simulated in the square window [-1,1]^2 from the
#' general linear model (GLM)
#' \deqn{Y(r) = \exp(-10\cdot ||r||) \cdot (g + z) + \epsilon(r),}{Y(r) = exp(-10 ||r||) (g + z) + e(r),}
#' where ||r|| denotes the Euclidean distance of the pixel to the origin and
#' the error stems from an inhomogeneous distribution over $I$ with the normal and
#' bimodal errors in the middle and periphery of the image:
#' \deqn{\epsilon(r) = \mathbb I (\|r\| \leq 0.5) G(r) + \mathbb I (\|r\| > 0.5) \frac{1}{2}G(r)^{1/5},}{e(r) = 1 (||r|| <= 0.5) G(r) + 1(||r|| > 0.5) 0.5 G(r)^{1/5},}
#' where G(r) is a Gaussian random field with the exponential correlation structure
#' with scale parameter 0.15 and standard deviation 0.2.
#' @format A list of the \code{image_set} containing the simulated images,
#' the discrete group factor in the list component \code{Group},
#' and the continuous factor z in the list component \code{z}.
#'
#' @usage data(imageset1)
#' @references
#' Mrkvička, T., Myllymäki, M. and Narisetty, N. N. (2019) New methods for multiple testing in permutation inference for the general linear model.
#' @keywords datasets
#' @keywords curves
#' @name imageset1
#' @docType data
#' @seealso \code{\link{graph.fanova2d}}, \code{\link{graph.fglm2d}}
#' @examples
#' data(imageset1)
#' par(mfrow=c(2,5))
#' plot(imageset1$image_set, idx=c(1:5, 11:15))
NULL
