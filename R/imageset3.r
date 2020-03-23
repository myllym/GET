#' A simulated set of images
#'
#' A simulated set of images with a categorical factor
#'
#'
#' We considered a categorical factor \code{Group} obtaining the values 0, 1 or 2
#' according to the group to which the image belongs to (10 images in each of the three
#' groups). The images were simulated in the square window [-1,1]^2 from the
#' general linear model (GLM)
#' \deqn{Y(r) = \exp(-10\cdot ||r||) \cdot (1 + \mathbf{1}(g=2)) + \epsilon(r),}{Y(r) = exp(-10 ||r||) (1 + 1(g=2)) + e(r),}
#' where ||r|| denotes the Euclidean distance of the pixel to the origin, g is the group and
#' the error stems from an inhomogeneous distribution over $I$ with the normal and
#' bimodal errors in the middle and periphery of the image:
#' \deqn{\epsilon(r) = \mathbf{1}(\|r\| \leq 0.5) G(r) + \mathbf{1}(\|r\| > 0.5) \frac{1}{2}G(r)^{1/5},}{e(r) = 1(||r|| <= 0.5) G(r) + 1(||r|| > 0.5) 0.5 G(r)^{1/5},}
#' where G(r) is a Gaussian random field with the exponential correlation structure
#' with scale parameter 0.15 and standard deviation 0.2.
#' @format A list of the \code{image_set} containing the simulated images, and
#' the discrete group factor in the list component \code{Group}.
#'
#' @usage data(imageset3)
#' @references
#' Mrkvička, T., Myllymäki, M. and Narisetty, N. N. (2019) New methods for multiple testing in permutation inference for the general linear model. arXiv:1906.09004 [stat.ME]
#' @keywords datasets
#' @keywords curves
#' @name imageset3
#' @docType data
#' @seealso \code{\link{graph.fanova}}, \code{\link{frank.fanova}}
#' @examples
#' \donttest{
#' data(imageset3)
#' plot(create_image_set(imageset3$image_set), idx=c(1:5, 11:15, 21:25), max_ncols_of_plots = 5)
#' }
NULL
