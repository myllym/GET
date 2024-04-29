#' Simulated data set
#'
#' Simulated data set mimicing stand age distributions in natural, near-natural
#' non-natural forests in the setup of Mrkvička et al. (2024).
#'
#'
#' The numbers of observations in the naturalness groups and in the dominant
#' species groups correspond to the numbers of plots in on rich mineral soils in
#' Finnish Lapland in the data set analysed in Myllymäki et al. (2023).
#' The stand age values are simulated values from the quantile regression model
#' of Mrkvička et al. (2024).
#'
#' @format A \code{data.frame} containing the dominant species (three categories:
#' Broadleaf, Conifer, Mixed), naturalness (Natural, Near-natural and Non-natural)
#' and simulated stand age.
#'
#' @usage data("naturalness")
#' @references
#' Myllymäki, M., Tuominen, S., Kuronen, M., Packalen, P. and Kangas, A. (2023)
#' The relationship between forest structure and naturalness in the Finnish national forest inventory.
#' Forestry: An International Journal of Forest Research, cpad053.
#' DOI: https://doi.org/10.1093/forestry/cpad053
#'
#' Mrkvička, T., Konstantinou, K., Kuronen, M. and Myllymäki, M. (2024)
#' Global quantile regression. arXiv:2309.04746 [stat.ME]
#' DOI: https://doi.org/10.48550/arXiv.2309.04746
#'
#' @keywords datasets
#' @keywords quantile regression
#' @name naturalness
#' @docType data
#' @seealso \code{\link{global_rq}}
#'
NULL
