#' Local brain activity at resting state
#'
#' Imaging measurements for local brain activity at resting state
#'
#'
#' The data are a small part of ABIDE fALFF data available at
#' ABIDE: http://fcon_1000.projects.nitrc.org/indi/abide/
#' fALFF: http://fcp-indi.github.io/docs/user/alff.html
#' and distributed under the CC BY-NC-SA 3.0 license,
#' https://creativecommons.org/licenses/by-nc-sa/3.0/.
#'
#' The data are fractional Amplitude of Low Frequency Fluctuations (fALFF) (Zou et al. 2008)
#' for Autism Brain Imaging Data Exchange collected resting state functional magnetic resonance
#' imaging (R-fMRI) datasets (Di Martino et al. 2013).
#' This data set in \pkg{GET} contains only a tiny part of the whole brain, namely
#' the region 9002 (the right Cerebelum Crus 1) at slice 23
#' (see Figure 2 in Mrkvicka et al., 2019) for 514 individuals with the autism spectrum
#' disorder (ASD) and 557 typical controls (TC) as specified in the given Group variable.
#' Further the sex and age of each subject is given.
#'
#' @format A list of the \code{curve_set} containing the data,
#' coordinates (x,y) where the data have been observed (third dimension is 23),
#' the discrete factor \code{Group} (1=Autism; 2=Control),
#' the discrete factor \code{Sex} (1=Male; 2=Female),
#' and the continuous factor \code{Age}.
#'
#' @usage data("abide_9002_23")
#' @references
#' Di Martino, A., Yan, C., Li, Q., Denio, E., Castellanos, F., Alaerts, K., Anderson, J., Assaf, M., Bookheimer, S., Dapretto, M., et al. (2013) The autism brain imaging data exchange: towards a large-scale evaluation of the intrinsic brain architecture in autism. Molecular psychiatry.
#'
#' Tzourio-Mazoyer, N., Landeau, B., Papathanassiou, D., Crivello, F., Etard, O., Delcroix, N., Mazoyer, B., and Joliot, M. (2002), Automated anatomical labeling of activations in SPM using a macroscopic anatomical parcellation of the MNI MRI single-subject brain. Neuroimage, 15, 273-289.
#'
#' Zou, Q.-H., Zhu, C.-Z., Yang, Y., Zuo, X.-N., Long, X.-Y., Cao, Q.-J., Wang, Y.-F., and Zang, Y.-F. (2008), An improved approach to detection of amplitude of low-frequency fluctuation (ALFF) for resting-state fMRI: fractional ALFF. Journal of neuroscience methods, 172, 137-141.
#'
#' Mrkvička, T., Myllymäki, M., Kuronen, M. and Narisetty, N. N. (2020) New methods for multiple testing in permutation inference for the general linear model. arXiv:1906.09004 [stat.ME]
#' @keywords datasets
#' @name abide_9002_23
#' @docType data
NULL
