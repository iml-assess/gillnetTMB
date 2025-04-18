#' @name gillnet
#'
#' @title Gillnet data
#'
#' @description Data of an experiment with several gillnets with different mesh sizes.
#'    Data can be analysed with function \code{\link{gillnetfit}}. Data copied from TropFishR package.
#'
#' @docType data
#'
#' @format A list consiting of:
#' \itemize{
#'   \item \code{midLengths}  the midlengths of size classes,
#'   \item \code{meshSizes}  the meshsizes,
#'   \item \code{catchPerNet_mat}  a matrix with the numbers in catch of the
#'     corresponding mesh sizes (same order),
#' }
#'
#' @source Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3), 471-477.
#'
#'  Holt, S. J. 1963. A method for determining gear selectivity and its application.
#'  \emph{ICNAF Special Publication}, 5: 106-115.
#'  
#'  Mildenberge, T. K., Taylor, M. H., Wolff, M. 2024. TropfishR. R package. https://github.com/tokami/TropFishR/tree/master
#'
#' @usage data(gillnet)
#'
#' @keywords data dataset selectivity gillnet
#'
#' @examples
#' data(gillnet)
#' str(gillnet)
#' summary(gillnet)
#'
NULL