#' @title Run a bivariate correlation on an edsurvey.data.frame
#'
#' @description Implements weighted polyserial, weighted polychoric, weighted spearman, and weighted pearson  correlations
#'
#' @param x          a numeric (or numeric factor in case of polychoric and polyserial) vector or an object can be forced to a numeric or factor vector
#' @param y          a numeric vector (or numeric factor in case of polychoric) or an object can be forced to a numeric or factor vector
#' @param method     a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "spearman", "polychoric", or "polyserial".
#' @param weights    a numeric vector of weights
#' 
#' @details 
#' In case of polyserial, x must be the observed ordinal variable, and y the observed continuous variable.
#' the correlation methods are calculated as described in seperate documentaiton.
#   \item{Pearson}{Calculate the Pearson product moment correlation coefficient. When a variable is categorical the levels are converted to numeric 
#                  and the values are displayed in the output. For subject scales or subscales, the plausible values are used. The 95 percent confidence
#                  interval is calculated using Fisher's Z-transformation.}
#'   
#'
#'
#'
#'
#'
#'
#'
#' @return
#' Computed value of the correlation
#' 
#'
#' @seealso \ifelse{latex}{\code{cor}}{\code{\link[stats]{cor}}}
#'
#'
#' @export
weightedCorr <- function(x, y, method = c("pearson", "spearman", "polyserial", "polychoric"), weights=rep(1,length(x))) {
  x <- as.numeric(as.character(x))
  y <- as.numeric(as.character(y))
  w <- as.numric(as.character(w))
  if(!is.vector(x)) stop(paste0("The argument ",sQuote("x"), " must be a vector."))
  if(!is.vector(y)) stop(paste0("The argument ",sQuote("y"), " must be a vector."))
  if(!is.vector(w)) stop(paste0("The argument ",sQuote("w"), " must be a vector."))
  if(length(x) != length(y)) stop(paste0("The vectors ", sQuote("x"), ", ", sQuote("y"), ", and ", sQuote("w") ," must all be of the same length."))
  if(length(x) != length(w)) stop(paste0("The vectors ", sQuote("x"), ", ", sQuote("y"), ", and ", sQuote("w") ," must all be of the same length."))
  
  if (method == "polyserial") {
    value <- polys(x, y, weightVar)
  }
  
  if (method == "polychoric") {
    value <- polyc(x, y, w=weights)
  }

  if (method == "pearson" | method == "spearman") {
    value <- contCorr(x, y, w = weights, method=method)
  }
  value
}