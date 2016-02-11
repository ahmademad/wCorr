#' @title Run a bivariate correlation on an edsurvey.data.frame
#'
#' @description Implements weighted polyserial, weighted polychoric, weighted spearman, and weighted pearson  correlations
#'
#' @param x          a numeric (or numeric factor in case of polychoric and polyserial) vector or an object can be forced to a numeric or factor vector
#' @param y          a numeric vector (or numeric factor in case of polychoric) or an object can be forced to a numeric or factor vector
#' @param method     a character string indicating which correlation coefficient (or covariance) is to be computed. One of "Pearson" (default), "Spearman", "Polychoric", or "Polyserial".
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
weightedCorr <- function(x, y, method = c("Pearson", "Spearman", "Polyserial", "Polychoric"), weights=rep(1,length(x))) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  weights <- as.numeric(weights)
  if(!is.vector(x)) stop(paste0("The argument ",sQuote("x"), " must be a vector."))
  if(!is.vector(y)) stop(paste0("The argument ",sQuote("y"), " must be a vector."))
  if(!is.vector(weights)) stop(paste0("The argument ",sQuote("weights"), " must be a vector."))
  if(length(x) != length(y)) stop(paste0("The vectors ", sQuote("x"), ", ", sQuote("y"), ", and ", sQuote("w") ," must all be of the same length."))
  if(length(x) != length(weights)) stop(paste0("The vectors ", sQuote("x"), ", ", sQuote("y"), ", and ", sQuote("weights") ," must all be of the same length."))
  
  value <- 0
  if (method == "Polyserial") {
    value <- polys(x, y, weights)
  }
  
  if (method == "Polychoric") {
    value <- polyc(x, y, w=weights)
  }

  if (method == "Pearson" | method == "Spearman") {
    value <- contCorr(x, y, w = weights, method=method)
  }
  value
}