#' @title Calculates bivariate correlation on diescrete or continious variables
#'
#' @description Implements weighted polyserial, weighted polychoric, weighted spearman, and weighted pearson correlations
#'
#' @param x          a numeric (or numeric factor in case of polychoric and polyserial) vector or an object can be forced to a numeric or factor vector
#' @param y          a numeric vector (or numeric factor in case of polychoric) or an object can be forced to a numeric or factor vector
#' @param method     a character string indicating which correlation coefficient (or covariance) is to be computed. One of "Pearson" (default), "Spearman", "Polychoric", or "Polyserial".
#' @param weights    a numeric vector of weights
#' @param ML         a boolean value indicating if full ML is to be used (polyserial and polychoric only, has no effect on Pearson or Spearman results). This substantially increases the compute time and has a very small change in the the result.
#' 
#' @details 
#' In case of polyserial, x must be the observed ordinal variable, and y the observed continuous variable. For polychoric, both must be categorical.
#' the correlation methods are calculated as described in seperate documentaiton. For Spearman the data is first ranked and then a Pearson type correlation
#' coefficient is calculated on the result. The ranking method gives averages for ties.
#'
#' The details of computation are given in the vignette.
#'
#' @return
#' A scalar that is the estimated correlation.
#'
#' @examples
#' # run a polyserial correlation
#' attach(mtcars)
#' weightedCorr(gear, x=cyl, method="polyserial")
#' # weight by MPG
#' weightedCorr(y=gear, x=cyl, method="polyserial", weights=mpg)
#'
#' # run a polychoric correlation
#' weightedCorr(gear, x=cyl, method="polychoric")
#' # weight by MPG
#' weightedCorr(y=gear, x=cyl, method="polychoric", weights=mpg)
#' detach(mtcars)
#' @seealso \ifelse{latex}{\code{cor}}{\code{\link[stats]{cor}}}
#'
#' @export
weightedCorr <- function(x, y, method = c("Pearson", "Spearman", "Polyserial", "Polychoric"), weights=rep(1,length(x)), ML=FALSE) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  weights <- as.numeric(weights)
  if(!is.vector(x)) stop(paste0("The argument ",sQuote("x"), " must be a vector."))
  if(!is.vector(y)) stop(paste0("The argument ",sQuote("y"), " must be a vector."))
  if(!is.vector(weights)) stop(paste0("The argument ",sQuote("weights"), " must be a vector."))
  if(length(x) != length(y)) stop(paste0("The vectors ", sQuote("x"), ", ", sQuote("y"), ", must be the same length."))
  if(length(x) != length(weights)) stop(paste0("The vectors ", sQuote("x"), ", ", sQuote("y"), ", and ", sQuote("weights") ," must all be of the same length."))

  value <- 0
  foundMethod <- FALSE
  #if (method == "Polyserial") {
  #  value <- polys(x, y, weights)
  #}
  method0 <- method
  method <- tolower(method)

  if(method == "polyserial") {
    if(length(unique(y)) == length(y) & length(unique(x)) < length(x)) {
      stop(paste0("Check argument definitions for ", sQuote("y"), " and ", sQuote("x") ,". The number of levels in the discrete variable ",sQuote("y")," is equal to the number of observations while the number of levels in continious variable ",  sQuote("x")," is less than the number of observations. Try transposing these two arguments."))
    }
    if(length(unique(y)) > length(unique(x))) {
      warning(paste0("Check argument definitions for ", sQuote("y"), " and ", sQuote("x") ,". The number of levels in the discrete variable ",sQuote("y")," is larger than the number of levels in continious variable ",  sQuote("x")," indicating a possible transposition of the arguments."))
    }
    if(is.factor(x)) {
      stop(paste0("The argument ", sQuote("X"), " is a factor but must be continious in the Polyserial correlation."))
    }
    value <- polys(x, y, weights, ML=ML)
    foundMethod <- TRUE
  }

  if (method == "polychoric") {
    value <- polyc(x, y, w=weights, ML=ML)
    foundMethod <- TRUE
  }

  if (method == "pearson" | method == "spearman") {
    value <- contCorr(x, y, w=weights, method=method)
    foundMethod <- TRUE
  }

  if(!foundMethod) {
    stop(paste0("Could not find method ",sQuote(method0), " see help for available methods."))
  }
  value
}