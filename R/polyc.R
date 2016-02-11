# based losely on Olsson, Ulf (1979), "Maximum Likelihood Estimation of the Polychoric Correlation Coefficient", Psychometrica, 44(4), 443-460.
polyc <- function(x,y,w) {
  polychoric(data.frame(x, y), weight=w)$rho[1,2]
}
