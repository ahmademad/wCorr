#' @importFrom minqa bobyqa
polysFast <- function(x, M, w, ML=FALSE) {
  if(ML) {
    values = mainF(x, M, w)
    bob <- bobyqa(par=c(atanh(cor(x,M)),imapThetaFast2(theta0)),
                  fn=optFcFast, x=unlist(values[5]),  w=w, temp1=unlist(values[2]),temp2=unlist(values[3]),
                                temp3=unlist(values[4]))
    return(  tanh(bob$par[1]))
  } else {
    values = mainF(x, M, w)
     #opt <- optimize(mainF, interval = imapCorFast(cor(x,M)) + c(-3,3), x,w,temp1, temp2, temp3)
    opt <- optimize(optFcFast, interval=unlist(values[1]), unlist(values[5]), w, unlist(values[2]),unlist(values[3]), unlist(values[4]))
    return( tanh(opt$minimum) )
  }
}


