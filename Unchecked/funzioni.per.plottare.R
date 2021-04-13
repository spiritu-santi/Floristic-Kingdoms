rescale <- function(x, newrange) {
  if (nargs() > 1 && is.numeric(x) && is.numeric(newrange)) {
    if (newrange[1] > newrange[2]) {
      newmin <- newrange[2]
      newrange[2] <- newrange[1]
      newrange[1] <- newmin
    }
    xrange <- range(x)
    if (xrange[1] == xrange[2]) stop("can't rescale a constant vector!")
    mfac <- (newrange[2] - newrange[1]) / (xrange[2] - xrange[1])
    return(newrange[1] + (x - xrange[1]) * mfac)
  }
  else {
    cat("Usage: rescale(x,newrange)\n")
    cat("\twhere x is a numeric object and newrange is the min and max of the new range\n")
  }
}


hexcols2 <- function(x) {
testcol <- hex(mixcolor( seq(0,199) / 199, polarLUV(70, 50, 30),
                           polarLUV(70, 50, 120)))
  testcol2 <- hex(mixcolor( seq(0,199) / 199, polarLUV(70, 50, 300),
                            polarLUV(70, 50, 210)))
  testcols <- matrix(NA, 200, 200)
  for (i in seq_len(200)) {
    testcols[i, ] <- hex(mixcolor( (i-1) / 199, hex2RGB(testcol), hex2RGB(testcol2)))
  }
  k <- x$nobj
  plotcols <- rep(0, k)
  names(plotcols) <- as.character(1:k)
  rans <- c(max(x$points[, 1]) - min(x$points[, 1]),
    max(x$points[, 2]) - min(x$points[, 2]))
 stan <- rans / max(rans)
reord <- cbind(rescale(x$points[, 1], c(stan[1], 0)), rescale(x$points[, 2], c(stan[2], 0)))  # rescale function in r stuff
for (i in rownames(reord)) {
    plotcols[i] <- testcols[ceiling(reord[i, 1] * 199) + 1,
                            ceiling(reord[i, 2] * 199) + 1]
  }
  plotcols
}
