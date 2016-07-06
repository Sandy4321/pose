## weight and theory plots
getit <- function(rho, s2n, decay, what, mod){
    fname <- sprintf("results/%s/sim-rho%g-s2n%g-decay%g-%s%s.txt",
            dgp, rho, s2n, decay, what, mod)
    #print(fname)
    lines <- readLines(fname)
    parsed <- sapply(lines, strsplit, split="\\|",USE.NAMES=F)
    lens <- sapply(parsed,length)
    nc<- 100
    misfits <- which(lens!=nc)
    if(length(misfits)>0){
        warning(length(misfits)," overwritten line dropped")
        #print(misfits)
        parsed <- parsed[-misfits] }
    parsed <- lapply(parsed, as.numeric)
    mat <- do.call(rbind, parsed)
    colnames(mat) <- paste("seg",1:100,sep=".")
    return(mat)
}


getavg <- function(rho, s2n, decay,what,mod){
  tab=getit(rho,s2n,decay,what,mod)
  tab[tab<=0] <- NA
  avg <- colMeans(tab, na.rm=TRUE)
  exist <- colSums(!is.na(tab))
  return(avg/exist)
}


plotweights <- function(rho,s2n,decay){
  Lgl0 <- getavg(rho,s2n,decay,"L","gl0")
  Lgl1 <- getavg(rho,s2n,decay,"L","gl1")
  Lgl10 <- getavg(rho,s2n,decay,"L","gl10")
  Lmrg <- getavg(rho,s2n,decay,"L","mrg")
  
  lam <- getavg(rho,s2n,decay,"lambda","")
  wsnormgl0 <- getavg(rho,s2n,decay,"wsnorm","gl0")/lam
  wsnormgl1 <- getavg(rho,s2n,decay,"wsnorm","gl1")/lam
  wsnormgl10 <- getavg(rho,s2n,decay,"wsnorm","gl10")/lam
  wsnormmrg <- getavg(rho,s2n,decay,"wsnorm","mrg")/lam
  
  wmingl0 <- getavg(rho,s2n,decay,"wmin","gl0")
  wmingl1 <- getavg(rho,s2n,decay,"wmin","gl1")
  wmingl10 <- getavg(rho,s2n,decay,"wmin","gl10")
  wminmrg <- getavg(rho,s2n,decay,"wmin","mrg")
  
  pdf(sprintf("weights-r%g-s%g-d%g.pdf",rho,s2n,decay), width=9, height=3)
  par(mfrow=c(1,3),mai=c(.7,.7,.1,.1))
  
  plot(wmingl0, ylim=range(c(wmingl0, wmingl1, wmingl10, wminmrg),na.rm=T), bty="n",
    xlab="", ylab="w-min",type="l",col=1, lty=2, xlim=c(1,100), cex.lab=1.5, log="y", lwd=1.5)
  lines(wmingl1, col=2, lwd=1.5)
  lines(wmingl10, col="gold", lwd=1.5)
  lines(wminmrg, col=4, lwd=1.5)
  
  plot(wsnormgl0, ylim=range(c(wsnormgl0, wsnormgl1, wsnormgl10, wsnormmrg),na.rm=T), bty="n",
    xlab="", ylab="w-norm",type="l",col=1, log="y", lty=2, xlim=c(1,100), cex.lab=1.5, lwd=1.5)
  lines(wsnormgl1, col=2, lwd=1.5)
  lines(wsnormgl10, col="gold", lwd=1.5)
  lines(wsnormmrg, col=4, lwd=1.5)
  
  plot(Lgl0, ylim=range(c(Lgl0, Lgl1, Lgl10, Lmrg),na.rm=T), bty="n",
    xlab="", ylab="L",type="l",col=1,log="y", lty=2, xlim=c(1,100), cex.lab=1.5, lwd=1.5)
  lines(Lgl1, col=2, lwd=1.5)
  lines(Lgl10, col="gold", lwd=1.5)
  lines(Lmrg, col=4, lwd=1.5)
  
  legend("topright",legend=c("Lasso","GL1", "GL10", "AL"),
      lwd=2, lty=c(2,1,1,1),col=c(1,2,"gold",4),bty="n", cex=1.3)
  
  mtext(side=1, "path segment", outer=TRUE, line=-2)
  dev.off()
}

dgp <- "n1000-binary-dense"
for(s2n in 1){
  for(decay in c(10,100)){
    for(rho in 0.9){
        plotweights(rho,s2n,decay)
    }
  }
}
