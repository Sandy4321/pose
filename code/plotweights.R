plotfunction=function(rho, s2n, prob, decay){
png(filename=sprintf("Plots/rho%g-s2n%g-s%g-decay%g.png",rho, s2n, prob, decay))
op=par(oma=c(0,0,3,0),mfrow=c(1,3),xpd=NA)
  Lgl0=read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-Lgl0.txt",
                         rho, s2n, prob, decay),sep="|")
  for(i in 1:dim(Lgl0)[1]){
    for(j in 1:dim(Lgl0)[2]){
      if (Lgl0[i,j]==0)
        Lgl0[i,j]=NA
    }
  }
  Lgl0=colMeans(Lgl0, na.rm=T)
  
  Lgl2=read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-Lgl2.txt",
                          rho, s2n, prob, decay),sep="|")
  for(i in 1:dim(Lgl2)[1]){
    for(j in 1:dim(Lgl2)[2]){
      if (Lgl2[i,j]==0)
        Lgl2[i,j]=NA
    }
  }
  Lgl2=colMeans(Lgl2, na.rm=T)
  
  Lgl10=read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-Lgl10.txt",
                           rho, s2n, prob, decay),sep="|")
  for(i in 1:dim(Lgl10)[1]){
    for(j in 1:dim(Lgl10)[2]){
      if (Lgl10[i,j]==0)
        Lgl10[i,j]=NA
    }
  }
  Lgl10=colMeans(Lgl10, na.rm=T)

  Lmrg=read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-Lmrg.txt",
                          rho, s2n, prob, decay),sep="|")
  for(i in 1:dim(Lmrg)[1]){
    for(j in 1:dim(Lmrg)[2]){
      if (Lmrg[i,j]==0)
        Lmrg[i,j]=NA
    }
  }
  Lmrg=colMeans(Lmrg, na.rm=T)
plot(Lgl0, ylim=range(c(Lgl0, Lgl2, Lgl10, Lmrg),na.rm=T), xlab="", ylab="",type="l",col=1,main="L")
par(new=T)
plot(Lgl2, ylim=range(c(Lgl0, Lgl2, Lgl10, Lmrg),na.rm=T),axes = FALSE, xlab = "", ylab = "",type="l",col=2)
par(new=T)
plot(Lgl10, ylim=range(c(Lgl0, Lgl2, Lgl10, Lmrg),na.rm=T),axes = FALSE, xlab = "", ylab = "",type="l",col=3)
par(new=T)
plot(Lmrg, ylim=range(c(Lgl0, Lgl2, Lgl10, Lmrg),na.rm=T),axes = FALSE, xlab = "", ylab = "",type="l",col=4)
legend("topleft",legend=c("Lasso","GL2", "GL10", "AL"),lty=c(1,1,1,1),col=c(1,2,3,4),bty="n")


wmingl0=read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-wmingl0.txt",
                           rho, s2n, prob, decay),sep="|")
wmingl0=colMeans(wmingl0, na.rm=T)

wmingl2=read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-wmingl2.txt",
                           rho, s2n, prob, decay),sep="|")
wmingl2=colMeans(wmingl2, na.rm=T)

wmingl10=read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-wmingl10.txt",
                            rho, s2n, prob, decay),sep="|")
wmingl10=colMeans(wmingl10, na.rm=T)

wminmrg=read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-wminmrg.txt",
                           rho, s2n, prob, decay),sep="|")
wminmrg=colMeans(wminmrg, na.rm=T)
plot(wmingl0, ylim=range(c(wmingl0, wmingl2, wmingl10, wminmrg),na.rm=T), xlab="", ylab="",type="l",col=1,main="w min")
par(new=T)
plot(wmingl2, ylim=range(c(wmingl0, wmingl2, wmingl10, wminmrg),na.rm=T),axes = FALSE, xlab = "", ylab = "",type="l",col=2)
par(new=T)
plot(wmingl10, ylim=range(c(wmingl0, wmingl2, wmingl10, wminmrg),na.rm=T),axes = FALSE, xlab = "", ylab = "",type="l",col=3)
par(new=T)
plot(wminmrg, ylim=range(c(wmingl0, wmingl2, wmingl10, wminmrg),na.rm=T),axes = FALSE, xlab = "", ylab = "",type="l",col=4)

cpineqgl0=read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-cpineqgl0.txt",
                             rho, s2n, prob, decay),sep="|")
cpineqgl0=colMeans(cpineqgl0, na.rm=T)

cpineqgl2=read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-cpineqgl2.txt",
                             rho, s2n, prob, decay),sep="|")
cpineqgl2=colMeans(cpineqgl2, na.rm=T)

cpineqgl10=read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-cpineqgl10.txt",
                              rho, s2n, prob, decay),sep="|")
cpineqgl10=colMeans(cpineqgl10, na.rm=T)

cpineqmrg=read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-cpineqmrg.txt",
                             rho, s2n, prob, decay),sep="|")
cpineqmrg=colMeans(cpineqmrg, na.rm=T)
plot(cpineqgl0, ylim=range(c(cpineqgl0, cpineqgl2, cpineqgl10, cpineqmrg),na.rm=T), xlab="", ylab="",type="l",col=1,main="Cp")
par(new=T)
plot(cpineqgl2, ylim=range(c(cpineqgl0, cpineqgl2, cpineqgl10, cpineqmrg),na.rm=T),axes = FALSE, xlab = "", ylab = "",type="l",col=2)
par(new=T)
plot(cpineqgl10, ylim=range(c(cpineqgl0, cpineqgl2, cpineqgl10, cpineqmrg),na.rm=T),axes = FALSE, xlab = "", ylab = "",type="l",col=3)
par(new=T)
plot(cpineqmrg, ylim=range(c(cpineqgl0, cpineqgl2, cpineqgl10, cpineqmrg),na.rm=T),axes = FALSE, xlab = "", ylab = "",type="l",col=4)
mtext(sprintf("rho=%s s2n=%s s=%s decay=%s",rho, s2n, prob, decay), outer = TRUE)
dev.off()
}


