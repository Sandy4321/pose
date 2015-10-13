rho=.5
s2n=1
prob=1
decay=50

getLavg <- function(mod, rho, s2n, prob, decay){
  Ltab=as.matrix(read.table(sprintf("results/sim-rho%g-s2n%g-prob%g-decay%g-L%s.txt",
                         rho, s2n, prob, decay, mod),sep="|"))
  Ltab[Ltab==0] <- NA
  Lavg <- colMeans(Ltab, na.rm=TRUE)
  Lexist <- colSums(!is.na(Ltab))
  return(Lavg/Lexist)
}


Lgl0 <- getLavg("gl0",rho,s2n,prob,decay)
Lgl2 <- getLavg("gl2",rho,s2n,prob,decay)
Lgl10 <- getLavg("gl10",rho,s2n,prob,decay)
Lmrg <- getLavg("mrg",rho,s2n,prob,decay)

plot(Lgl0, ylim=range(c(Lgl0, Lgl2, Lgl10, Lmrg),na.rm=T), 
  xlab="", ylab="",type="l",col=1,log="y", xlim=c(1,60))
lines(Lgl2, col=2)
lines(Lgl10, col=3)
lines(Lmrg, col=4)
legend("topleft",legend=c("Lasso","GL2", "GL10", "AL"),lty=c(1,1,1,1),col=c(1,2,3,4),bty="n")



#pdf(file=sprintf("plots/rho%g-s2n%g-s%g-decay%g.pdf",rho, s2n, prob, decay), width=8, height=4)

par(oma=c(0,0,3,0),mfrow=c(1,3),xpd=NA)

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
#dev.off()
}


