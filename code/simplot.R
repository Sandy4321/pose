# simplot


pdf(file="sim_paths.pdf", width=7, height=2.5)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.3,.2), 
	omi=c(.2,.2,0,0))
for(mod in list(gl0,gl2,gl10)){
	plot(mod$g, xlab="", ylab="", select=FALSE, col=rgb(.5,.5,.75,.75))
	abline(v=log(mod$g$lambda[which.min(AICc(mod$g))]), lty=2, lwd=1.5) }
mtext(side=1, "log lambda", 
	font=3, outer=TRUE, cex=.7)
mtext(side=2, "coefficient", 
	font=3, outer=TRUE, cex=.7)
dev.off()

pdf(file="sim_cv.pdf", width=7, height=2.5)
ylim<-c(10,40)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.3,.1), 
	omi=c(.2,.2,0,0))
for(mod in list(gl0,gl2,gl10)){
	plot(mod, xlab="", ylab="", ylim=ylim) 
}
mtext(side=1, "log lambda", 
	font=3, outer=TRUE, cex=.7)
mtext(side=2, "mean square error", 
	font=3, outer=TRUE, cex=.7)
dev.off()

sum(coef(gl0$gamlr)[-1,]!=0)
sum(coef(gl2$gamlr)[-1,]!=0)
sum(coef(gl10$gamlr)[-1,]!=0)


pdf(file="sim_ic.pdf", width=7, height=2.5)
#ylim<-c(2.65,3.9)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.3,.1), 
	omi=c(.2,.2,0,0))
for(mod in list(gl0,gl2,gl10)){
	plot(log(mod$g$lam), AIC(mod$g)/n, 
		xlab="", ylab="", pch=20, col="grey75",ylim=c(2.8,3.7))
	points(log(mod$g$lam), BIC(mod$g)/n, pch=20, col=rgb(.25,.4,.7))
	if(mod$g$gamma==0) legend("topleft",h=TRUE,pch=20,
						legend=c("aic","bic"),text.col="grey25",
						col=c("grey75",rgb(.25,.4,.7)),bty="n")
	#points(log(mod$g$lam), ebf(mod,d$x)/n, pch=20, col="pink")
}
mtext(side=1, "log lambda", 
	font=3, outer=TRUE, cex=.7)
mtext(side=2, "BIC / n", 
	font=3, outer=TRUE, cex=.7)
dev.off()

