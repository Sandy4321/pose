library(Matrix)
library(gamlr)

data(hockey)
x <- cBind(config,onice)
y <- as.numeric(goal$who=="HOME")



times <- list()
models <- list()
gamma = c(0,10^(-5:5))

for(v in gamma){
	cat(sprintf("\n v=%g: \n",v))
	vt <- system.time( 
		vf <- gamlr(x=x, y=y, 
				free=1:ncol(config), thresh=1e-8, 
				gamma=v, lambda.min.ratio=0.01,
				family="binomial", 
				standardize=FALSE) 
		) 
	times <- c(times, list(vt))
	models <- c(models, list(vf))
	print(vt)
}

names(models) <- names(times) <- gamma

tm <- sapply(times,function(t) t["user.self"])
names(tm) <- names(times)

L <- length(tm)
v <- as.numeric(names(times))
v[1] <- 1e-6
vaxis = 10^c(-4,-2,0,2,4)

lambda <- models[[1]]$lambda
nlambda <- length(lambda)
mbic <- matrix(nrow=100,ncol=L)
maic <- matrix(nrow=100,ncol=L)
for(l in 1:L){
	mbic[nlambda:1,l] <- BIC(models[[l]]) 
	maic[nlambda:1,l] <- AIC(models[[l]]) }

## and the winner is...
data(hockey)
bm <- ceiling(which.min(mbic)/100)
am <- ceiling(which.min(maic)/100)

bmins <- sapply(models,function(m) min(BIC(m)))
nullbic <-  models[[1]]$dev[1] + log(models[[1]]$nobs)*(ncol(config)+1)
amins <- sapply(models,function(m) min(AIC(m)))
nullaic <-  models[[1]]$dev[1] + 2*(ncol(config)+1)

pdf(file="nhl_time.pdf",width=6,height=3)
par(xpd=NA,mai=c(1,.8,.4,.1),mfrow=c(1,1))
plot(v[-1], tm[-1], 
	xlim=range(v)*c(0.3,1.2), ylim=c(5,45),
	log="xy", type="l", lwd=2, 
    xlab="gamma", ylab="seconds",
    bty="n", xaxt="n",yaxt="n")
axis(1, at=vaxis)
axis(2,las=2,at=c(5,15,30,45))
lines(v[1:2],tm[1:2],lty=3, lwd=2)
points(v[1],tm[1],pch=20)
text(v[1]*.975,tm[1]*1.01,sprintf("lasso"),pos=3,font=3)
dev.off()


## selection
pdf(file="nhl_ic.pdf",width=8,height=2.75)
par(xpd=NA,mai=c(.7,.7,.4,.1),mfrow=c(1,4))

image(log(lambda[nlambda:1]),log(c(v[2]/2,v[-1])),-log(maic),
		yaxt="n",col=grey(seq(0,1,length=255)^10),
		ylab="gamma", xlab="log lambda", cex.lab=1.2)
axis(2,at=log(vaxis),labels=vaxis, las=2)
text(x=-10, y=max(log(v))+4, 
	"log AIC: ", font=3, cex=1.2)
legend(x=-9.25,y=max(log(v))+6,h=TRUE,bty="n",
	fill=c(0,1),legend=c("lowest","highest"))
points(log(lambda[which.min(AIC(models[[am]]))]),
		log(v[am]),pch=4,lwd=1.5,col=3)

par(xpd=FALSE)
plot(v[-1], amins[-1], 
	xlim=range(v)*c(0.3,1.2), ylim=c(min(amins),nullaic+100),
	log="x", type="l", lwd=2, cex.lab=1.2, 
    xlab="gamma", ylab="minimum AIC", 
    bty="n", xaxt="n")
axis(1, at=vaxis)
abline(h=nullaic,lty=2)
lines(v[1:2],amins[1:2],lty=3, lwd=2)
points(v[1],amins[1],pch=20)

par(xpd=NA)

image(log(lambda[nlambda:1]),log(c(v[2]/2,v[-1])),-log(mbic),
		yaxt="n",col=grey(seq(0,1,length=255)^10),
		ylab="gamma", xlab="log lambda", cex.lab=1.2)
axis(2,at=log(vaxis),labels=vaxis, las=2)
text(x=-10, y=max(log(v))+4, 
	"log BIC: ", font=3, cex=1.2)
legend(x=-9.25,y=max(log(v))+6,h=TRUE,bty="n",
	fill=c(0,1),legend=c("lowest","highest"))
points(log(lambda[which.min(BIC(models[[bm]]))]),
		log(v[bm]),pch=4,lwd=1.5,col=3)

par(xpd=FALSE)
plot(v[-1], bmins[-1], 
	xlim=range(v)*c(0.3,1.2), ylim=c(min(bmins),nullbic+100),
	log="x", type="l", lwd=2, 
    xlab="gamma", ylab="minimum BIC", cex.lab=1.2, 
    bty="n", xaxt="n")
axis(1, at=vaxis)
abline(h=nullbic,lty=2)
lines(v[1:2],bmins[1:2],lty=3, lwd=2)
points(v[1],bmins[1],pch=20)

dev.off()

M <- 6:8
gr=.5
poscol <- c(rgb(gr,1,gr),rgb(gr,gr,1),"grey70",rep(rgb(1,gr,gr),2))

pdf(file="nhl_paths.pdf",width=7,height=2.5)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.3,.1), 
	omi=c(.2,.2,0,0))
for(m in M){
	models[[m]]$lambda <- models[[m]]$lambda
	plot(models[[m]],xlab="",ylab="", 
		col=c(rep(0,ncol(config)),poscol[player$pos]), 
		select=FALSE,df=FALSE,
		ylim=c(-3,3), bty="n",
		main=sprintf("gamma = %g",v[m],
			font.main=1))
	abline(v=log(models[[m]]$lambda[which.min(BIC(models[[m]]))]),lwd=1.25,lty=2)
	abline(v=log(models[[m]]$lambda[which.min(AIC(models[[m]]))]),lwd=1.25,lty=2)
}
mtext(side=1, "log lambda", 
	font=3, outer=TRUE, cex=.8)
mtext(side=2, "coefficient", 
	font=3, outer=TRUE, cex=.8)
dev.off()




n <- nrow(x)
nfold <- 20
rando <- sample.int(n)
chunks <- round(seq.int(0,n,length=nfold+1))
foldid <- rep.int(1:nfold,times=diff(chunks))[rando]

# system("R CMD SHLIB ~/project/packages/gamlr/src/*.c -o gamlr.so")
# system("rm ~/project/packages/gamlr/src/*.o")
# dyn.load("gamlr.so")
# for(f in Sys.glob("~/project/packages/gamlr/R/*.R")) source(f)

cva <- cv.gamlr(x=x, y=y, family="binomial", foldid=foldid, gamma=0,
	free=1:ncol(config),standardize=FALSE,verb=TRUE)

cvb <- cv.gamlr(x=x, y=y, family="binomial", foldid=foldid, gamma=1,
	free=1:ncol(config),standardize=FALSE,verb=TRUE)

cvc <- cv.gamlr(x=x, y=y, family="binomial", foldid=foldid, gamma=10,
	free=1:ncol(config),standardize=FALSE,verb=TRUE)

lambda <- cva$gamlr$lambda

pdf(file="nhl_cv.pdf",width=7,height=2.5)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.3,.1), 
	omi=c(.2,.2,.2,0))

plot(cva)
mtext("gamma = 0 (lasso)", line=2, cex=.8,font=2)
abline(v=log(lambda[which.min(AIC(cva$gamlr))]),lty=2,col="darkorange")
abline(v=log(lambda[which.min(BIC(cva$gamlr))]),lty=2,col="darkorange")
plot(cvb)
abline(v=log(lambda[which.min(AIC(cvb$gamlr))]),lty=2,col="darkorange")
abline(v=log(lambda[which.min(BIC(cva$gamlr))]),lty=2,col="darkorange")
mtext("gamma = 1", line=2, cex=.8,font=2)
plot(cvc)
abline(v=log(lambda[which.min(AIC(cvc$gamlr))]),lty=2,col="darkorange")
abline(v=log(lambda[which.min(BIC(cva$gamlr))]),lty=2,col="darkorange")
mtext("gamma = 10", line=2, cex=.8,font=2)
mtext(side=1, "log lambda", 
	font=3, outer=TRUE, cex=.8)
mtext(side=2, "binomial deviance", 
	font=3, outer=TRUE, cex=.8)

dev.off()

#B <- coef(cva,select="min")
B <- coef(cva$gamlr,k=2)
B <- B[-(1:8),]
length(B <- B[B!=0])
print(round(B[order(-abs(B))[1:20]],3))


save.image(file="nhl.rda",compress="xz")

