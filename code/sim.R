

n <- 50
p <- 5
rho <- 0.98
xvar <- matrix(ncol=p,nrow=p)
for(i in 1:p) for(j in i:p) xvar[i,j] <- rho^{abs(i-j)}

x <- matrix(rnorm(n*p),ncol=p)%*%chol(xvar)
x <- scale(x)*sqrt(n/(n-1))
gi <- solve(crossprod(x))
e <- rnorm(n)
e <- e-mean(e)
y <- rowSums(x)+e

H <- x%*%gi%*%t(x)
B <- gi%*%t(x)%*%y
BS <- B + gi%*%t(x)%*%y

fs <- x%*%BS
cor(fs,(diag(n)-H)%*%e)





s <- 5
g <- crossprod(x[,1:s])/n
d <- t(x[,-(1:s)])%*%(x[,1:s])/n
one <- rep(1,s)

## bickle assumption 2: this < 1/2c0
max(abs(d%*%one))/min(eigen(g)$val)

## irrepresentable
max(abs(d%*%solve(g)%*%one))






p <- %*%solve(g)%*%t(x)-diag(20)

b <- c(3,2,1,0,0)
e <- rnorm(20)
y <- x%*%b + e
y <- y-mean(y)

lm(y~x)

xxi <- solve(g[1:3,1:3])
xxxi <- x[,1:3]%*%xxi 
dj <- t(x[,4])%*%xxxi



library(Matrix)
library(gamlr)
library(snow)

smv=2
bgv=10
n <- 1e3
p <- 2*n
s2n <- 3/4

B <- 1e3
NC <- 16

## fixed 
xvar <- matrix(ncol=p,nrow=p)
for(i in 1:p) 
	for(j in i:p) 
		xvar[i,j] <- 0.9^{abs(i-j)}
C = chol(xvar)
beta <- matrix( (-1)^(1:p)*exp(-(1:p)/10) )

## random
dgp <- function(nobs){
	x <- rnorm(p*nobs)
	z <- rbinom(nobs*p,size=1,prob=.5)
	x <- matrix(x, nrow=nobs)%*%C
	x <- Matrix(x*z, sparse=TRUE)
	mu = x%*%beta
	y <- mu + rnorm(nobs,sd=sd(as.vector(mu))/s2n)
	list(x=x,y=y,mu=mu)
}

## one fit
set.seed(5807)
d <- dgp(n)

lasso <- cv.gamlr(d$x, d$y, verb=1, lambda.min.ratio=0.1)
smvar <- cv.gamlr(d$x, d$y, gamma=smv, verb=1, lambda.min.ratio=0.1)
bgvar <- cv.gamlr(d$x, d$y, gamma=bgv, verb=1, lambda.min.ratio=0.1)

pdf(file="sim_paths.pdf", width=7, height=2.5)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.3,.2), 
	omi=c(.2,.2,0,0))
for(mod in list(lasso,smvar,bgvar)){
	plot(mod$g, xlab="", ylab="", select=FALSE, col=rgb(.5,.5,.75,.75))
	abline(v=log(mod$g$lambda[which.min(BIC(mod$g))]), lty=2, lwd=1.5) }
mtext(side=1, "log lambda", 
	font=3, outer=TRUE, cex=.7)
mtext(side=2, "coefficient", 
	font=3, outer=TRUE, cex=.7)
dev.off()

pdf(file="sim_cv.pdf", width=7, height=2.5)
ylim<-c(2.65,3.9)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.3,.1), 
	omi=c(.2,.2,0,0))
for(mod in list(lasso,smvar,bgvar)){
	plot(mod, xlab="", ylab="", ylim=ylim) 
}
mtext(side=1, "log lambda", 
	font=3, outer=TRUE, cex=.7)
mtext(side=2, "mean square error", 
	font=3, outer=TRUE, cex=.7)
dev.off()


ebf <- function(fit, x){
	phat <- apply(fit$g$b!=0,2,sum)	
	ebf <- (n/2+p+1) + 0.5*fit$g$deviance
	ebf <- ebf - 0.5*phat*log(2*pi)
	H <- as.matrix(tcrossprod(t(x)))
	lDH <- apply(fit$g$b!=0,2,
		function(bi){ 
			if(!any(bi)){ return(0) }
			else if(sum(bi)==1){ return(log(H[bi,bi])) }
			else{ return(determinant(H[bi,bi])$mod) }
		})
	ebf <- ebf + 0.5*lDH
	gam <- fit$g$gamma
	lam <- fit$g$lam
	ebf <- ebf - phat*log(n*lam/2)
	if(gam > 0){ 
		b <- as.matrix(apply(x,2,sd)*abs(fit$g$b))
		r <- lam/gam
		gl <- (n*r*lam+1)*log(1+r*b)
		ebf <- ebf + colSums(gl)
	}
	return(ebf)
}

pdf(file="sim_ic.pdf", width=7, height=2.5)
#ylim<-c(2.65,3.9)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.3,.1), 
	omi=c(.2,.2,0,0))
for(mod in list(lasso,smvar,bgvar)){
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

## monte carlo
metrics <- c("AIC","BIC","CV1se","CVmin","EBF")
select <- function(fit,x){
	preds <- predict(fit$g,x,select=NULL)
	preds <- preds[,c(which.min(AIC(fit$g)),
					which.min(BIC(fit$g)),
					fit$seg.1se,
					fit$seg.min,
					which.min(ebf(fit,x)) )]
	colnames(preds) <- metrics
	as.matrix(preds)
}

eval <- function(y, f){
	f <- as.matrix(f)
	sst <- sum( (y-mean(y))^2 )
	sse <- colSums( (y-f)^2 )
	r2 <- 1-sse/sst
	return(r2)
}

cl <- makeCluster(spec=rep("localhost", NC))
clusterSetupRNG(cl)
clusterExport(cl, ls())

dofit <- function(b){

	require(Matrix)
	require(gamlr)

	d <- dgp(n*2)
	x <- d$x[1:n,]
	y <- d$y[1:n]

	## train
	lasso <- cv.gamlr(x, y)
	smvar <- cv.gamlr(x, y, gamma=smv)
	bgvar <- cv.gamlr(x, y, gamma=bgv)

	## test
	x <- d$x[n + 1:n,]
	y <- d$y[n + 1:n]

	lassopred <- select(lasso,x)
	smvarpred <- select(smvar,x)
	bgvarpred <- select(bgvar,x)

	oos <- c(eval(y, lassopred),
			eval(y, smvarpred),
			eval(y, bgvarpred))

	return(oos)
}

oos <- clusterApplyLB(cl,1:B,dofit)
stopCluster(cl)

oos <- matrix(unlist(oos),byrow=TRUE,ncol=length(metrics)*3)
colnames(oos)<- paste(
	rep(c("lasso","smvar","bgvar"),each=length(metrics)),
	metrics,sep=".")
oos <- as.data.frame(oos)
save(oos,lasso,smvar,bgvar,d,file="../results/sim.rda", compress="xz")

print(lapply(oos,mean))

crit <- c("EBF","AIC","BIC","CV1se","CVmin")
critcol <- c(rep("dodgerblue",3),rep("lawngreen",2))
pdf(file="sim_oos.pdf", width=9, height=2.5)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.2,.1), 
	omi=c(0,.2,0,0))
for(m in c("lasso","smvar","bgvar")){
	lab=paste(m,crit,sep=".")
	boxplot(oos[,lab], col=critcol,
		xaxt="n", xlab="",ylab="",ylim=c(0,.4))# ,main=m)
	axis(1,at=1:length(crit), labels=crit)
}
mtext(side=2, "R2", 
	font=3, outer=TRUE, cex=.9)
dev.off()












