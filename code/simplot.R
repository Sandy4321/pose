# simplot
# for f in sim95*; do rename $(echo $f | grep -o "^sim[0-9]*-") "sim-" $f; done

segs <- c("CV.1se","CV.min","AICc","AIC","BIC")
getit <- function(f, rho, s2n){
	fname <- sprintf("results/sim-rho%g-s2n%g/%s.txt",
			rho, s2n, f)
	#print(fname)
	lines <- readLines(fname)
	parsed <- sapply(lines, strsplit, split="\\|",USE.NAMES=F)
	lens <- sapply(parsed,length)
	parsed <- lapply(parsed, as.numeric)
	nc<- median(lens)
	misfits <- which(lens!=nc)
	if(length(misfits)>0){
		warning(length(misfits)," overwritten line dropped")
		print(misfits)
		parsed <- parsed[-misfits] }
	mat <- do.call(rbind, parsed)
	if(nc==23){ colnames(mat) <- c(
		"Cp","snet1se","snetmin",
		paste("mrg",segs,sep="."),
		paste("gl0",segs,sep="."),
		paste("gl2",segs,sep="."),
		paste("gl10",segs,sep="."))
	 }
	if(nc==6) colnames(mat) <- c(
		"gl0","gl2","gl10","mrg","snet")
	if(nc==400) colnames(mat) <- paste(
		rep(c("mrg","gl0","gl2","gl10"),each=100), 1:100, sep=".")
	if(nc==100) colnames(mat) <- paste("seg",1:100,sep=".")
	return(mat)
}

printit <- function(means, rho, s2n, cpline){
	for(s in segs){
		l <- paste(c(s, sprintf("%s",  
			means[paste(c("gl0","gl2","gl10","mrg"),s,sep=".")])),
		collapse=" & ")
		if(s=="CV.1se") 
			cat(l,sprintf("& %s &\\\\\n", means["snet1se"]))
		if(s=="CV.min") cat(l,
			sprintf("& %s &  $\\mr{sd}(\\bm{\\mu})/\\sigma=%g$ \\\\\n", 
				means["snetmin"], s2n))
		if(s=="AICc") cat(l,
			sprintf("& & $\\rho=%g$ \\\\\n", 
				rho))
		if(s=="AIC") cat(l,
			sprintf("& & %s \\\\\n", cpline))
		if(s=="BIC") cat(l,"& & \\\\\n \\hline \n")
	}
	invisible()
}


printmse <- function(rho, s2n){
	means <- round(colMeans(getit("mse", rho=rho, s2n=s2n),na.rm=TRUE),1)
	minm <- paste(min(means[-1]))
	means <- sapply(means,as.character)
	means[means==minm] <- sprintf("{\\bf %s}",minm)
	cpline = sprintf("\\multirow{2}{*}{$C_p$ mse = %s}", means["Cp"])
	printit(means, rho, s2n, cpline)
}

for(s2n in c(2,1,0.5))
	for(rho in c(0,0.5,0.9))
	 printmse(rho=rho, s2n=s2n)


printr2 <- function(rho, s2n){
	means <- round(colMeans(getit("r2", rho=rho, s2n=s2n),na.rm=TRUE),2)
	minm <- paste(max(means[-1]))
	means <- sapply(means,as.character)
	means[means==minm] <- sprintf("{\\bf %s}",minm)
	cpline = sprintf("\\multirow{2}{*}{$C_p ~ R^2$ = %s}", means["Cp"])
	printit(means, rho, s2n, cpline)
}

for(s2n in c(2,1,0.5))
	for(rho in c(0,0.5,0.9))
	 printr2(rho=rho, s2n=s2n)

printfdrsens <- function(rho, s2n){
	fdr <- round(colMeans(getit("fdr", rho=rho, s2n=s2n),na.rm=TRUE),1)
	sens <- round(colMeans(getit("sens", rho=rho, s2n=s2n),na.rm=TRUE),1)
	means <- sprintf( "%.02f $\\mid$ %.02f", fdr/100, sens/100)
	names(means) <- names(fdr)
	avg <- fdr+(100-sens)
	sCp <- mean(getit("s",rho=rho,s2n=s2n)[,1])
	#means[avg==min(avg[-1])] <- sprintf("{\\bf %s}",means[avg==min(avg[-1])])
	cpline = sprintf("\\multirow{2}{*}{$\\bar{s}_{C_p}$ = %.01f}",sCp)
	printit(means, rho, s2n, cpline)
}

for(s2n in c(2,1,0.5))
	for(rho in c(0,0.5,0.9))
	 printfdrsens(rho=rho, s2n=s2n)

colMeans(getit("sgn", rho=.5, s2n=1),na.rm=TRUE)
colMeans(getit("times",0.5,1))

## long format


s2n <- 0.5
rho <- 0.5
fdrlong <- colMeans(getit("fdrlong",rho=rho,s2n=s2n))/100
senslong <- colMeans(getit("senslong",rho=rho,s2n=s2n))/100

pdf(file="fdr.pdf", width=7, height=2.5)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.2,.2), 
	omi=c(.2,.2,0,0))
for(mod in c("gl0","gl2","gl10")){
	irrep <- colMeans(getit(sprintf("irrep%s",mod),rho=rho,s2n=s2n))/100
	plot(irrep, type="l", ylim=c(0,1), bty="n", ylab="",xlab="",lwd=1.5)
	lines(fdrlong[grep(mod,names(fdrlong))], col=2,lwd=1.5)
	lines(senslong[grep(mod,names(senslong))], col="green",lwd=1.5)
}
mtext(side=1, "path segment", font=3, outer=TRUE, cex=.7)
legend("right",lwd=2, col=c(1,2,"green"),bty="n",
	legend=c("irrep.","FDR","sensitivity"))
dev.off()

rho <- .5
r2long <- colMeans(getit("r2long",rho=rho,s2n=s2n))
pdf(file="r2.pdf", width=7, height=2.5)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.2,.2), 
	omi=c(.2,.2,0,0))
for(mod in c("gl0","gl2","gl10")){
	cpi <- getit(sprintf("cpineq%s",mod),rho=rho,s2n=s2n)
	L <- getit(sprintf("L%s",mod),rho=rho,s2n=s2n)
	if(mod=="gl2" & s2n==1) cpi <- cpi[-c(277,278),]
	L <- colMeans(L*cpi)
	lnz <- which(L!=0)
	plot(lnz,L[lnz]/max(L),lwd=1.5,col="gold", 
		xlim=c(1,100), ylim=c(0,1), bty="n", ylab="",xlab="",type="l")
	lines(r2long[grep(mod,names(r2long))], col=4, lwd=1.5)  
	lines(colMeans(cpi),lwd=1.5)
}
mtext(side=1, "path segment", font=3, outer=TRUE, cex=.7)
legend("topright",lwd=2, col=c(1,"gold",4),bty="n",
	legend=c("min.weight ok","L/max(L)","R2"))
dev.off()





## single run examples
id <- 0
rho <- 0.5
s2n <- 1

source("code/simdata.R")
library(gamlr)

## draw the data
d <- dgp(id=id, s2n=s2n, rho=rho)

## gamma lasso
gl0 <- cv.gamlr(d$x, d$y.train)
gl2 <- cv.gamlr(d$x, d$y.train, gamma=2)
gl10 <- cv.gamlr(d$x, d$y.train, gamma=10)


pdf(file="sim_paths.pdf", width=7, height=2.5)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.3,.2), 
	omi=c(.2,.2,0,0))
for(mod in list(gl0,gl2,gl10)){
	plot(mod$g, xlab="", ylab="", select=FALSE, 
		col=rgb(.5,.5,.5,.5), xlim=c(-3.5,-.5))
	 }
mtext(side=1, "log lambda", 
	font=3, outer=TRUE, cex=.7)
mtext(side=2, "coefficient", 
	font=3, outer=TRUE, cex=.7)
dev.off()

pdf(file="sim_cv.pdf", width=7, height=2.5)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.3,.1), 
	omi=c(.2,.2,0,.2))
for(mod in list(gl0,gl2,gl10)){
	plot(mod, xlab="", ylab="", col="grey50",
		ylim=c(11,20), xlim=c(-3.5,-.5), select=FALSE) 
	lines( log(mod$g$lam), exp(AICc(gl0$g)/nrow(d$x)), lwd=2)
}
mtext(side=1, "log lambda", 
	font=3, outer=TRUE, cex=.7)
mtext(side=2, "error", 
	font=3, outer=TRUE, cex=.7)
legend("topright",bty="n", lwd=3,
	col=c("grey75","black"), 
	legend=c("CV MSE","exp(AICc/n)"))
dev.off()
