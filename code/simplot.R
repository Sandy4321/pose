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
		parsed <- parsed[-misfits] }
	mat <- do.call(rbind, parsed)
	if(nc==24){ colnames(mat) <- c(
		"Cp","snet1se","snetmin","scad",
		paste("mrg",segs,sep="."),
		paste("gl0",segs,sep="."),
		paste("gl2",segs,sep="."),
		paste("gl10",segs,sep="."))
	 	return(mat[,colnames(mat)!="scad"])
	 }
	if(nc==6) colnames(mat) <- c(
		"gl0","gl2","gl10","mrg","snet","scad")
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

# single run examples
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
