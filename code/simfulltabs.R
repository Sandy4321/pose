
segs <- c("CV.1se","CV.min","AICc","AIC","BIC")
sparse <- "sparse-"

getit <- function(f, rho, s2n, decay){
	fname <- sprintf("results/sim-%srho%g-s2n%g-decay%g-%s.txt",
			sparse, rho, s2n, decay, f)
	#print(fname)
	lines <- readLines(fname)
	parsed <- sapply(lines, strsplit, split="\\|",USE.NAMES=F)
	lens <- sapply(parsed,length)
	nc<- median(lens)
	misfits <- which(lens!=nc)
	if(length(misfits)>0){
		warning(length(misfits)," overwritten line dropped")
		#print(misfits)
		parsed <- parsed[-misfits] }
	parsed <- lapply(parsed, as.numeric)
	mat <- do.call(rbind, parsed)
	if(nc==23){ colnames(mat) <- c(
		"Cp","snet1se","snetmin",
		paste("mrg",segs,sep="."),
		paste("gl0",segs,sep="."),
		paste("gl1",segs,sep="."),
		paste("gl10",segs,sep="."))
	 }
	if(nc==6) colnames(mat) <- c(
		"gl0","gl1","gl10","mrg","snet")
	if(nc==400) colnames(mat) <- paste(
		rep(c("mrg","gl0","gl1","gl10"),each=100), 1:100, sep=".")
	if(nc==100) colnames(mat) <- paste("seg",1:100,sep=".")
	return(mat)
}

printit <- function(means, rho, s2n, decay, cpline){
	for(s in segs){
		l <- paste(c(s, sprintf("%s",  
			means[paste(c("gl0","gl1","gl10","mrg"),s,sep=".")])),
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


printr2 <- function(decay){
	for(s2n in c(2,1,0.5)){
		for(rho in c(0,0.5,0.9)){
			means <- round(colMeans(getit("r2", rho=rho, s2n=s2n, decay),na.rm=TRUE),2)
			ismax <- which(means==max(means[-1]))
			means[-ismax] <- sprintf("%.2f",means[-ismax])
			means[ismax] <- sprintf("{\\bf %s}",means[ismax])
			cpline = sprintf("\\multirow{2}{*}{$C_p ~ R^2$ = %s}", means["Cp"])
			printit(means, rho, s2n, decay, cpline)
		}
	}
}


printsens <- function(decay){
	for(s2n in c(2,1,0.5)){
		for(rho in c(0,0.5,0.9)){
			fdr <- round(colMeans(getit("fdr", rho=rho, s2n=s2n, decay),na.rm=TRUE),1)
			sens <- round(colMeans(getit("sens", rho=rho, s2n=s2n, decay),na.rm=TRUE),1)
			means <- sprintf( "%.02f $\\mid$ %.02f", fdr/100, sens/100)
			names(means) <- names(fdr)
			avg <- fdr+(100-sens)
			sCp <- mean(getit("s",rho=rho,s2n=s2n, decay)[,1])
			#means[avg==min(avg[-1])] <- sprintf("{\\bf %s}",means[avg==min(avg[-1])])
			cpline = sprintf("\\multirow{2}{*}{$\\bar{s}_{C_p}$ = %.01f}",sCp)
			printit(means, rho, s2n, decay, cpline)
		}
	}
}

