# simplot
binary <- "binary-"
getit <- function(f, rho, s2n, decay){
	segs <- c("CV.1se","CV.min","AICc","AIC","BIC")
	fname <- sprintf("results/sim-%srho%g-s2n%g-decay%g-%s.txt",
			binary, rho, s2n, decay, f)
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

printline <- function(f, rho, s2n, decay){
	if(rho==0.5 & decay==50) cat("\\it ", s2n, " & \\it ", decay, " & \\it ", rho, " & ")
	else if(rho==0.5) cat(" & \\it ", decay, " & \\it ", rho," & ")
	else cat("& & \\it ", rho, " & ")
	m <- colMeans(getit(f,rho,s2n,decay))
	cp <- m["Cp"]
	m <- m[
		c("gl0.AICc","gl0.CV.min",
		"gl1.AICc","gl1.CV.min",
		"gl10.AICc","gl10.CV.min",
		"mrg.AICc","mrg.CV.min","snetmin")]
	wtb <- round( (1-m/cp)*100 )

	best <- max(m)
	closetobest <- which( (1-m/best)<0.01 )

	isbest <- which(wtb==min(wtb))
	wtb <- paste(wtb)
	wtb[isbest] <- sprintf("{\\bf %s}",wtb[isbest])
	#wtb[setdiff(closetobest,isbest)] <- sprintf("{\\bf\\color{black!75} %s}",wtb[setdiff(closetobest,isbest)])
	#wtb[-closetobest] <- sprintf("{\\color{black!75} %s}",wtb[-closetobest])
	cat(wtb, sep = " & ") 
	cat(" & \\it ", sprintf("%0.2f",cp))
	if(rho==0.9) cat(" \\\\[1ex]\n\\cline{2-3}")
	else cat(" \\\\\n")
}

printtab <- function(){
	for(s2n in c(2,1,0.5)){
		cat("\\hline\\rule{0pt}{3ex}\n")
		for(decay in c(10,50,100,200)){
			for(rho in c(0,0.5,.9))
				printline("r2",rho,s2n,decay)
			if(decay!=200) cat("\\rule{0pt}{3ex}\n")
		}
	}
}