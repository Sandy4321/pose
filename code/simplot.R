# simplot
nobs <- 100
design <- "binary"
support <- "sparse"
rho <- 0.5
s2n <- 1
decay <- 10

getit <- function(f, nobs, design, support, rho, s2n, decay){
	segs <- c("CV.1se","CV.min","AICc","AIC","BIC")
	fname <- sprintf("results/n%d-%s-%s/sim-rho%g-s2n%g-decay%g-%s.txt",
			nobs, design, support, rho, s2n, decay, f)
	#print(fname)
	lines <- readLines(fname)
	parsed <- sapply(lines, strsplit, split="\\|",USE.NAMES=F)
	lens <- sapply(parsed,length)
	nc<- median(lens)
	misfits <- which(lens!=nc)
	if(length(misfits)>0){
		#warning(length(misfits)," overwritten line dropped")
		#print(misfits)
		parsed <- parsed[-misfits] }
	parsed <- lapply(parsed, as.numeric)
	mat <- do.call(rbind, parsed)
	if(nc==23){ colnames(mat) <- c(
		"Oracle","snet1se","snetmin",
		paste("mrg",segs,sep="."),
		paste("gl0",segs,sep="."),
		paste("gl1",segs,sep="."),
		paste("gl10",segs,sep="."))
	 }
	if(nc==4) colnames(mat) <- c("CVgam","CVerr","ICgam","ICerr")
	if(nc==6) colnames(mat) <- c("gl0","gl1","gl10","mrg","snet")
	if(nc==400) colnames(mat) <- paste(
		rep(c("mrg","gl0","gl1","gl10"),each=100), 1:100, sep=".")
	if(nc==100) colnames(mat) <- paste("seg",1:100,sep=".")
	return(mat)
}

printR2line <- function(nobs, design, support, rho, s2n, decay, fname=""){
	if(rho==0.5 & decay==50) cat("\\it ", s2n, " & \\it ", decay, " & \\it ", rho, " & ", file=fname, append=TRUE)
	else if(rho==0.5) cat(" & \\it ", decay, " & \\it ", rho," & ", file=fname, append=TRUE)
	else cat("& & \\it ", rho, " & ", file=fname, append=TRUE)
	m <- colMeans(getit("r2",nobs, design, support, rho,s2n,decay))
	cp <- m["Oracle"]
	m <- m[
		c("gl0.AICc","gl0.CV.min",
		"gl1.AICc","gl1.CV.min",
		"gl10.AICc","gl10.CV.min",
		"mrg.AICc","mrg.CV.min","snetmin")]
	wtb <- round( (cp-m)/abs(cp)*100 )

	best <- min(wtb)
	closetobest <- which( (best-m)/abs(best) <0.01 )

	isbest <- which(wtb==min(wtb))
	wtb <- paste(wtb)
	wtb[isbest] <- sprintf("{\\bf %s}",wtb[isbest])
	#wtb[setdiff(closetobest,isbest)] <- sprintf("{\\bf\\color{black!75} %s}",wtb[setdiff(closetobest,isbest)])
	#wtb[-closetobest] <- sprintf("{\\color{black!75} %s}",wtb[-closetobest])
	cat(wtb, sep = " & ", file=fname, append=TRUE) 
	cat(" & \\it ", sprintf("%0.2f",cp), file=fname, append=TRUE)
	if(rho==0.9 & decay!=200) cat(" \\\\[1ex]\n\\cline{2-3}", file=fname, append=TRUE)
	else cat(" \\\\\n", file=fname, append=TRUE)
}


printR2tab <- function(nobs, design, support, fname=""){
	otype <- c("$C_p$ optimal","true nonzero")
	names(otype) <- c("dense","sparse")

	preamble <- 
	sprintf(
"
\\begin{table}
\\vspace{-.2cm}
\\footnotesize
\\caption{ {\\bf  %d observations, %s design with %s covariates.}
  Average out-of-sample $R^2$, reported as  \\%% worse than the Oracle 
  -- MLE fit on the %s covariates -- 
  across 1000 samples.}
\\begin{center}
\\begin{tabular}{ccc|cc|cc|cc|cc|c|c}
\\hline &&&\\multicolumn{9}{|c|}{~}\\\\[-1ex]
\\multicolumn{3}{c}{~}&\\multicolumn{9}{|c|}{\\bf \\%% Worse than Oracle } &   \\\\[1ex]
& &
& \\multicolumn{2}{c}{lasso} 
& \\multicolumn{2}{c}{GL $\\gamma=1$} 
& \\multicolumn{2}{c}{GL $\\gamma=10$} 
& \\multicolumn{2}{c}{marginal AL} 
& \\multicolumn{1}{c|}{~} & \\\\[-0.5ex]
$\\mathrm{sd}(\\boldsymbol{\\eta})/\\sigma$ & {\\sf d} & $\\rho$ 
& ~~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~~}
& ~~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~~}
& ~~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~~}
& ~~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~~} 
& \\multicolumn{1}{c|}{ MCP} & Oracle $R^2$ \\\\[.5ex]
", nobs, design, support, otype[support])

	cat(preamble,file=fname, append=TRUE)
	for(s2n in c(2,1,0.5)){
		cat("\\hline\\rule{0pt}{3ex}\n", file=fname, append=TRUE)
		for(decay in c(10,50,100,200)){
			for(rho in c(0,0.5,.9))
				printR2line("r2",nobs, design, support, rho,s2n,decay, fname=fname)
			if(decay!=200) cat("\\rule{0pt}{3ex}\n", file=fname, append=TRUE)
		}
	}
	cat("\\hline\\end{tabular}
\\end{center}
\\end{table}
\n", file=fname, append=TRUE)
}


printdetail <- function(means, rho, s2n, decay, cpline, bg=NULL){
	for(s in c("CV.1se","CV.min","AICc","AIC","BIC")){
		l <- paste(c(s, sprintf("%s",  
			means[paste(c("gl0","gl1","gl10","mrg"),s,sep=".")])),
		collapse=" & ")
		if(s=="CV.1se") 
			cat(l,sprintf("& %s & & & \\\\\n", means["snet1se"]))
		if(s=="CV.min"){
			if(!is.null(bg))
			 cat(l,
				sprintf("& %s & %s & & $\\mr{sd}(\\bm{\\mu})/\\sigma=%g$ \\\\\n", 
					means["snetmin"], bg["CVerr"], s2n))
			else
			 cat(l,
				sprintf("& %s & & &  $\\mr{sd}(\\bm{\\mu})/\\sigma=%g$ \\\\\n", 
					means["snetmin"], s2n))
		}
		if(s=="AICc"){
		 	if(!is.null(bg))
		 		cat(l, sprintf("& & & %s &  $\\rho=%g$ \\\\\n", bg["ICerr"], rho))
		 	else
		 		cat(l, sprintf("& & & &  $\\rho=%g$ \\\\\n", rho))
		}
		if(s=="AIC") cat(l,
			sprintf("& & & &  %s \\\\\n", cpline))
		if(s=="BIC") cat(l,"& & & &  \\\\\n \\hline \n")
	}
	invisible()
}


printMSE <- function(nobs, design, support, decay){
	preamble <- sprintf(
"\n
\\begin{table}[p]\\vspace{-.5cm}
\\caption[l]{ Predictive MSE, for %d observations, %s design, %s covariates, and  \\textsf{d}=%d.}
\\vspace{-.5cm}
\\small\\setstretch{1}
\\begin{center}
\\begin{tabular}{l*{7}{c}|r}
 & lasso & GL $\\gamma=1$ & GL $\\gamma=10$ & AL & MCP  & ICbest & CVbest  \\\\
\\cline{1-7}
", nobs, design, support, decay)
	cat(preamble)
	for(s2n in c(2,1,0.5)){
		for(rho in c(0,0.5,0.9)){
			bg <- round(colMeans(getit("bestgamma", nobs, design, support, rho, s2n, decay))[c(2,4)],2)
			bgs <- sprintf("%.2f", bg)
			names(bgs) <- names(bg)
			means <- round(colMeans(getit("mse", 
				nobs=nobs, design=design, support=support, rho=rho, s2n=s2n, decay),na.rm=TRUE),2)
			ismin <- which(means==min(means[-1]))
			means[-ismin] <- sprintf("%.2f",means[-ismin])
			means[ismin] <- sprintf("{\\bf %s}",means[ismin])
			cpline = sprintf("\\multirow{2}{*}{$Oracle: $ = %s}", means["Oracle"])
			bgmin <- which(bg==min(means[-1]))
			bgs[bgmin] <- sprintf("{\\bf %s}",bgs[bgmin])
			printdetail(means, rho, s2n, decay, cpline, bgs)
		}
	}
	cat(
"\\end{tabular}
\\end{center}
\\vspace{-1cm}
\\end{table}
\n")
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



printsupp <- function(fname){
	doctop <- 
"\\documentclass[12pt]{article}
\\usepackage{amssymb,amsmath,setspace,anysize,times,dsfont}
\\usepackage{caption}
\\captionsetup{
  font=small,
  labelfont=normalfont,
  singlelinecheck=false,
  justification=justified
}
\\marginsize{1.1in}{.9in}{.3in}{1.4in}
\\pdfminorversion=4
\\begin{document}
\\setcounter{page}{25}
\\setcounter{equation}{20}
\\setcounter{section}{6}
\\setcounter{table}{2}
\\setcounter{figure}{3}
\\setstretch{1.3}\n\n"	

cat(doctop, file=fname)
for(nobs in c(100,1000))
	for(support in c("dense","sparse"))
		for(design in c("binary","continuous")){
			printtab(nobs,design,support,fname=fname)
			cat("\n\n",file=fname, append=TRUE)
		}
cat("\\end{document}\n",file=fname, append=TRUE)
}

printsupp("paper/simulations.tex")
