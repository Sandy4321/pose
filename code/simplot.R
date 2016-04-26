# simplot

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
\\clearpage
\\begin{table}
\\vspace{-.2cm}
\\footnotesize
\\caption{ 
	{\\bf  Predictive $\\boldsymbol{R^2}$ for %d observations, 
	%s design with %s covariates.}
  Reported as  \\%% worse than the Oracle 
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
				printR2line(nobs, design, support, rho,s2n,decay, fname=fname)
			if(decay!=200) cat("\\rule{0pt}{3ex}\n", file=fname, append=TRUE)
		}
	}
	cat("\\hline\\end{tabular}
\\end{center}
\\end{table}
\n", file=fname, append=TRUE)
}

printdetail <- function(fname, means, rho, s2n, decay, cpline, bg=NULL){
	for(s in c("CV.1se","CV.min","AICc","AIC","BIC")){
		l <- paste(c(s, sprintf("%s",  
			means[paste(c("gl0","gl1","gl10","mrg"),s,sep=".")])),
		collapse=" & ")
		if(!is.null(bg)) amps <- "& & &"
		else amps <- "&"
		if(s=="CV.1se") 
			cat(l,sprintf("& %s %s \\\\\n", means["snet1se"], amps), file=fname, append=TRUE)
		if(s=="CV.min"){
			if(!is.null(bg))
			 cat(l,
				sprintf("& %s & %s & & $\\mathrm{sd}(\\mathbf{\\mu})/\\sigma=%g$ \\\\\n", 
					means["snetmin"], bg["CVerr"], s2n), file=fname, append=TRUE)
			else
			 cat(l,
				sprintf("& %s %s  $\\mathrm{sd}(\\mathbf{\\mu})/\\sigma=%g$ \\\\\n", 
					 means["snetmin"], amps, s2n), file=fname, append=TRUE)
		}
		if(s=="AICc"){
		 	if(!is.null(bg))
		 		cat(l, sprintf("& & & %s &  $\\rho=%g$ \\\\\n", bg["ICerr"], rho), file=fname, append=TRUE)
		 	else
		 		cat(l, sprintf("& %s $\\rho=%g$ \\\\\n", amps, rho), file=fname, append=TRUE)
		}
		if(s=="AIC") cat(l,
			sprintf("& %s  %s \\\\\n", amps, cpline), file=fname, append=TRUE)
		if(s=="BIC") cat(l,sprintf("& %s  \\\\\n \\hline \n",amps), file=fname, append=TRUE)
	}
	invisible()
}


printMSE <- function(fname="", nobs, design, support, decay){
	preamble <- sprintf(
"\n\\clearpage
\\begin{table}\\vspace{-.5cm}
\\caption[l]{ { \\bf Prediction MSE for n=%d, %s design, 
%s covariates, and  decay  %d}.}
\\vspace{-.5cm}
\\footnotesize\\setstretch{1}
\\begin{center}
\\begin{tabular}{l*{7}{c}|r}
 & lasso & GL $\\gamma=1$ & GL $\\gamma=10$ & AL & MCP  & CVbest & ICbest  \\\\
\\cline{1-9}
", nobs, design, support, decay)
	cat(preamble, file=fname, append=TRUE)
	for(s2n in c(2,1,0.5)){
		for(rho in c(0,0.5,0.9)){
			bg <- round(colMeans(getit("bestgamma", 
				nobs, design, support, rho, s2n, decay))[c(2,4)],2)
			bgs <- sprintf("%.2f", bg)
			names(bgs) <- names(bg)
			means <- round(colMeans(getit("mse", 
				nobs=nobs, design=design, support=support, rho=rho, s2n=s2n, decay),
				na.rm=TRUE),2)
			ismin <- 1+ which(means[-1]==min(means[-1]))
			means[-ismin] <- sprintf("%.2f",means[-ismin])
			means[ismin] <- sprintf("{\\bf %s}",means[ismin])
			cpline = sprintf("\\multirow{2}{*}{$Oracle: $ %s}", means["Oracle"])
			bgmin <- which(bg==min(means[-1]))
			bgs[bgmin] <- sprintf("{\\bf %s}",bgs[bgmin])
			printdetail(fname, means, rho, s2n, decay, cpline, bgs)
		}
	}
	cat(
"\\end{tabular}
\\end{center}
\\vspace{-1cm}
\\end{table}
\n", file=fname, append=TRUE)
}



printEMSE <- function(fname="", nobs, design, support, decay){
	preamble <- sprintf(
"\n\\clearpage
\\begin{table}\\vspace{-.5cm}
\\caption[l]{ { \\bf Estimation MSE for n=%d, %s design, 
%s covariates, and  decay  %d}.}
\\vspace{-.5cm}
\\footnotesize\\setstretch{1}
\\begin{center}
\\begin{tabular}{l*{5}{c}|r}
& lasso & GL $\\gamma=1$ & GL $\\gamma=10$ & marginal AL & sparsenet MCP  & \\\\
 \\cline{1-7}
", nobs, design, support, decay)
	cat(preamble, file=fname, append=TRUE)
	for(s2n in c(2,1,0.5)){
		for(rho in c(0,0.5,0.9)){
			means <- round(colMeans(getit("estmse", 
				nobs=nobs, design=design, support=support, rho=rho, s2n=s2n, decay),
				na.rm=TRUE),3)
			ismin <- 1+ which(means[-1]==min(means[-1]))
			means[-ismin] <- sprintf("%.3f",means[-ismin])
			means[ismin] <- sprintf("{\\bf %s}",means[ismin])
			cpline = sprintf("\\multirow{2}{*}{$Oracle: $ %s}", means["Oracle"])
			printdetail(fname, means, rho, s2n, decay, cpline)
		}
	}
	cat(
"\\end{tabular}
\\end{center}
\\vspace{-1cm}
\\end{table}
\n", file=fname, append=TRUE)
}



printS <- function(fname="", nobs, design, support, decay){
	preamble <- sprintf(
"\n\\clearpage
\\begin{table}\\vspace{-.5cm}
\\caption[l]{ { \\bf Nonzero coefficients at n=%d, %s design, 
%s covariates, and  decay  %d}.}
\\vspace{-.5cm}
\\footnotesize\\setstretch{1}
\\begin{center}
\\begin{tabular}{l*{5}{c}|r}
& lasso & GL $\\gamma=1$ & GL $\\gamma=10$ & marginal AL & sparsenet MCP  & \\\\
 \\cline{1-7}
", nobs, design, support, decay)
	cat(preamble, file=fname, append=TRUE)
	for(s2n in c(2,1,0.5)){
		for(rho in c(0,0.5,0.9)){
			means <- round(colMeans(getit("s", 
				nobs=nobs, design=design, support=support, rho=rho, s2n=s2n, decay),
				na.rm=TRUE),2)
			cpline = sprintf("\\multirow{2}{*}{$Oracle: $ %s}", means["Oracle"])
			printdetail(fname, means, rho, s2n, decay, cpline)
		}
	}
	cat(
"\\end{tabular}
\\end{center}
\\vspace{-1cm}
\\end{table}
\n", file=fname, append=TRUE)
}


printsens <- function(fname="", nobs, design, support, decay){
	preamble <- sprintf(
"\n\\clearpage
\\begin{table}\\vspace{-.5cm}
\\caption[l]{ {\\it }
{ \\bf FDR $\\boldsymbol{\\mid}$ Sensitivity for n=%d, %s design, %s covariates, and  decay  %d}.}
\\vspace{-.5cm}
\\footnotesize\\setstretch{1}
\\begin{center}
\\begin{tabular}{l*{5}{c}|r}
 & lasso & GL $\\gamma=1$ & GL $\\gamma=10$ & marginal AL & sparsenet MCP  & \\\\
 \\cline{1-7}
", nobs, design, support, decay)
	cat(preamble, file=fname, append=TRUE)
	for(s2n in c(2,1,0.5)){
		for(rho in c(0,0.5,0.9)){
			fdr <- round(colMeans(getit("fdr", nobs=nobs, design=design, support=support, rho=rho, s2n=s2n, decay),na.rm=TRUE),1)
			sens <- round(colMeans(getit("sens", nobs=nobs, design=design, support=support, rho=rho, s2n=s2n, decay),na.rm=TRUE),1)
			means <- sprintf( "%.02f $\\mid$ %.02f", fdr/100, sens/100)
			names(means) <- names(fdr)
			avg <- fdr+(100-sens)
			sCp <- mean(getit("s", nobs=nobs, design=design, support=support, rho=rho,s2n=s2n, decay)[,1])
			#means[avg==min(avg[-1])] <- sprintf("{\\bf %s}",means[avg==min(avg[-1])])
			cpline = sprintf("\\multirow{2}{*}{$\\bar{s}_{Oracle}$ = %.01f}",sCp)
			printdetail(fname, means, rho, s2n, decay, cpline)
		}
	}
cat(
"\\end{tabular}
\\end{center}
\\vspace{-1cm}
\\end{table}
\n", file=fname, append=TRUE)
}


printsumline  <- function(fname="", nobs, support, s2n, decay, root=FALSE){
	if(decay=="fast") dset <- c(10,50)
	else dset <- c(100,200)

	if(s2n == 1) cat(sprintf("\\it n=%d", nobs), " & \\it ", s2n, " & ", file=fname, append=TRUE)
	else cat("& \\it ", s2n, " & ", file=fname, append=TRUE)

	r2 <- bg <- mse <- c()
	for(design in c("binary","continuous"))
		for(rho in c(0,0.5,0.9))
			for(decay in dset){
				r2 <- rbind(r2, getit("r2",nobs, design, support, rho,s2n,decay))
				bg <- rbind(bg, getit("bestgamma",nobs, design, support, rho,s2n,decay))
				mse <- rbind(mse, getit("mse",nobs, design, support, rho,s2n,decay))
		}

	mr2 <- apply(r2,2,mean)
	or2 <- mr2["Oracle"]

	if(root){
		mse <- sqrt(mse)
		bg <- sqrt(bg) }

	m <- apply(mse,2,mean)
	cp <- m["Oracle"]
	m <- m[
		c("gl0.AICc","gl0.CV.min",
		"gl1.AICc","gl1.CV.min",
		"gl10.AICc","gl10.CV.min",
		"mrg.AICc","mrg.CV.min","snetmin")]
	wtb <- round( (m-cp)/abs(cp)*100 )

	bgm <- apply(bg,2,mean)
	bgwtb <- round( (bgm[c(4,2)]-cp)/abs(cp)*100)
	wtb <- c(wtb[1:6],bgwtb,wtb[7:9])

	best <- min(wtb)
	isbest <- which( (wtb-best) < 0.01 )

	wtb <- paste(wtb)
	wtb[isbest] <- sprintf("{\\bf %s}",wtb[isbest])
	cat(wtb, sep = " & ", file=fname, append=TRUE) 
	cat(" & \\it ", sprintf("%0.2f",or2), file=fname, append=TRUE)

	if(s2n==0.5 & nobs==1000) cat(" \\\\[1ex]\n\\cline{2-2}\\rule{0pt}{3ex}", 
		file=fname, append=TRUE)
	else if(s2n==0.5 & nobs==100) cat(" \\\\[1ex]\n\\hline", 
		file=fname, append=TRUE)
	else cat(" \\\\\n", file=fname, append=TRUE)
}

printSummary <- function(fname="", root=TRUE){
	cat("%%!TEX root = pose.tex\n\n", file=fname)
	preamble <- "
\\begin{table}
\\footnotesize
\\begin{center}
\\begin{tabular}{cc|cc|cc|cc|cc|cc|c|c}
& & \\multicolumn{11}{l|}{\\bf \\% worse than oracle } & \\\\[1ex]
& \\multirow{2}{*}{$\\displaystyle\\frac{\\mathrm{sd}(\\boldsymbol{\\eta})}{\\sigma}$} 
& \\multicolumn{2}{c}{lasso} 
& \\multicolumn{2}{c}{GL $\\gamma=1$} 
& \\multicolumn{2}{c}{GL $\\gamma=10$} 
& \\multicolumn{2}{c}{GL select} 
& \\multicolumn{2}{c}{ adapt. lasso} 
& \\multicolumn{1}{c|}{~} & \\it Oracle \\\\[-0.5ex]
& 
& ~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~}
& ~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~}
& ~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~}
& ~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~}
& ~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~} 
& \\multicolumn{1}{c|}{ MCP} & $R^2$ \\\\[1ex]
\\hline
"
	cat(preamble,file=fname, append=TRUE)
	for(support in c("dense", "sparse"))
		for(decay in c("fast","slow")){
			cat(sprintf(
				"\\multicolumn{2}{l|}{\\it %s model,} &&&&&&&&&&&\\\\
\\multicolumn{2}{l|}{\\it %s decay} &&&&&&&&&&&\\\\", support, decay), 
				file=fname, append=TRUE)
			for(nobs in c(1000,100))
				for(s2n in c(2,1,1/2)){
					printsumline(fname=fname, nobs,support,s2n, decay, root=root)
			}
	}

	oracle <- c("{\\it MLE oracle on support}", "{\\it MLE oracle on true support}")
	names(oracle) <- c("dense","sparse")

	if(root) msepre <- "R"
	else msepre <- ""

	cat(sprintf("\\end{tabular}
\\end{center}
\\caption{\\label{tab:sumtables} Out-of-sample predictive %sMSE, 
reported as  \\%% worse than Oracle (corresponding $R^2$ on far right),
averaged over 1000  samples from various configurations of (\\ref{simdgp}).
The Oracle is MLE fit either on  $C_p$-optimal support for the dense model or
on the true sparse support. Each row of this table corresponds to average
performance across many data generating processes; see the supplement for more
detailed results.  Lasso (GL $\\gamma=0$), GL, and AL routines were executed
in {\\tt gamlr}.  MCP denotes results from the {\\tt sparsenet} MCP solver. GL
`select' chooses amongst $\\gamma \\in \\{0,1,10\\}$ using either AICc or CV.
The best results are bolded.}
\\end{table}
\n", msepre), file=fname, append=TRUE)
}


#printSummary("paper/sumtables.tex")



printestimsumline  <- function(fname="", nobs, support, s2n, decay, root=TRUE){
	if(decay=="fast") dset <- c(10,50)
	else dset <- c(100,200)

	if(s2n == 1) cat(sprintf("\\it n=%d", nobs), " & \\it ", s2n, " & ", file=fname, append=TRUE)
	else cat("& \\it ", s2n, " & ", file=fname, append=TRUE)

	mse <- c()
	for(design in c("binary","continuous"))
		for(rho in c(0,0.5,0.9))
			for(decay in dset){
				mse <- rbind(mse, getit("estmse",nobs, design, support, rho,s2n,decay))
		}


	if(root){
		mse <- sqrt(mse)
		}

	m <- apply(mse,2,mean)
	cp <- m["Oracle"]
	m <- m[
		c("gl0.AICc","gl0.CV.min",
		"gl1.AICc","gl1.CV.min",
		"gl10.AICc","gl10.CV.min",
		"mrg.AICc","mrg.CV.min","snetmin")]
	m <- round( m, 3)

	best <- min(m)
	isbest <- which( m==best )

	m <- paste(m)
	m[isbest] <- sprintf("{\\bf %s}",m[isbest])
	cat(m, sep = " & ", file=fname, append=TRUE) 
	cat(" & \\it ", sprintf("%0.2f",cp), file=fname, append=TRUE)

	if(s2n==0.5 & nobs==1000) cat(" \\\\[1ex]\n\\cline{2-2}\\rule{0pt}{3ex}", 
		file=fname, append=TRUE)
	else if(s2n==0.5 & nobs==100) cat(" \\\\[1ex]\n\\hline", 
		file=fname, append=TRUE)
	else cat(" \\\\\n", file=fname, append=TRUE)
}

printEstimSummary <- function(fname="", root=TRUE){
	if(root) msepre <- "R"
	else msepre <- ""
	preamble <- sprintf("
\\begin{table}[h!]
\\footnotesize
\\caption{\\label{tab:esttables}Summary of estimation %sMSE against the true coefficients,
averaged over 1000  samples from our simulation model under different designs and $\\rho$.
The Oracle is MLE fit either on  $C_p$-optimal support for the dense model or
on the true sparse support.  
The best results are bolded.}
\\begin{center}
\\vskip -.5cm
\\begin{tabular}{cc|cc|cc|cc|cc|c|c}
& & \\multicolumn{9}{l}{\\bf %sMSE} & \\\\[1ex]
& \\multirow{2}{*}{$\\displaystyle\\frac{\\mathrm{sd}(\\boldsymbol{\\eta})}{\\sigma}$} 
& \\multicolumn{2}{c}{lasso} 
& \\multicolumn{2}{c}{GL $\\gamma=1$} 
& \\multicolumn{2}{c}{GL $\\gamma=10$} 
& \\multicolumn{2}{c}{ adapt. lasso} 
& \\multicolumn{1}{c}{~} & \\\\[-0.5ex]
& 
& ~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~}
& ~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~}
& ~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~}
& ~~\\scriptsize\\it AICc & \\multicolumn{1}{c}{\\scriptsize\\it CV~~}
& \\multicolumn{1}{c}{ MCP} & \\it Oracle \\\\[1ex]
\\hline
", msepre)
	cat(preamble,file=fname, append=TRUE)
	for(support in c("dense", "sparse"))
		for(decay in c("fast","slow")){
			cat(sprintf(
				"\\multicolumn{2}{l|}{\\it %s model,} &&&&&&&&&\\\\
\\multicolumn{2}{l|}{\\it %s decay} &&&&&&&&&\\\\", support, decay), 
				file=fname, append=TRUE)
			for(nobs in c(1000,100))
				for(s2n in c(2,1,1/2)){
					printestimsumline(fname=fname, nobs,support,s2n, decay, root=root)
			}
	}

	cat(sprintf("\\end{tabular}
\\end{center}
\\end{table}
\n", msepre), file=fname, append=TRUE)
}



printsupp <- function(fname){

cat("%%!TEX root = supplemental.tex\n\n", file=fname)

print("summary")
cat("\\noindent {\\bf\\large Detailed simulation results}

\\vskip .25cm
\\noindent
Table \\ref{tab:esttables} is a summary of estimation error 
across various configurations of our simulation model, 
analogous to the predictive RMSE table in the main draft.
This is followed by detailed tabulation of prediction MSE (Tables 4-35), 
estimation MSE (36-67),  estimated model dimension (68-99), 
and sensitivity/FDR across (100-131) all 
simulation models and selection methods.
\\vskip .25cm

", file=fname, append=TRUE)

# #cat("\\subsection{Summary out-of-sample prediction results}\n\n", file=fname, append=TRUE)
# for(nobs in c(100,1000))
# 	for(support in c("dense","sparse"))
# 		for(design in c("binary","continuous")){
# 			printR2tab(nobs,design,support,fname=fname)
# 			cat("\n\n",file=fname, append=TRUE)
# 		}
printEstimSummary(fname)

print("mse")
#cat("\\subsection{Detailed out-of-sample prediction results}\n\n", file=fname, append=TRUE)
for(nobs in c(1000,100))
	for(support in c("dense","sparse"))
		for(design in c("binary","continuous"))
			for(decay in c(10,50,100,200)){		
				printMSE(fname=fname, nobs,design,support,decay)
				cat("\n\n",file=fname, append=TRUE)
		}

print("estimation")
#cat("\\subsection{Estimation error relative to true parameters}\n\n", file=fname, append=TRUE)
for(nobs in c(1000,100))
	for(support in c("dense","sparse"))
		for(design in c("binary","continuous"))
			for(decay in c(10,50,100,200)){		
				printEMSE(fname=fname, nobs,design,support,decay)
				cat("\n\n",file=fname, append=TRUE)
		}

print("support")
#cat("\\subsection{Number of estimated nonzero coefficients}\n\n", file=fname, append=TRUE)
for(nobs in c(1000,100))
	for(support in c("dense","sparse"))
		for(design in c("binary","continuous"))
			for(decay in c(10,50,100,200)){		
				printS(fname=fname, nobs,design,support,decay)
				cat("\n\n",file=fname, append=TRUE)
		}

print("sensitivity")
#cat("\\subsection{Sensitivity and false discovery relative to true parameters
# (sparse covariates) or the $\\boldsymbol{C_p}$ Oracle (dense covariates). }\n\n", file=fname, append=TRUE)
for(nobs in c(1000,100))
	for(support in c("dense","sparse"))
		for(design in c("binary","continuous"))
			for(decay in c(10,50,100,200)){		
				printsens(fname=fname, nobs,design,support,decay)
				cat("\n\n",file=fname, append=TRUE)
		}
}


printsupp("paper/simulations.tex")
