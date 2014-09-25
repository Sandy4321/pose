library(gamlr)

# load("~/project/hockey/data/nhldesign.rda")
# team <- XT
# team <- team[,colSums(team^2)!=0]
# player <- XP
# config <- XS
# goal <- goal[,c("whoscored","season","awayteam","hometeam",
# 				"period","differential","session")]
# goal$whoscored <- as.numeric(goal$whoscored=="HOME")
# goal$session <- as.numeric(goal$session=="Playoffs")
# names(goal) <- c("homegoal","season",
# 				"team.away","team.home","period",
# 				"differential","playoffs")
# rownames(team) <- rownames(config) <- rownames(player) <- rownames(goal) <- NULL
# save(goal, config, player, team, 
# 	file="~/packages/gamlr/data/hockey.rda",compress="xz")

## design 
data(hockey)
x <- cBind(config,team,player)
for(s in unique(goal$season)){
	xps <- player*(goal$season==s)
	colnames(xps) <- paste(colnames(player),s,sep=".")
	x <- cBind(x,xps)
	print(s)
}
x <- x[,colSums(x^2)!=0]
dim(x) #  69449 13007

y <- goal$homegoal
unpen <- 1:(ncol(config)+ncol(team))

####### full study
times <- list()
models <- list()
gamma = c(0,10^(-5:5))

for(v in gamma){
	cat(sprintf("\n v=%g: \n",v))
	vt <- system.time( 
		vf <- gamlr(x, y, gamma=v, tol=1e-8,
  				free=unpen, verb=0, standardize=FALSE, family="binomial") 
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
fill <- function(v) c(v,rep(tail(v,1),100-length(v)))
for(l in 1:L){
	mbic[nlambda:1,l] <- fill(BIC(models[[l]]))
	maic[nlambda:1,l] <- fill(AICc(models[[l]])) }

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
plot(v[2:9], tm[2:9], 
	xlim=c(min(v)*0.3,100), ylim=c(10,360),
	log="xy", type="l", lwd=2, 
    xlab="gamma", ylab="seconds",
    bty="n", xaxt="n",yaxt="n")
axis(1, at=c(1e-6,1e-4,1e-2,1,1e2))
axis(2,las=2,at=c(15,30,60,120,360))
lines(v[1:2],tm[1:2],lty=3, lwd=2)
points(v[1],tm[1],pch=20)
text(v[1]*.975,tm[1]*1.01,sprintf("lasso"),pos=3,font=3)
dev.off()


## selection
pdf(file="nhl_ic.pdf",width=7,height=3)
par(mfrow=c(1,2),mai=c(.9,.3,.6,0),omi=c(0,.75,0,1.2))

image(log(lambda[nlambda:1]),log(c(v[2]/2,v[2:9])),-log(maic[,1:9]),
		yaxt="n",col=grey(seq(0,1,length=255)^20),
		ylab="", xlab="", main="log AICc")
axis(2,at=log(vaxis),labels=vaxis,las=2)

image(log(lambda[nlambda:1]),log(c(v[2]/2,v[2:9])),-log(mbic[,1:9]),
		yaxt="n",col=grey(seq(0,1,length=255)^20),
		ylab="", xlab="", main="log BIC")

par(xpd=NA)
legend(x=-6.5, y=4, bty="n", fill=c(0,1), legend=c("lowest","highest"))
mtext(side=1,"log lambda",font=3,cex=1.2,outer=TRUE,line=-1.5)
mtext(side=2,"gamma",font=3,cex=1.2,outer=TRUE,line=2.5)
dev.off()

## look at estimated player [career] effects
gvec <- c(0,1,10)
for(gamma in gvec){
	fit <- gamlr(x, y, gamma=gamma, 
	  free=unpen, standardize=FALSE, family="binomial")
	plot(fit)
	
	B <- coef(fit)[-c(1,unpen+1),]
	sum(B!=0) # number of measurable effects (AICc selection)
	B[order(-B)[1:10]] # 10 biggest
	
	## grab current player effects
	now <- goal$season=="20132014"
	p1314 <- B[grep("20132014",names(B))]
	names(p1314) <- sub(".20132014","",names(p1314))
	Bnow <- B[names(B)%in%colnames(player)[
		colSums(player[now,]^2)!=0]]
	Bnow[names(p1314)] <- Bnow[names(p1314)]+p1314
	Bnow[order(-Bnow)[1:10]] # 10 biggest
	
	
	pm <- colSums(player[now,names(Bnow)]) # traditional plus minus
	ng <- colSums(abs(player[now,names(Bnow)])) # total number of goals
	# The individual effect on probability that a
	# given goal is for vs against that player's team
	p <- 1/(1+exp(-Bnow)) 
	# multiply ng*p - ng*(1-p) to get expected plus-minus
	ppm <- ng*(2*p-1)
	
	# organize the data together and print top 20
	effectg <- data.frame(player= names(Bnow),
	  beta=round(Bnow,2),ppm=round(ppm,1),pm=pm)
	rownames(effectg) <- NULL
	effectg <- effectg[order(-effectg$ppm),]
	cat("\n\ngamma = ",gamma,"number nz = ",sum(Bnow!=0),"\n")
	print(effectg[1:5,])
	if(gamma==0) effect <- effectg[1:50,]
	else effect <- cbind(effect,effectg[1:50,])
}
for(j in c(1,5,9)) effect[,j] <- as.character(effect[,j])
for(i in 1:25){
	cat(paste(as.character(effect[i,-c(2,6,10)]),collapse=" & "), "\\\\\n") }


n <- nrow(x)
nfold <- 10
rando <- sample.int(n)
chunks <- round(seq.int(0,n,length=nfold+1))
foldid <- rep.int(1:nfold,times=diff(chunks))[rando]

cva <- cv.gamlr(x=x, y=y, family="binomial", foldid=foldid, gamma=0,
	free=unpen,standardize=FALSE,verb=1)

cvb <- cv.gamlr(x=x, y=y, family="binomial", foldid=foldid, gamma=1,
	free=unpen,standardize=FALSE,verb=1)

cvc <- cv.gamlr(x=x, y=y, family="binomial", foldid=foldid, gamma=10,
	free=unpen,standardize=FALSE,verb=1)

lambda <- cva$gamlr$lambda

pdf(file="nhl_cv.pdf",width=7,height=2.5)
par(mfrow=c(1,3), 
	mai=c(.4,.4,.3,.1), 
	omi=c(.2,.2,.2,0.1))

plot(cva,xlim=c(-10,-6.75),ylim=c(1.155,1.2))
legend("topright", lty=c(2,1), col=c(1,"darkorange"),legend=c("CV","AICc"),bty="n")
mtext("gamma = 0 (lasso)", line=2, cex=.8,font=2)
abline(v=log(lambda[which.min(AICc(cva$gamlr))]),col="darkorange")
plot(cvb,xlim=c(-10,-6.75),ylim=c(1.155,1.2))
abline(v=log(lambda[which.min(AICc(cvb$gamlr))]),col="darkorange")
mtext("gamma = 1", line=2, cex=.8,font=2)
plot(cvc,xlim=c(-10,-6.75),ylim=c(1.155,1.2))
abline(v=log(lambda[which.min(AICc(cvc$gamlr))]),col="darkorange")
mtext("gamma = 10", line=2, cex=.8,font=2)
mtext(side=1, "log lambda", 
	font=3, outer=TRUE, cex=.8)
mtext(side=2, "binomial deviance", 
	font=3, outer=TRUE, cex=.8)

dev.off()
