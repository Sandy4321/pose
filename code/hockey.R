library(gamlr)
data(hockey)
x <- cBind(config,team,player)
for(s in unique(goal$season)){
	xps <- player*(goal$season==s)
	colnames(xps) <- paste(colnames(player),s,sep=".")
	x <- cBind(x,xps)
	print(s)
}
x <- x[,colSums(x^2)!=0]

y <- goal$homegoal
unpen <- 1:(ncol(config)+ncol(team))

## change gamma to get different levels of sparsity
fit <- gamlr(x, y, gamma=2, verb=1,
  free=unpen, standardize=FALSE, family="binomial")
plot(fit)

B <- coef(fit)[-c(1,unpen+1),]
cat(sum(B!=0), "nonzero effects\n") 
B[order(-B)[1:10]] # 10 biggest

## grab current player effects
now <- goal$season=="20132014"
p1314 <- B[grep("20132014",names(B))]
names(p1314) <- sub(".20132014","",names(p1314))
Bnow <- B[names(B)%in%colnames(player)[
	colSums(player[now,]^2)!=0]]
Bnow[names(p1314)] <- Bnow[names(p1314)]+p1314
Bnow[order(-Bnow)[1:10]] # 10 biggest


pm <- colSums(player[now,names(Bnow)]*c(-1,1)[y[now]+1]) # traditional plus minus
ng <- colSums(abs(player[now,names(Bnow)])) # total number of goals
# The individual effect on probability that a
# given goal is for vs against that player's team
p <- 1/(1+exp(-Bnow)) 
# multiply ng*p - ng*(1-p) to get expected plus-minus
ppm <- ng*(2*p-1)

# organize the data together and print top 20 by PPM
effect <- data.frame(player= names(Bnow),
  beta=round(Bnow,2),ppm=round(ppm,1),pm=pm)
rownames(effect) <- NULL
effect <- effect[order(-effect$ppm),]
print(effect[1:10,])


