## Reference for the gap statistic method: Tibshirani, Walther, and Hastie. J. R. Statist. Soc. B (2001) 63, Part 2, pp. 411-423

## Inputs
## dat: data frame with an abundance column named N and either a trait column named trait or multiple trait columns for higher-dimensional tests.
## nozeros: logical. Discard species with zero abundance from the null communities?
## multiD: logical. Perform analysis using all columns of dat other than N (multitrait analysis)? If FALSE, use just the one named "trait".
## numnulls: integer. Number of null communities to test against.
## mink: integer. Minimum number of clusters to search for.
## maxk: integer. Maximum number of clusters to search for.
## nstartingpoints: integer. fed to argument "centers" of function kmeans(). Number of different random starting points for the clusters.
## weighting: integer. possible values are 0, 1, 2. Type of weight applied to within-cluster dispersal. 0 corresponds to no weighting. For options 1 and 2, see 
## 				Tibshirani et al. 2001 and Yan & Ye 2007, respectively.
## plot: logical. Plot results as a a gap curve?
## plotquant90: logical. If plotting results, plot 90th quantile of the null communities?
## verbose: logical. Print dots on console indicating which number of clusters is currently being tested?


## Outputs
## data: data frame. Rows show values corresponding to each number of clusters tested. Columns are as follows
## 		k: number of clusters tested (all integers between mink and maxk)
## 		gap: gap index = difference in log dispersal between observed community and mean of null communities
## 		Egap: mean gap index across null communities 
## 		sdgap: standard deviation of the gap index across null communities 
## 		nullquant: 95th quantile of the gap index across null communities
## 		nullquant90: 90th quantile of the gap index across null communities
##		kmaxnullquant: 95th quantile of gap index, taken only among those null communities whose max gap occurred at k clusters
##		kmaxnullquant90: 90th quantile of gap index, taken only among those null communities whose max gap occurred at k clusters
## 		logWk: log of the within-cluster dispersal returned from kmeans()
##		ElogWk: mean of the values of logWk across the null communities
## khat: integer. Number of clusters estimated for the observed community = number of clusters at which the gap index was maximal
## maxgap: real. Gap statistic = maximum value of the gap index across all number of clusters tested in the observed community
## maxnullquant: real. 95th quantile of the gap statistics of the null communities
## maxnullquant90: real. 90th quantile of the gap statistics of the null communities
## mink: integer. Minimum number of clusters tested
## maxk: integer. Maximum number of clusters tested
## z.score: real. Defined as (maxgap-mean(maxnullgap))/sd(maxnullgap)
## p.value: real. Fraction of null communities with a higher gap statistic than the observed community.
## dispersal: character indicating which weighting was used. 'D' corresponds to no weighting, 'D/2n' is the weighting in Tibshirani et al. 2001, and 
##		'D/(2n(n-1))' is the weighting in Yan & Ye 2007, where D is the within-cluster dispersal from the kmeans function.

KmeansGap=function(dat,nozeros=FALSE,multiD=FALSE,numnulls=100,mink=1,maxk=NULL,nstartingpoints=100,weighting=0,plot=FALSE,plotquant90=TRUE,verbose=TRUE){
	library(plyr)
	set.seed(0)
	
	if(class(dat)!='data.frame') stop('Input to KmeansGap must be a data frame')
	if(!multiD & !'trait'%in%names(dat)) stop('Input data frame for the 1d case must have a "trait" column')
	if(!'N'%in%names(dat)) stop('Input data frame must have a "N" column')
	
	## prepare data for analysis
	dat=dat[apply(dat,1,function(v) all(!is.na(v))),] 	# remove NAs	
	if(nozeros) dat=dat[dat$N>0,]						# remove species with zero abundance
	dat$N=dat$N/min(1,min(dat$N[dat$N>0]))				# normalize abundances such that the rarest species with positive abundance has abundance 1
	if(!multiD) dat=dat[order(dat$trait),]				# order species by trait
	
	traitdata=if(multiD) subset(dat,select=-N) else data.frame(trait=dat$trait)
	
	## If maxk not provided, set minimum of 20 or half the number of species
	if(is.null(maxk)) maxk=round(min(20,sum(dat$N>0)/2))
	
	## vector of number of clusters to test
	kvec=mink:maxk	
	
	## kms() applies kmeans() to a vector where each element is an individual in the community.
	## will first try kmeans algorithm with random initial locations for the centers of the clusters. If that fails, will place those initial centers
	## at evenly spaced locations on the trait axis.
	rep.row=function(r,a) colwise(function(x) rep(x,a))(r)
	kms=function(mat,n,k){
		if(ncol(mat)>1) com=rep.row(r=mat,a=n) else com=rep(mat$trait,n)
		if(k==1){ 
			z=kmeans(com,centers=1,iter.max=100)
		}else{
			nullcenters=apply(unique(mat),2,function(v) sample(seq(min(v),max(v),l=k)))
			z=try(kmeans(com,centers=nullcenters,iter.max=100),silent=TRUE)
			if(!is.null(attr(z,'condition'))) z=kmeans(com,centers=k,nstart=nstartingpoints,iter.max=100)
		}
		return(z)
	}
	
	## read within-cluster dispersion (sum of squared trait distances) and apply desired weighting
	Dispersion=function(dtf,weight){
		if(weight==0) return(dtf$tot.withinss)							## Unweighted sum of dist^2
		if(weight==1) return(with(dtf,sum(withinss/(2*size))))			## Tibshirani et al 2001
		if(weight==2) return(with(dtf,sum(withinss/(2*size*(size-1)))))	## Correction from Yan & Ye 2007
	}
	
	## perform the kmeans test on the observed community
	Wk=sapply(kvec,function(k){ 	## Wk = within-cluster sum of squares (ie squared distances); D in Tibshirani et al 2001
		mod=kms(mat=traitdata,n=dat$N,k)
		Dispersion(dtf=mod,weight=weighting)
	})
	
	## perform the kmeans test on the null communities
	nullWk=sapply(kvec,function(k){
		if(verbose) cat(".",if(k%%10==0 | k==maxk) paste(k,"\n"))
		
		ndat=sapply(seq(numnulls),function(run) sample(dat$N))
		return(apply(ndat,2,function(N){ 
			mod=kms(mat=traitdata,n=N,k)
			Dispersion(dtf=mod,weight=weighting)
		}))
		
	})
	
	logWk=log(Wk)										## vector length(kvec) 
	lognullWk=log(nullWk)								## matrix numnulls x length(kvec)
	ElogWk=apply(lognullWk,2,mean,na.rm=TRUE)			## vector length(kvec)
	Gapk=ElogWk-logWk									## vector length(kvec)
	nullGapk=t(ElogWk-t(lognullWk))						## matrix numnulls x length(kvec)
	
	maxgap=max(Gapk,na.rm=TRUE)							## scalar
	maxnullgap=apply(nullGapk,1,max,na.rm=TRUE)			## vector numnulls
	kmaxnullgap=apply(nullGapk,1,which.max)				## vector numnulls
	kmaxnullquant=merge(data.frame(kmax=kvec),ddply(data.frame(kmax=kmaxnullgap,max=maxnullgap),.(kmax),function(v) quantile(v$max,.95)),by='kmax',all.x=TRUE)
	kmaxnullquant90=merge(data.frame(kmax=kvec),ddply(data.frame(kmax=kmaxnullgap,max=maxnullgap),.(kmax),function(v) quantile(v$max,.9)),by='kmax',all.x=TRUE)
	khat=kvec[which.max(Gapk)]							## scalar
	z.score=(maxgap-mean(maxnullgap))/sd(maxnullgap)	## scalar
	p.value=sum(maxnullgap>maxgap)/length(maxnullgap)	## scalar
	
	mngap=apply(nullGapk,2,mean,na.rm=TRUE); 			## mean  of the null gaps --- vector length(kvec)
	sdgap=apply(nullGapk,2,sd,na.rm=TRUE); 				## stdev of the null gaps --- vector length(kvec)
	nullquant=apply(nullGapk,2,function(vec) as.numeric(quantile(vec,.95,na.rm=TRUE)))	## 95% quantile of the null gaps --- vector length(kvec)
	nullquant90=apply(nullGapk,2,function(vec) as.numeric(quantile(vec,.9,na.rm=TRUE)))	## 90% quantile of the null gaps --- vector length(kvec)
	maxnullquant=quantile(maxnullgap,.95)				## 95% quantile of the max null gaps --- scalar
	maxnullquant90=quantile(maxnullgap,.9)				## 90% quantile of the max null gaps --- scalar
	
	if(plot){		
		plot(0,t='n',xlim=c(min(kvec),max(kvec)),ylim=c(0,max(Gapk,maxnullquant)),las=1,xlab='No. clusters',ylab='Gap')
		polygon(x=c(kvec,rev(kvec)),y=c(mngap,rev(nullquant)),col='grey90',border=NA)
		lines(kvec,rep(maxnullquant,length(kvec)),col=2)
		if(plotquant90){ 
			polygon(x=c(kvec,rev(kvec)),y=c(mngap,rev(nullquant90)),col='grey80',border=NA)
			lines(kvec,rep(maxnullquant90,length(kvec)),col=2,lty=2)
		}
		points(kvec,Gapk,t='o',pch=20,lwd=2)
		polygon(x=c(kvec,rev(kvec)),y=c(mngap,rep(-maxnullquant,length(kvec))),col='white',border=NA)
		box()
	}
	
	if(weighting==0) disp='D'; if(weighting==1) disp='D/2n'; if(weighting==2) disp='D/(2n(n-1))'
	
	return(list(
		data=data.frame(
			k=kvec,
			gap=Gapk,
			Egap=mngap,
			sdgap=sdgap,
			nullquant=nullquant,
			nullquant90=nullquant90,
			kmaxnullquant=kmaxnullquant[,2],
			kmaxnullquant90=kmaxnullquant90[,2],
			logWk=logWk,
			ElogWk=ElogWk
		),
		khat=khat,
		maxgap=maxgap,
		maxnullquant=maxnullquant,
		maxnullquant90=maxnullquant90,
		mink=mink,
		maxk=maxk,
		z.score=z.score,
		p.value=p.value,
		dispersion=disp
	))
}

PlotKmeansGap=function(gap,ymin=NULL,ymax=NULL,linewidth=1,ylabb='Gap',xlabb='Number of clusters',line90=1,whiten=FALSE){
	kvec=gap$data$k; Gapk=gap$data$gap; mngap=gap$data$Egap
	nullquant=gap$data$nullquant; nullquant90=gap$data$nullquant90
	maxnullquant=gap$maxnullquant; maxnullquant90=gap$maxnullquant90
	if(is.null(ymin)) ymin=0
	if(is.null(ymax)) ymax=max(Gapk,maxnullquant)
	polymin=rep(ymin,length(kvec))
	plot(0,t='n',xlim=c(min(kvec),max(kvec)),ylim=c(ymin,ymax),las=1,xlab=xlabb,ylab=ylabb)
	polygon(x=c(kvec,rev(kvec)),y=c(polymin,rev(nullquant)),col='grey90',border=NA)
	polygon(x=c(kvec,rev(kvec)),y=c(polymin,rev(nullquant90)),col='grey80',border=NA)
	lines(kvec,rep(maxnullquant,length(kvec)),col=2)
	if(line90) lines(kvec,rep(maxnullquant90,length(kvec)),col=2,lty=2)
	points(kvec,Gapk,t='o',pch=20,lwd=linewidth)
	if(whiten) polygon(x=c(kvec,rev(kvec)),y=c(polymin,rep(-maxnullquant,length(kvec))),col='white',border=NA)
	box()
	return(NULL)
}

PlotClusters=function(trait,N,K,plot=TRUE,colors=TRUE,nstartingpoints=10000,xlab=NULL,...){
	N=N/min(1,min(N[N>0])) 
	k=kmeans(rep(trait,N),centers=K,nstart=nstartingpoints)
	dtf=unique(data.frame(trait=rep(trait,N),N=rep(N,N),cluster=k$cluster))
	dtf=dtf[order(dtf$trait),]
	dtf$cluster=match(dtf$cluster,unique(dtf$cluster))
	cols=if(colors) dtf$cluster else ifelse(dtf$cluster%%2==1,1,2)
	if(plot) with(dtf,plot(trait,N,t='h',col=cols,las=1,xlab=ifelse(is.null(xlab),'Trait',xlab),ylab='Abundance',...))
	return(dtf)	
}