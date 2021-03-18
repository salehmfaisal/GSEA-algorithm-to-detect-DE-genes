#setwd("C:/Users/Faisal/Desktop/Project fall")

spec <- read.csv('NRAP_JGI_Species.dat')
spec1 <- as.matrix(spec[,8:11])  # this is for SAE only
    head(spec)
   or <- levels(spec$Order)
   n.or <- length(or)


nspc <- rep(0,n.or)    # number of Species in each order
   for (i in 1:n.or){
     nspc[i] <- sum(spec$Order == or[i])
   }
nspc





or.keep <- rep(0,n.or)   # deciding which Families to display 	  
	  for (i in 1:n.or){ 
	      if(sum(spec$Order == or[i]) > 5){  
            or.keep[i] <- i 	   
		  }	
		} 

	  ors <- or[or.keep]
	  ndisp <- length(ors)
ors
ndisp



spec$X1 = log(spec[,8]+1)/mean(log(spec[,8]+1))
spec$X2 = log(spec[,9]+1)/mean(log(spec[,9]+1))
spec$X3 = log(spec[,10]+1)/mean(log(spec[,10]+1))
spec$X4 = log(spec[,11]+1)/mean(log(spec[,11]+1))


pval=rep(1,n.or)
zz <- file("CIor.txt", open = "wt")
 sink(zz)

for (i in 1 : n.or)  {
	subdata=spec[spec$Order==or[i],]
	if(nrow(subdata)>=5){

		newd1=matrix(c(subdata$X1,subdata$X2,subdata$X3,subdata$X4),ncol=4,byrow=FALSE)
		newd1

		L1=mean(newd1[,1])-mean(newd1)
		L2=mean(newd1[,2])-mean(newd1)
		L3=mean(newd1[,3])-mean(newd1)
		L4=mean(newd1[,4])-mean(newd1)

		sp.name = subdata[,7]
		years=c(0,5,23,37)
		nrow=nrow(newd1)
		newdata1=NULL

		sp.name = matrix(sp.name,nrow,1)
		Y = matrix(t(newd1),nrow*4,1)  
		month = rbind(matrix(0,nrow,1), matrix(5,nrow,1),matrix(23,nrow,1), matrix(37,nrow,1)) 
		species = rbind(sp.name,sp.name, sp.name, sp.name)
		dat = data.frame(Y = Y, month = month,species= species)
		dat

		#library("lme4")
		#lm1=lmer(Y~factor(month)+(1|species), data=dat, REML = F)

		#lm0=lmer(Y~(1|species), data=dat, REML = F)
		#an1=anova(lm0,lm1)
lm1=lm(Y~factor(month)+species, data=dat)
lm0=lm(Y~species, data=dat)
an=anova(lm1,lm0)

		pval[i]=an[[6]][2]
#fm1W <- confint(lm1, method="Wald")
# print(fm1W)

	}

}

 sink()
unlink("CI1.txt")



sum(pval < 1)
w1 = which(pval<0.05/n.or)             ## sum(pval<1)) ?
namelist1=or[w1]


spec <- read.csv('NRAP_JGI_Species.dat')
spec1 <- as.matrix(spec[,8:11])  # this is for SAE only
   head(spec)
   fams <- levels(spec$Family)
   n.fam <- length(fams)


nspc <- rep(0,n.fam)    # number of Species in each order
   for (i in 1:n.fam){
     nspc[i] <- sum(spec$Family == fams[i])
   }
nspc





fam.keep <- rep(0,n.fam)   # deciding which Families to display 	  
	  for (i in 1:n.fam){ 
	      if(sum(spec$Family == fams[i]) > 5){  
            fam.keep[i] <- i 	   
		  }	
		} 

	  fms <- fams[fam.keep]
	  ndisp <- length(fms)
fms
ndisp



spec$X1 = log(spec[,8]+1)/mean(log(spec[,8]+1))
spec$X2 = log(spec[,9]+1)/mean(log(spec[,9]+1))
spec$X3 = log(spec[,10]+1)/mean(log(spec[,10]+1))
spec$X4 = log(spec[,11]+1)/mean(log(spec[,11]+1))


pval=rep(1,n.fam)
zz <- file("CI.txt", open = "wt")
 sink(zz)

for (i in 1 : n.fam)  {
	subdata=spec[spec$Family==fams[i],]
	if(nrow(subdata)>=5){

		newd1=matrix(c(subdata$X1,subdata$X2,subdata$X3,subdata$X4),ncol=4,byrow=FALSE)
		newd1

		L1=mean(newd1[,1])-mean(newd1)
		L2=mean(newd1[,2])-mean(newd1)
		L3=mean(newd1[,3])-mean(newd1)
		L4=mean(newd1[,4])-mean(newd1)

		sp.name = subdata[,7]
		years=c(0,5,23,37)
		nrow=nrow(newd1)
		newdata1=NULL

		sp.name = matrix(sp.name,nrow,1)
		Y = matrix(t(newd1),nrow*4,1)  
		month = rbind(matrix(0,nrow,1), matrix(5,nrow,1),matrix(23,nrow,1), matrix(37,nrow,1)) 
		species = rbind(sp.name,sp.name, sp.name, sp.name)
		dat = data.frame(Y = Y, month = month,species= species)
		dat

		#library("lme4")
		#lm1=lmer(Y~factor(month)+(1|species), data=dat, REML = F)

		#lm0=lmer(Y~(1|species), data=dat, REML = F)
		#an1=anova(lm0,lm1)
lm1=lm(Y~factor(month)+species, data=dat)
lm0=lm(Y~species, data=dat)
an=anova(lm1,lm0)

		pval[i]=an[[6]][2]
#fm1W <- confint(lm1, method="Wald")
# print(fm1W)

	}

}

 sink()
unlink("CI1.txt")



sum(pval < 1)
w1 = which(pval<0.05/n.fam)             ## sum(pval<1)) ?
namelist2=fams[w1]

##Merging Two mem and Name List
nameList=NULL

for(i in 1:length(namelist1)){
nameList[[i]]=namelist1[[i]]
for (j in 1:length(namelist2)){
nameList[[i+j]]=namelist2[[j]]}
}
print(nameList)


