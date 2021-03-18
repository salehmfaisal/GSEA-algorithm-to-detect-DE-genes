library(EnrichmentBrowser)
install.packages("statmod")
library(statmod)
library(edgeR)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EDASeq")
library(EDASeq)

##Making Data Set
xx=data[,8:11]
write.table(xx,"exprs.tab",row.names=F, sep="\t")
x=seq(1,length(data$Species), by=1)
write.csv(x,"rowData1.tab",row.names=F)
data.dir =  "C:/Users/afais/Documents"
  exprs.file <- file.path(data.dir, "exprs.tab")
  cdat.file <- file.path(data.dir, "colData.tab")
  rdat.file <- file.path(data.dir, "rowData1.tab")   # need to update this file to have ENTREZ id's instead 
  se <- readSE(exprs.file, cdat.file, rdat.file)
gmt.file <- file.path(data.dir, "microbialset.gmt")
  go.filt.ms <- getGenesets(gmt.file)
 length(go.filt.ms)
  assay(se)[1:4,1:4]
  ourSE = se


#Differentially expression Analysis
  ourSE1 <- deAna(ourSE, de.method="edgeR")
  rowData(ourSE1, use.names=TRUE)


##Set based enrichedment analysis  
sbea.res <- sbea(method="camera", se= ourSE1, gs= go.filt.ms, perm=0, padj.method="BH", alpha=0.1) 
# "perm" is irrelevant for CAMERA 
gsRanking(sbea.res)   
  
sbea.res1 <- sbea(method="gsea", se= ourSE1, gs= go.filt.ms, perm=100, padj.method="BH", alpha=0.1)  
gsRanking(sbea.res1)   
list1 = gsRanking(sbea.res1)$GENE.SET 
list2=gsRanking(sbea.res)$GENE.SET
intersect(list1, list2)



##Check With previous result
##Making the nameList for Family and order
spec <- read.csv('NRAP_JGI_Species.dat')
spec1 <- as.matrix(spec[,8:11])  # this is for SAE only
    head(spec)
   or <- levels(spec$Order)
   n.or <- length(or)
or.keep <- rep(0,n.or)   # deciding which Orders to display 	  
	  for (i in 1:n.or){ 
	      if(sum(spec$Order == or[i]) > 5){  
            or.keep[i] <- i 	   
		  }	
		} 

	  ors <- or[or.keep]
	  ndisp <- length(ors)

spec$X1 = log(spec[,8]+1)/mean(log(spec[,8]+1))
spec$X2 = log(spec[,9]+1)/mean(log(spec[,9]+1))
spec$X3 = log(spec[,10]+1)/mean(log(spec[,10]+1))
spec$X4 = log(spec[,11]+1)/mean(log(spec[,11]+1))


pval=rep(1,n.or)
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
lm1=lm(Y~factor(month)+species, data=dat)
lm0=lm(Y~species, data=dat)
an=anova(lm1,lm0)

		pval[i]=an[[6]][2]
	}

}

sum(pval < 1)
w1 = which(pval<0.05/n.or)             ## sum(pval<1)) ?
namelist1=or[w1]


spec <- read.csv('NRAP_JGI_Species.dat')
   fams <- levels(spec$Family)
   n.fam <- length(fams)


fam.keep <- rep(0,n.fam)   # deciding which Families to display 	  
	  for (i in 1:n.fam){ 
	      if(sum(spec$Family == fams[i]) > 5){  
            fam.keep[i] <- i 	   
		  }	
		} 

	  fms <- fams[fam.keep]
spec$X1 = log(spec[,8]+1)/mean(log(spec[,8]+1))
spec$X2 = log(spec[,9]+1)/mean(log(spec[,9]+1))
spec$X3 = log(spec[,10]+1)/mean(log(spec[,10]+1))
spec$X4 = log(spec[,11]+1)/mean(log(spec[,11]+1))


pval=rep(1,n.fam)

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
lm1=lm(Y~factor(month)+species, data=dat)
lm0=lm(Y~species, data=dat)
an=anova(lm1,lm0)

		pval[i]=an[[6]][2]
	}

}



sum(pval < 1)
w2 = which(pval<0.05/n.fam)             ## sum(pval<1)) ?
namelist2=fams[w2]

##Merging Two mem and Name List
nameList=NULL

for(i in 1:length(namelist1)){
nameList[[i]]=namelist1[[i]]
for (j in 1:length(namelist2)){
nameList[[i+j]]=namelist2[[j]]}
}
print(nameList)

###Checking The Common Orders and Families in between GSEA method and F-test
intersect(list1, nameList)
#Camera vs F-test
intersect(list2,nameList)


