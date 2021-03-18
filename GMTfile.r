######my data 
##Starting of the coding
data=read.csv("NRAP_JGI_Species.dat")
###Data Preparation
####For Order
h.order=as.factor(data$Order)
Orders=levels(h.order)
head(Orders)
L=length(Orders)
mem1 = NULL
nameList1=NULL 
nameArr1=NULL 
for (k in 1:L){
	   id = Orders[k]
       mem1[[k]]=which(data$Order==id)
nameList1[[k]]=id
   nameArr1[[k]]=paste("order:",id)  
}


###For Family
h.family=as.factor(data$Family)
Families=levels(h.family)
head(Families)
M=length(Families)
mem2 = NULL
nameList2=NULL 
nameArr2=NULL 
for (d in 1:M){   
	   id1 = Families[d]
       mem2[[d]]=which(data$Family==id1)
nameList2[[d]]=id1
   nameArr2[[d]]=paste("family:", id1) 
          
     }


##Merging Two mem and Name List
nameList=NULL
mem=NULL
for(i in 1:length(mem1)){
mem[[i]]=mem1[[i]]
nameList[[i]]=nameList1[[i]]
for (j in 1:length(mem2)){
mem[[i+j]]=mem2[[j]]
nameList[[i+j]]=nameList2[[j]]}
}
print(nameList)



 N = length(mem)
	 
	 sink("microbialset.gmt")	 
	 for (i in 1:N){  
			L0 = paste(mem[[i]], collapse = "\t")
			L = paste(nameList[[i]],"NA",L0,sep="\t")
			cat(L)
			cat("\n")		
	 }
	 sink()

