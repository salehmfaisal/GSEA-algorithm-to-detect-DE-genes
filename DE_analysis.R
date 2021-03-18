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
##Preparing xx data matrix
xx=data[,8:11]
#Normalization
xx$SAE3_00=log(xx[,1]+1)/mean(log(xx[,1]+1))
xx$SAE3_05=log(xx[,2]+1)/mean(log(xx[,2]+1))
xx$SAE3_23=log(xx[,3]+1)/mean(log(xx[,3]+1))
xx$SAE3_37=log(xx[,4]+1)/mean(log(xx[,4]+1))
head(xx)


#####Function
loglkhMVR <- function(b, nu){    # b is the parameter vector  
      m <- Nuse
	  s0s <- (b[1] + b[2]/(A.use - b[3]))^4
      out <- nu/2*sum(log(s0s*nu/2)) - m*lgamma(nu/2) - (2 + nu)*sumslog - nu/2*sum(s0s/sigmA^2)
      out
    }  

###sum(log(s0s*nu/2)) gives me INf.

bacid = rownames(data)
##Implimentation of Code to xx and mem
## start implementing "OM"
G=length(mem)
N=length(bacid)
nameArr=as.character(nameList)
n1=2
n2=2
set1 <- (1:n1)
   set2 <- (1+n1):(n1+n2) 
A.use <- rowMeans(xx)              
   in12 <- 1/n1 + 1/n2
   Nuse <- N
   ng = numeric(i+j) 
for (g in 1:G){
      ng[g] = length(mem[[g]])##Number of member in a mem.
   }

sd1 <- apply(xx[,set1],1,sd)
   sd2 <- apply(xx[,set2],1,sd)
S2x <- sd1^2  #@@@ another option here: try pooled variance estimates (Sx, jointly) ==> this will also help with n_B = 1   
	  S2y <- sd2^2       # @ @ @ @ @ @ 
 M <- rowMeans(xx[,set1]) - rowMeans(xx[,set2])  
 plot(A.use, M)  


 # initialize sample variances for each gene     
   sigmA <- rep(0.2, Nuse)
   tau <- 0.9            # st.dev. of the Delta = "jump" distribution
   
   MC <- 1000    # number of Monte Carlo runs. MC = 1000, nskip = 5  takes about 3 min to run. 
   nskip <- 5       ## >>   do not have highly correlated samples here? But need to average the results still, to get hgavg!
   burnin <- 100
   
   h <- rep(0,Nuse) 
   havg <- h                # how often the gene was in the "significant" portion
   hgavg = rep(0,G)
   pDE <- 0.1
   pDEg <- rep(0.1, G)
   alphap <- 1 
   betap <- 9         # beta prior parameters for pDE distribution
   df0.tau <- 10
   v0.tau <- 1          # inv-Xi-sq. prior for tau 
    s0sq <- 0.17^2; nu <- 3.5  	#  initial parameters for the hyper-distribution of sigmas
	b <- c(0.5, 1, 0)               # initialize MVR parameters b
	
	   jump.nu <- 0.1; jump.s0 <- 0.2    # parameters for Metropolis
	   jump.b0 <- 0.05;  jump.b1 <- 0.05; jump.b2 <- 0.05;
	acc.ctr0 <- 0; acc.ctr1 <- 0;  acc.ctr2 <- 0; acc.ctrnu <- 0;
   
   Sigmahist <- matrix(0, Nuse, MC)   # storing the MC results 
   s0hist <- rep(0,MC) 
   nuhist <- rep(0,MC) 
   pDEhist <- rep(0, MC);   pDEghist <- matrix(0,G,MC)
   tauhist <- rep(0, MC) 
   bhist <- matrix(0, length(b), MC) 
   
   prior.df <- 3;  prior.var <- 0.04   # for sigma distribution

for (mc in 1:MC){
   
     for (iskip in 1:nskip){
   
     # Gibbs step for h_k
	 
	  pDEvec <- rep(0,Nuse);   Nmem <- rep(0,Nuse) 	  
	  for (g in 1:G){
	    memg <- mem[[g]]            
	     pDEvec[memg] <- pDEvec[memg]  + pDEg[g]
		 Nmem[memg] <- Nmem[memg] + 1   		 # + one way to go about averaging
	  }	 
	  Nmem[pDEvec ==0] <- 1
	  pDEvec[pDEvec ==0] <- pDE     # all "unaffiliated" genes
	   pDEvec <- pDEvec/Nmem
	 
	  sigtau <- sqrt(sigmA^2*in12 + tau^2)
	  ph0 <- dnorm(M, 0 , sigmA*sqrt(in12))*(1-pDEvec)  + 1E-12
	  ph1 <- dnorm(M, 0 , sigtau)*pDEvec + 1E-12
	  p1post <- ph1/(ph0 + ph1)                    # a problem if both are 0
	  h <- (runif(Nuse) < p1post)
	   
	 # FCP for D_k
       Dvar <- 1/(1/in12/sigmA^2 + 1/tau^2)
       D <- h * rnorm(Nuse, M/in12/sigmA^2*Dvar, sqrt(Dvar)) 	 
	      # D = 0 if h = 0, D is the "jump size"
	   
	 # FCP for sigma 
        ssq <- S2x*(n1-1) +  S2y*(n2-1)  + (M-D)^2/in12 
		s0sq <- (b[1] + b[2]/(A.use - b[3]))^4
		sigmA <- sqrt((nu * s0sq + ssq)/rchisq(Nuse,(n1 + n2 - 1) + nu))
	   
# ** # FCP for pDE     			====>  a flaw here since not all genes are "unaffil."?
      alpha1 <- alphap + sum(h)
	  beta1 <- betap + sum(1-h)
	  pDE <- rbeta(1, alpha1, beta1)   # this is now the "baseline" pDE
	  
	 # FCP for pDEg			==========> can't improve?
	 
	 for (g in 1:G){
	    memg <- mem[[g]]         # ng[g] = length(memg)?
        if (ng[g] > 0){
		  alpha1 <- alphap + sum(h[memg])
	      beta1 <- betap + sum(1-h[memg])
	      pDEg[g] <- rbeta(1, alpha1, beta1) 
        } else {
 		 pDEg[g] <- pDE
		}
     }		 
	    
	 # FCP for tau (variance of Delta_i)
	  tau <- sqrt((sum(D^2) + v0.tau*df0.tau)/rchisq(1,sum(h) + df0.tau))
	   
	  # FCP for nu and s0sq, variance hyperparameters  

	 sumslog <- sum(log(sigmA))
	 lnb0 <- log(b[1])
     lnb1 <- log(b[2])
	 lnb2 <- log(- b[3])
	 lnnu <- log(nu)
     nu.prop  <- exp(lnnu + jump.nu*runif(1,-1,1) )
	 b0.prop  <- exp(lnb0 + jump.b0*runif(1,-1,1) )
	 b1.prop  <- exp(lnb1 + jump.b1*runif(1,-1,1) )
	 b2.prop  <-  - exp(lnb2 + jump.b2*runif(1,-1,1) )    # later, replace this 2 by a "thresh" or something 

	 
     # accept/reject for log(sigma_0^2) and log(nu)
	    b.prop <- c(b0.prop, b[2], b[3])
        acc.prob1 <- exp(min(0, loglkhMVR(b.prop,nu) - loglkhMVR(b,nu) ))
		if (runif(1) < acc.prob1){
		    b <- b.prop
			acc.ctr0 <- acc.ctr0 + 1
 	     }
	    b.prop <- c(b[1], b1.prop, b[3])
		acc.prob1 <- exp(min(0, loglkhMVR(b.prop,nu) - loglkhMVR(b,nu) ))
		if (runif(1) < acc.prob1){
		    b <- b.prop
			acc.ctr1 <- acc.ctr1 + 1
 	     }
		b.prop <- c(b[1], b[2], b2.prop)
		acc.prob1 <- exp(min(0, loglkhMVR(b.prop,nu) - loglkhMVR(b,nu) ))
		if (runif(1) < acc.prob1){
		    b <- b.prop
			acc.ctr2 <- acc.ctr2 + 1
 	     }
				
		acc.prob2 <- exp(min(0, loglkhMVR (b,nu.prop) - loglkhMVR (b,nu) ))
		if (runif(1) < acc.prob2){
		    nu <- nu.prop
			acc.ctrnu <- acc.ctrnu + 1
 	     } 
	  if (mc > burnin){	
        havg <- havg + h/((MC-burnin)*nskip) 	
        hgavg = hgavg + (pDEg > pDE)/((MC-burnin)*nskip) 		
	  }	
	 } # end skip  
	   
	 # output bloc

	 pDEhist[mc] <- pDE
	 pDEghist[,mc] <- pDEg - pDE    # recording the difference, to get a q-value 
	 tauhist[mc] <- tau
	 Sigmahist[,mc] <- sigmA
	 #s0hist[mc] <- sqrt(s0sq)
	 bhist[,mc] <- b
	 nuhist[mc] <- nu
	 
   }  # end of Gibbs
 par(mfrow=c(1,4))
  plot(havg)

   
  acc.f1 <- acc.ctr1/(MC*nskip);     acc.f2 <- acc.ctr2/(MC*nskip);
	acc.f0 <- acc.ctr0/(MC*nskip);    acc.nu <- acc.ctrnu/(MC*nskip) 
		c(acc.f0, acc.f1, acc.f2, acc.nu)

 if (0==1){		
	windows(12,7)
	par(mfrow=c(1,2))
	plot(s0hist, type="l") 
    plot(nuhist, type="l") 
	
	windows(12,7); 	par(mfrow=c(1,2))
	plot(pDEhist, type="l") 
    plot(tauhist, type="l") 
 }	
  
if (0==1){	  qval = 1 - havg
	 	#qval <- matrix(0,N,1) 
 		#qval[g.use] <- 1 - havg     #> for the filtering version: adjust them so that the filtered-out genes get q-value of 1
		#qval[!g.use] <- 1                    		
		
		logq <- -log10(qval + 1E-6)	
	
  th <- seq((min(qval)+1e-11)*1.001, max(qval)*0.999, 0.01)      
          #  use the same th as for edgeR
nth <- length(th)
SelOM <- numeric(nth); FdOM <- numeric(nth)

for (i in 1:nth){
  rej <- (qval < th[i])   # the number of "Selected" (rejected) genes
  XX <- table(rej, de.list) 
  SelOM[i] <- sum(XX[2,]) 
  FdOM[i] <- XX[2,1]
}  
}


print(proc.time() - ptm)

# -----------------------------------------------------------
# Producing plots for all the methods 
#windows(8,6); 



if (0==1){                             # turn these off for now
plot(SelOM,FdOM+0.001,log="y")
  title(paste("Run no.",irep))
  lines(SelERT,FdERT+0.001, col="blue")
  lines(SelDE2,FdDE2+0.001, col="green")

SelOM.all <- c(SelOM.all, SelOM); SelERT.all <- c(SelERT.all, SelERT);  SelDE2.all <- c(SelDE2.all, SelDE2) 
FdOM.all <- c(FdOM.all, FdOM); FdERT.all <- c(FdERT.all, FdERT);  FdDE2.all <- c(FdDE2.all, FdDE2) 
qvalOM <- c(qvalOM, qval)    

sOM <- smooth.spline(FdOM.all, SelOM.all, df=50)   
ycOM <- c(seq(0.01,1,0.01), seq(1,N,10))
xcOM <- predict(sOM, ycOM)$y

sERT <- smooth.spline(FdERT.all, SelERT.all, df=50)   
ycERT <- c(seq(0.01,1,0.01), seq(1,N,10))
xcERT <- predict(sERT, ycERT)$y

sDE2 <- smooth.spline(FdDE2.all, SelDE2.all, df=50)   
ycDE2 <- c(seq(0.01,1,0.01), seq(1,N,10))
xcDE2 <- predict(sDE2, ycDE2)$y


plot(SelOM.all,FdOM.all+0.001,log="yx",col="grey70", xlab="All positives", ylab = "False Positives")
 lines(xcOM, ycOM, lwd=2, col = "red")
 lines(xcERT, ycERT, lwd=2, col = "blue", lty = 2)
  lines(xcDE2, ycDE2, lwd=2, col = "magenta", lty=3)
   legend(200,1,c("OM", "deSeq2","edgeR"), col=c("red", "blue", "magenta"), lty = c(1,2,3), lwd=2)
   

  runId <- "run1nb-b"
  write.csv(round(qvalOM,6),paste("qvalOM",runId,".csv",sep=""), row.names = FALSE)
  write.csv(round(qvalERT,6),paste("qvalERT",runId,".csv",sep=""), row.names = FALSE)
   write.csv(round(qvalDE2,6),paste("qvalBS",runId,".csv",sep=""), row.names = FALSE)
  }
  
 #points(Sel.all,Fd.all+0.001, col="yellow")

 meanPDEg <- apply(pDEghist,1, mean)
 qvalG <- 1 - apply((pDEghist >0),1, mean)    ## >> redo this based on hgavg
 plot(qvalG)
  
	
	# checking how "real" the FDR is 

qvFDR <- function(qval){	
   # th <- seq((min(qval)+1e-11)*1.001, max(qval)*0.999, 0.01)      
          #  use the same th as for edgeR
	nth <- length(th)
	Sel <- numeric(nth); Fd <- numeric(nth)
	FDR <- Fd;   recal <- Fd

	for (i in 1:nth){
	  rej <- (qval < th[i])   # the number of "Selected" (rejected) genes
	  rej <- c(rej,TRUE)             # to ensure no blank row in XX table
		de.list1 <- c(de.all,1)
	  XX <- table(rej, de.list1) 	
	  Sel[i] <- sum(XX[2,]) - 1 
	  Fd[i] <- XX[2,1]
	  recal[i] <- mean(qval[qval < th[i]])  
  	  # sum over all q-values smaller than threshold: "integrated OM". Does not work that well: may underestimate FDR. 
	}  
	 FDR <- Fd/(Sel+0.001)
	 func <- NULL
	 func$FDR <- FDR
	 func$rec <- recal
    return(func)
  }	
 
 th <- seq(0.0005,0.5,0.001)
 de.all <- rep(de.list, Nrep)   
 
 if (0==1){
  FDom <- qvFDR(qvalOM)
  FDRom <- FDom$FDR
 windows(5,5);   par(mar=rep(3,4))
 plot(th, FDRom, xlab="nominal", ylab="actual", col="red", asp=1, xlim=c(0,0.1), ylim=c(0,0.2))
 
   plot(th, FDRom, xlab="nominal", ylab="actual", col="red", asp=1, xlim=c(0,0.5), ylim=c(0,0.5))
    lines(c(0,1),c(0,1))
	points(th, qvFDR(qvalDE2)$FDR,  col="green",pch=16, cex=0.5)
 	points(th, qvFDR(qvalERT)$FDR,  col="blue",pch=16, cex=0.5)
	points(FDom$rec, FDRom, col="magenta", pch=6, cex=0.5)     # plots recalibrated q-values
 }
 ## recalculating the q-values 
 
                
  ## 32 bit: 166 sec/rep,  64 bit: 245 sec/rep ?? 

  # 	c(acc.f0, acc.f1, acc.f2, acc.nu)
  
  ## results for ALL dataset
  
    # sigG = which(qvalG < 0.05)
	
	plot(qvalG, 1 -hgavg)   
   sigG = which(1 - hgavg < 0.05)		
   nameArr[sigG]
	withORA = intersect(nameArr[sigG], set1)   # the list has 49 genesets in common!
	

   






