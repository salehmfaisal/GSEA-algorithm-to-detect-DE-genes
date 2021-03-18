# OM: full Bayesian version "my method"
#   
#   get ROC curves using edgeRoc1, baySroc1 etc.
# 
#   version C starts with implementing "filtering" (dropping low counts)
#
#   version D will fit a Mean/variance Relationship ("MVR") 
#     version DG will introduce Gene Sets!
#       version DG1 was used for grant18, comparison with GSEA and others.
#       version DG2 is used for real data ("ALL" dataset?), microarray -- no need for filtering and such. 
#        computationally difficult: 5499 gene sets (after filtering!) and 9010 genes. Can still have a matrix that big! 
#          version DG2a: checking p-values bu permuting sample labels 

     loglkhMVR <- function(b, nu){    # b is the parameter vector  
      m <- Nuse
	  s0s <- (b[1] + b[2]/(A.use - b[3]))^4
      out <- nu/2*sum(log(s0s*nu/2)) - m*lgamma(nu/2) - (2 + nu)*sumslog - nu/2*sum(s0s/sigmA^2)
      out
    }   
   
 
 #  setwd('C:/localH/research/gibbseq')  
 #  setwd('C:/local/Reiss/Daniel/OMmay5_2016')
  #source('fullBrocDG2.r')  #  source('fullBrocBnb.r')  previous version 
 
  # require(edgeR) 
  #require(MASS)   # contains "rnegbin" function
  #require(DESeq2) 
  require(limma)  
   
   SelBS.all <- NULL;    SelOM.all <- NULL;  SelERT.all <- NULL
   SelDE2.all <- NULL
   FdBS.all <- NULL;    FdOM.all <- NULL; FdERT.all <- NULL
   FdDE2.all <- NULL;
   
   qvalERT <- NULL;  qvalOM <- NULL; qvalBS <- NULL
   qvalDE2 <- NULL
      # keeping q-values for saving the results of multiple runs  
   
  # N and G already precomputed, gene results in xx, group ID's in allSE$GROUP
   

 ptm <- proc.time()
    
  ## start implementing "OM"
     
   set1 <- (1:n1)
   set2 <- (1+n1):(n1+n2)   # columns of sData to be compared

   A.use <- rowMeans(xx)              #A.use <- A[g.use]
   in12 <- 1/n1 + 1/n2
   Nuse <- N
   ng = numeric(g)       # number of genes in each set
   for (g in 1:G){
      ng[g] = length(mem[[g]])
   }

    # ------------------------------------------------------------  Group stuff
         # currently done in "tryDataset.r"
    #-------------------------------------------
 	  
   # xx <- xx[g.use,]
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
  mean(havg > 1 - 0.05)    # low q-values
  
 # mean(Sigmahist)
 # plot(colMeans(Sigmahist),type="l")
 # sort(havg, decreasing = T)[1:75]
  
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
	
 # [1] "GO:0001765_membrane_raft_assembly"                                                          
 # [2] "GO:0001960_negative_regulation_of_cytokine-mediated_signaling_pathway"                      
 # [3] "GO:0002576_platelet_degranulation"                                                          
 # [4] "GO:0006636_unsaturated_fatty_acid_biosynthetic_process"                                     
 # [5] "GO:0006909_phagocytosis"                                                                    
 # [6] "GO:0006915_apoptotic_process"                                                               
 # [7] "GO:0006919_activation_of_cysteine-type_endopeptidase_activity_involved_in_apoptotic_process"
 # [8] "GO:0007159_leukocyte_cell-cell_adhesion"                                                    
 # [9] "GO:0008360_regulation_of_cell_shape"                                                        
# [10] "GO:0008630_intrinsic_apoptotic_signaling_pathway_in_response_to_DNA_damage"                 
# [11] "GO:0009062_fatty_acid_catabolic_process"                                                    
# [12] "GO:0009615_response_to_virus"                                                               
# [13] "GO:0010737_protein_kinase_A_signaling"                                                      
# [14] "GO:0016446_somatic_hypermutation_of_immunoglobulin_genes"                                   
# [15] "GO:0019221_cytokine-mediated_signaling_pathway"                                             
# [16] "GO:0019732_antifungal_humoral_response"                                                     
# [17] "GO:0030035_microspike_assembly"                                                             
# [18] "GO:0031397_negative_regulation_of_protein_ubiquitination"                                   
# [19] "GO:0032020_ISG15-protein_conjugation"                                                       
# [20] "GO:0032060_bleb_assembly"                                                                   
# [21] "GO:0032507_maintenance_of_protein_location_in_cell"                                         
# [22] "GO:0033137_negative_regulation_of_peptidyl-serine_phosphorylation"                          
# [23] "GO:0033153_T_cell_receptor_V(D)J_recombination"                                             
# [24] "GO:0034113_heterotypic_cell-cell_adhesion"                                                  
# [25] "GO:0034138_toll-like_receptor_3_signaling_pathway"                                          
# [26] "GO:0034260_negative_regulation_of_GTPase_activity"                                          
# [27] "GO:0034340_response_to_type_I_interferon"                                                   
# [28] "GO:0035455_response_to_interferon-alpha"                                                    
# [29] "GO:0035456_response_to_interferon-beta"                                                     
# [30] "GO:0035666_TRIF-dependent_toll-like_receptor_signaling_pathway"                             
# [31] "GO:0035912_dorsal_aorta_morphogenesis"                                                      
# [32] "GO:0036120_cellular_response_to_platelet-derived_growth_factor_stimulus"                    
# [33] "GO:0043066_negative_regulation_of_apoptotic_process"                                        
# [34] "GO:0043123_positive_regulation_of_I-kappaB_kinase/NF-kappaB_signaling"                      
# [35] "GO:0043124_negative_regulation_of_I-kappaB_kinase/NF-kappaB_signaling"                      
# [36] "GO:0043297_apical_junction_assembly"                                                        
# [37] "GO:0043312_neutrophil_degranulation"                                                        
# [38] "GO:0043534_blood_vessel_endothelial_cell_migration"                                         
# [39] "GO:0045071_negative_regulation_of_viral_genome_replication"                                 
# [40] "GO:0045087_innate_immune_response"                                                          
# [41] "GO:0045629_negative_regulation_of_T-helper_2_cell_differentiation"                          
# [42] "GO:0045792_negative_regulation_of_cell_size"                                                
# [43] "GO:0046135_pyrimidine_nucleoside_catabolic_process"                                         
# [44] "GO:0046632_alpha-beta_T_cell_differentiation"                                               
# [45] "GO:0046825_regulation_of_protein_export_from_nucleus"                                       
# [46] "GO:0048505_regulation_of_timing_of_cell_differentiation"                                    
# [47] "GO:0048813_dendrite_morphogenesis"                                                          
# [48] "GO:0050778_positive_regulation_of_immune_response"                                          
# [49] "GO:0050798_activated_T_cell_proliferation"                                                  
# [50] "GO:0050832_defense_response_to_fungus"                                                      
# [51] "GO:0051017_actin_filament_bundle_assembly"                                                  
# [52] "GO:0051092_positive_regulation_of_NF-kappaB_transcription_factor_activity"                  
# [53] "GO:0051292_nuclear_pore_complex_assembly"                                                   
# [54] "GO:0051493_regulation_of_cytoskeleton_organization"                                         
# [55] "GO:0051607_defense_response_to_virus"                                                       
# [56] "GO:0051673_membrane_disruption_in_other_organism"                                           
# [57] "GO:0051764_actin_crosslink_formation"                                                       
# [58] "GO:0060337_type_I_interferon_signaling_pathway"                                             
# [59] "GO:0060993_kidney_morphogenesis"                                                            
# [60] "GO:0071222_cellular_response_to_lipopolysaccharide"                                         
# [61] "GO:0090382_phagosome_maturation"                                                            
# [62] "GO:0097190_apoptotic_signaling_pathway"                                                     
# [63] "GO:1900027_regulation_of_ruffle_assembly"                                                   
# [64] "GO:1900087_positive_regulation_of_G1/S_transition_of_mitotic_cell_cycle"                    
# [65] "GO:1901216_positive_regulation_of_neuron_death"                                             
# [66] "GO:1902715_positive_regulation_of_interferon-gamma_secretion"                               
# [67] "GO:1905244_regulation_of_modification_of_synaptic_structure"                                
# [68] "GO:2000096_positive_regulation_of_Wnt_signaling_pathway,_planar_cell_polarity_pathway"      
# [69] "GO:2000134_negative_regulation_of_G1/S_transition_of_mitotic_cell_cycle"                    
# [70] "GO:2000811_negative_regulation_of_anoikis"                                                  
# [71] "GO:2001205_negative_regulation_of_osteoclast_development" 
	
   # However, with a different prior (a=1, b=9) the list is only 27 long, pDE ~= 0.075 
   
 # [1] "GO:0002576_platelet_degranulation"                                                          
 # [2] "GO:0006909_phagocytosis"                                                                    
 # [3] "GO:0006915_apoptotic_process"                                                               
 # [4] "GO:0006919_activation_of_cysteine-type_endopeptidase_activity_involved_in_apoptotic_process"
 # [5] "GO:0007159_leukocyte_cell-cell_adhesion"                                                    
 # [6] "GO:0007229_integrin-mediated_signaling_pathway"                                             
 # [7] "GO:0008360_regulation_of_cell_shape"                                                        
 # [8] "GO:0009615_response_to_virus"                                                               
 # [9] "GO:0030035_microspike_assembly"                                                             
# [10] "GO:0031397_negative_regulation_of_protein_ubiquitination"                                   
# [11] "GO:0036120_cellular_response_to_platelet-derived_growth_factor_stimulus"                    
# [12] "GO:0038083_peptidyl-tyrosine_autophosphorylation"                                           
# [13] "GO:0043066_negative_regulation_of_apoptotic_process"                                        
# [14] "GO:0043123_positive_regulation_of_I-kappaB_kinase/NF-kappaB_signaling"                      
# [15] "GO:0043124_negative_regulation_of_I-kappaB_kinase/NF-kappaB_signaling"                      
# [16] "GO:0043312_neutrophil_degranulation"                                                        
# [17] "GO:0045071_negative_regulation_of_viral_genome_replication"                                 
# [18] "GO:0045087_innate_immune_response"                                                          
# [19] "GO:0046632_alpha-beta_T_cell_differentiation"                                               
# [20] "GO:0051017_actin_filament_bundle_assembly"                                                  
# [21] "GO:0051092_positive_regulation_of_NF-kappaB_transcription_factor_activity"                  
# [22] "GO:0051607_defense_response_to_virus"                                                       
# [23] "GO:0051764_actin_crosslink_formation"                                                       
# [24] "GO:0060337_type_I_interferon_signaling_pathway"                                             
# [25] "GO:0071222_cellular_response_to_lipopolysaccharide"                                         
# [26] "GO:0097190_apoptotic_signaling_pathway"                                                     
# [27] "GO:2000096_positive_regulation_of_Wnt_signaling_pathway,_planar_cell_pola
  
  
  
  ## running GSEA and others?  ===============================================================================================
  
  
  
  