	# looking for a nice dataset to analyze 
#####To install this package, start R (version "3.6") and enter: 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ALL", version = "3.10")
######To install this package, start R (version "3.6") and enter: 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("hgu95av2.db", version = "3.10")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EnrichmentBrowser", version = "3.10")
###Start
  require(EnrichmentBrowser)
  library(ALL)
  data(ALL)     # this is a human dataset 
  
  ind.bs <- grep("^B", ALL$BT)
  ind.mut <- which(ALL$mol.biol %in% c("BCR/ABL", "NEG"))
  sset <- intersect(ind.bs, ind.mut)
  all.eset <- ALL[, sset]
  
   dim(all.eset)
   exprs(all.eset)[1:4,1:4]
   allSE <- probe2gene(all.eset)  # installs  hgu95av2.db and other databases
    head(rownames(allSE))
	

    allSE$GROUP <- ifelse(allSE$mol.biol == "BCR/ABL", 1, 0)
    table(allSE$GROUP)	    # sample sizes 42 and 37 
  
    
 	 ## =======================================================================
	 ## ============  preparing for the List-making 
	 
	 xx = assay(allSE)
	 geneid = rownames(allSE)      # Entrez ID's from the data table 
	 
	 # dat = "C:/localH/research/gibbseq/data/hsa.filt.gmt"
   	 # cc = count.fields(dat, sep = '\t')	
     # dat1 = readLines(dat)
     # dat1 = strsplit(dat1, "\t")
	 # nL = length(dat1) 
	 # nameGS = rep('',nL)    # 5499 genesets in "hsa.filt.gmt"
	 # for (i in 1:nL){
	    # nxtline =  dat1[[i]]
	    # nameGS[i] = nxtline[1]  # just the name of geneset 
    # }
	 
	 hsa.filt.gs = getGenesets("hsa.filt.gmt") 
	 nameGS = names(hsa.filt.gs)  ## !! this replaces the above loop
	 nameList = as.list(nameGS) 
	 	 
	   # hsa.filt.gs[[10]]      "3692" "5901" "6059" "7514"    hsa.... = mem?  or need to translate into 1,..., G? 
	 
     #  as it stands, translate both *mem* and Gene names (*geneid*) into 1,..., G 
	 
 ### >>>>  begin procedure  ---------------------------- >-- 	
 
	 G = length(hsa.filt.gs)
	 N = length(geneid) 
	 mem = vector("list", G)     # null list
	 memID = vector("list", G)   # list of Geneset names 
	 for (k in 1:N){   # renaming genes and rewriting Gene set list "mem"
	   id = geneid[k]
       for (g in 1:G){
	     if (id %in% hsa.filt.gs[[g]]){ 
            mem[[g]] = c(mem[[g]], k)   # change from "id" to "k"
	     }  
	   }
     }
	setlen = numeric(G) 
    for (g in 1:G){
	   setlen[g] = length(mem[[g]])
	   
	}
	# note: not all ID's in "hsa.filt.gs" are present in the data set!
	# there are 35 "sets" that have no entries;   4380 sets have length > 3 
	 sum(setlen > 3)
	 
	for (g in G:1){                      # remove from list, have to do it from the back!
	   if( setlen[g] <= 3){
	     mem[[g]] = NULL
		 nameList[[g]] = NULL
	   } 
	}
	G = length(mem)	
	nameArr = as.character(nameList)
	
	# sort columns of xx according to group

   n1 = sum(allSE$GROUP == 0);   n2 = sum(allSE$GROUP == 1)
   #xx1 = xx[,allSE$GROUP == 0]
   #xx2 = xx[,allSE$GROUP == 1]    xx = cbind(xx1,xx2)      # not very economical, but easy 
    col1 = which(allSE$GROUP == 0);   col2 = which(allSE$GROUP == 1);
   xx = xx[,c(col1,col2)]    # this is better    
 
  
  ### <<<< end procedure----------------------------- >-- 	 
  
