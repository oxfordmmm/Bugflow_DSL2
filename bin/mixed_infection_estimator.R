#!/usr/bin/env Rscript

#################################################################################################################
# Mixed infection estimator                                                                                     #
# An algorithm for detection of mixed infections in bacterial whole genome sequence data.                       #
# The algorithm analyses a set of defined variable sites for evidence of mixed infection with up to 2 strains,  #
# if a mixed infection is present, the relative proportions of the dominant and minor haplotypes are estimated. #
# This version omits reporting of the dominant and minor haplotypes                                             #
#                                                                                                               #
# david.eyre@ndm.ox.ac.uk                                                                                       #
# 04 March 2013, update 21 June 2023                                                                                               #
#################################################################################################################

library("doMC")

args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
mixed_site_file <- args[2]
outlog <- args[3]
bs.iter <- strtoi(args[4]) #number of bootstrap iterations to run for parameter estimates

### OPTIONS ###

#number of cores to use for bootstrapping on Mac / Linux, ignored on Windows
bs.cores = 1

#base error probability, probability that any one base call represents an error
p.err = 2e-3 


### END OPTIONS ###



### CONSTANTS ###
DIPLOID = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
BASES = c("A", "C", "G", "T")


### FUNCTIONS ###
#functions for defining m, mixture proportion and d, between haplotype diversity by maximum likelihood

#define functions to keep values between 0.5 and 1
get.logit2 = function(p) log((2*p-1)/(1-(2*p-1)))
get.inv.logit2 = function(p) (  ( (exp(p) / (1 + exp(p))) /2 ) + 0.5)

#define functions to keep values between 0 and 1
get.logit = function(p) log(p/(1-p))
get.inv.logit = function(p) exp(p) / (1 + exp(p)) 

#several functions to generate the likelihood of m, the mixture proportion
get.p.b = function(A, m, b, e) {
	# returns p(b|A,m,e) where A = {AA, AC, ..., TT} and b is a base in a single read, and e is the read error probability
	p.b.1 = m * e/3 # p where major ≠ b 
	p.b.2 = (1-m) * e/3 # p where minor ≠ b
	if(b == substring(A,1,1)) p.b.1 = m * (1-e) # p where major = b
	if(b == substring(A,2,2)) p.b.2 = (1-m) * (1-e) # p where minor = b
	p.b = p.b.1 + p.b.2
	return(p.b)
	}

get.p.b.diploid = function(m) {
	# returns a matrix of p(b|A,m,e) for all values of b {A, C, G, T} - rows and values of A - columns
	# e is set by constant above
	mat = matrix(NA, nrow=4, ncol=16)
	for (i in 1:4) mat[i,] = sapply(DIPLOID, get.p.b, m=m, b=BASES[i], e=p.err)
	return(mat)
	}

get.p.site = function(site, m, bc, mat) {
	# returns p(B1 to Bn | A, m, e) where B1 to Bn are all reads at a given site
	# - product of p(b|A,m,e) across all reads - taken from base counts
	# where mat is a 16x4 matrix from get.p.b.diploid
	exp(colSums(log(mat)*bc[site,2:5]))
	}

get.d = function(d) {
	# returns a vector of probabilities for the 16 DIPLOID forms
	aa = (1-d)/4
	ab = d/12
	c(aa, ab, ab, ab, ab, aa, ab, ab, ab, ab, aa, ab, ab, ab, ab, aa)
	}

logliki.m = function(parm, bc) {
	# return the log likihood of m,d - product over all sites [ sum over all values of A [ p(B1-Bn | A, m, e)  * p(A) ] ]
	# p(A) from the vector given by get.d
	m = get.inv.logit2(parm[1])
	d = get.inv.logit(parm[2])
	mat = get.p.b.diploid(m)
	i = 1:nrow(bc)
	sum ( sapply(i, function(i,m,bc,mat,d) log(sum(get.p.site(i,m,bc,mat)*get.d(d))), m=m, bc=bc, mat=mat, d=d) )
	}

get.dev = function(ml, bc) {
	#compare the likelihood obtained to the likelihood under an approximation to the null hypothesis
	h0.liki = logliki.m(c(get.logit2(0.9999),get.logit(.00001)),bc)
	D = -2*(h0.liki-ml)
	return(D)
	}


sample.bc = function(bc) {
	#function for sampling base counts
	bc.samp = bc[sample(seq(1:nrow(bc)), nrow(bc), replace=TRUE),]
	return(bc.samp)
}








### ANALYSE THE DATA ###



#set up output header
cat("id\tML_m\tML_d\tML_sites_diff\tdeviance\tm0.025\tm0.975\td0.025\td0.975\n", file=outlog, append=FALSE)
#id, estimate of m, d, deviance statistic, lower bound m, upper bound, lower bound on d, upper bound, haplotype pairs if mixed


f = mixed_site_file
  
#check file is present
if (is.na(file.info(f)$size)) {
    print(paste("Missing file,", f, ", for id",id,"skipping...."))
} else {
    #print info message to screen
    print(paste("Starting ",id,"....",sep=""))
    
    #read in base count data
    base.count = as.matrix(read.table(f,header=T,sep="\t"))

    #if no base counts skip the file
    if(nrow(base.count)==0)
    {
        cat(id, 0, 0, 0, 0, 0, 0, 0, "0\n", sep="\t", file = outlog, append=TRUE)
        stop(paste("Base counts missing for id",id,"skipping...."))
    }
            
    #if any of the sites have a base count of zero skip the file
    if(min(rowSums(base.count[,2:5]))==0)
    {
        stop(paste("Base counts of 0 for at least one site in id",id,"skipping...."))
    }
    
    #ensure base counts are sorted in site order in genome
    base.count[rank(base.count[,1]),] = base.count[c(1:nrow(base.count)),]
    
    
    #obtain initial esimates of d and m
    #estimate of d from total heterozygous sites
    sites = nrow(base.count)
    d.obs = sum(rowSums(base.count[,2:5]==0)<3)
    if (d.obs==0) d.obs=0.00001
    (d.init = d.obs/sites)
    
    #estimate for starting value of m from mean across heterozygous sites
    m.init = 0.999
    if (sum(rowSums(base.count[,2:5]==0)<3)>0)
        {
        bc.het = base.count[rowSums(base.count[,2:5]==0)<3,2:5] #base counts with at least one heterozygote
        if (sum(rowSums(base.count[,2:5]==0)<3)==1)
        {
            (m.init = max(bc.het)/sum(bc.het))
        }
        else
        {
            bc.het.max = sapply(1:nrow(bc.het), function(i) max(bc.het[i,]))
            bc.het.tot = rowSums(bc.het)
            (m.init = mean(bc.het.max/bc.het.tot))
        }
        }
    
    #avoid exact 50/50 mix
    if (m.init == 0.5)
    {
        m.init = 0.50001
    }
    
    init_parm = c(get.logit2(m.init), get.logit(d.init)) #take initial parameters and express in terms of logit functions above to constrain optimisation values
    
    ###numerically optimise the values of m and d and obtain ML
    opt.md = optim(init_parm, logliki.m, control=list(fnscale=-1), bc=base.count)
    opt.md
    #ML value for m
    md.ML.m = get.inv.logit2(opt.md$par[1])
    #ML value for d
    md.ML.d = get.inv.logit(opt.md$par[2])
    
    #get deviance from null hypothesis
    md.dev = get.dev(ml=opt.md$value, bc=base.count)
    
    print(paste("got ML m:", md.ML.m, "and d:", md.ML.d, " (deviance: ",md.dev,")"))
    
    
    #run bootstrap, use multiple cores on Mac / Linux
    if (Sys.info()['sysname']!="Windows") {
        # multiple core support for bootstrap
        registerDoMC(cores=bs.cores)
        bs.md = foreach(i = 1:bs.iter, .combine=rbind) %dopar% {
            bc.samp = sample.bc(base.count)
            opt.md.samp = optim(init_parm, logliki.m, control=list(fnscale=-1), bc=bc.samp)
            bs.m = get.inv.logit2(opt.md.samp$par[1])
            bs.d = get.inv.logit(opt.md.samp$par[2])
            c(bs.m, bs.d)
            }
    } else {
        bs.md = foreach(i = 1:bs.iter, .combine=rbind) %do% {
            bc.samp = sample.bc(base.count)
            opt.md.samp = optim(init_parm, logliki.m, control=list(fnscale=-1), bc=bc.samp)
            bs.m = get.inv.logit2(opt.md.samp$par[1])
            bs.d = get.inv.logit(opt.md.samp$par[2])
            c(bs.m, bs.d)
        }		
    }
        
        
    bs.m.ci = quantile(bs.md[,1], probs=c(0.025, 0.975))
    bs.d.ci = quantile(bs.md[,2], probs=c(0.025, 0.975))
    print (paste("m (", bs.m.ci["2.5%"], " - ", bs.m.ci["97.5%"], ")"))
    print (paste("d (", bs.d.ci["2.5%"], " - ", bs.d.ci["97.5%"], ")"))
    
    cat(id, md.ML.m, md.ML.d, md.ML.d*nrow(base.count), md.dev, bs.m.ci["2.5%"], bs.m.ci["97.5%"], 
        bs.d.ci["2.5%"], bs.d.ci["97.5%"], sep="\t", file = outlog, append=TRUE)
    cat("\n", file = outlog, append=TRUE)
}
