### NMR
# INTERVAL NMR traits: /rds/rds-jmmh2-post_qc_data/interval/phenotype/nightingale_metabolomics/qc/nmr_qcgwas.csv
#contains 40,849 observations and 230 traits
#exclude the "ace" trait, so 229 traits 

# INTERVAL GWAS: /home/ja628/rds/rds-jmmh2-results/private/metabolomics/nmr/raw_results/interval/ 

# INTERVAL imputed genetic data: /home/ja628/rds/rds-jmmh2-post_qc_data/interval/imputed



### Factor analysis of observed traits
Y <- read.csv("traits.csv") # N x P matrix of trait measurements, allowing for missing measurements
P <- ncol(Y)
N <- nrow(Y)

#check the number of missing values
count_non_missing_values <- as.data.frame(colSums(!is.na(Y)))

# remove any traits with significantly fewer observations


#check pairwise missing values
count_non_missing_values_pairwise = matrix(0,P,P)
colnames(count_non_missing_values_pairwise) <- rownames(count_non_missing_values_pairwise) <-  colnames(Data2018_trait_all[,-1])
for(i in 1:P){
  for(j in 1:P){
    check_pair = Y[,c(rownames(count_non_missing_values_pairwise)[i],
                                           colnames(count_non_missing_values_pairwise)[j])]
    check_count = sum(complete.cases(check_pair))
    count_non_missing_values_pairwise[i,j] = check_count
  }
}

# trait covariance
covY <- cov(Y, use="pairwise.complete.obs", method='pearson')

# check if covY is positive semi-definite and if not, then make it
	minev <- min(eigen(covY)$values)
	A <- covY
	if(minev < 0) {
	 diag(A) <- diag(A)+ abs(minev)+10^(-10)
	 covY <- A
	}

write.csv(covY, file='Trait_cov.csv')


# determine number of latent factors
## NOTE: in fa.parallel function n.obs = how many cases were used to find the correlations
median_nonmissing2 <- median(count_non_missing_values_pairwise)
library(psych)
fa.parallel(covY, fm="ml", n.obs=median_nonmissing2,fa="fa",show.legend=F,main="") # use median non-missing pairwise

## set nfactors to the number of latent factors suggested above
nfactors <-  

pdf("scree.pdf") # jpeg("scree.jpg") tiff("scree.tiff") 
fa.parallel(covY, fm="ml", n.obs=median_nonmissing2,fa="fa",show.legend=F,main="")
abline(v=nfactors,col="red")
dev.off()


# using number of latent factors from above, run factor analysis using covariance matrix
faY <- fa(r=covY, nfactors=nfactors, fm="ml", rotate="varimax")
faYloading <- as.matrix(faY$loadings)


####rescaled factors to get contributions 
latent_terms <- factor_contributions(faYloading) # output is factor loadings and re-scaled factor loadings (contributions) matrices with observed traits ordered by maximum contributing latent factor
rescaledFAloading <- latent_terms[[1]]
FAloading <- latent_terms[[2]]

write.csv(rescaledFAloading, file='Re-scaled-Factor_loading_matrix.csv')
write.csv(FAloading, file='Factor_loading_matrix.csv')
write.csv(covY[rownames(FAloading),rownames(FAloading)], file='Trait_cov.csv')


#### region-wise

### input gwas data from observed traits for region of interest
## obsgwas should be a list with P components (one for each trait)
## each data.frame should have the following columns, in any order, and may have others as well but they will not be used
## c("chromosome", "BP", "CHR", "EA", "NEA", "EAF", "beta","p_value","rsID")

fnames <- read.txt("paths_with_filenames_of_gwas_region_files.txt")
obsgwas <- vector("list", P)
for(i in 1:P) obsgwas[[i]] <- read.csv(fnames[i])

# harmonise gwas so variants are aligned to the same alleles
source("dataprep.R") 
RPinfo <- obsgwas[[1]] # align to first gwas
# need allele1, allele2, rsID as column names
ind <- which(colnames(RPinfo) %in% c("EA","NEA")) # match these to names in gwas files
colnames(RPinfo)[ind] <- c("allele1","allele2") 

tmp <- obsgwas
obsgwas <- lapply(tmp,alignGWAS,RPinfo=RPinfo)
rm(tmp)

# set same snps in all GWAS
intsnps <- obsgwas[[1]]$rsID
for(i in 2:length(obsgwas)) intsnps <- intersect(intsnps,obsgwas[[i]]$rsID)
for(i in 1:length(obsgwas)) obsgwas[[i]] <- obsgwas[[i]][which(obsgwas[[i]]$rsID %in% intsnps),]


# QC steps, INFO>0.4, MAF>0.005 have some lines on this in the Aim 1 file

### estimate latent factor GWAS summary stats from observed trait GWAS summary stats #####

source("latentGWAS.R")

covY <- read.csv("Trait_cov.csv",row.names=1) # trait covariance matrix 
FAloadings <- read.csv("Factor_loading_matrix.csv",row.names=1) # factor loadings
#covY <- covY[rownames(FAloadings),rownames(FAloadings)] # observed traits in same order as in FAloadings
obsGWAS <- obgwas[rownames(FAloadings)]


# latent gwas
latent_gwas = latentGWAS(obsGWAS = obsGWAS,covY = covY,L = FAloadings, beta_colname="beta",se_colname="SE",
                       snpID_colname="rsID",EAF = obsGWAS[[1]]$EAF)

# identify latent traits with a GWS signal in region
check <- sapply(latent_gwas,function(x) any(x$p_value<5E-8))
fm_traits_latent <- names(check[check==TRUE]) # latent traits with a GWS signal in region

# identify traits with contribution >0.20 from latent traits to be fine-mapped and that have a GWAS signal in region
rescaledFAloadings <- read.csv("Re-scaled-Factor_loading_matrix.csv",row.names=1) 
tmp <- rescaledFAloadings[,fm_traits_latent] 
linked_traits <- rownames(tmp)[which(apply(tmp,1, function(x) any(x>0.2)))] # identify traits with contribution >0.20 from latent traits to be fine-mapped
check <- sapply(obsGWAS[linked_traits],function(x) any(x$p_value<5E-8))
fm_traits_obs <- names(check[check==TRUE]) # obs traits with a GWS signal in region and contrib > 0.2 with latent factors fine-mapped


### LD and align gwas to coding in LD
library(bigsnpr)
library(vroom)
library(propagate)



fstub <- "path_to_bgen_file"
bgenfiles = paste0(fstub,'/bgen_13_regions/reg',id_reg,'_chr',id_chr,'_',id_start,'_',id_end,'_date20240213.bgen')
backingfile = paste0('bgen_13_regions/reg',id_reg,'_output_2')
list_snp_id=list()
list_snp_id[[1]] = paste0(obsGWAS[[1]]$CHR,'_',
                          obsGWAS[[1]]$BP,'_',
                          obsGWAS[[1]]$EA,'_',
                          obsGWAS[[1]]$NEA)

rds = snp_readBGEN(
  bgenfiles,
  backingfile,
  list_snp_id,
  ind_row = NULL,
  bgi_dir = dirname(bgenfiles),
  read_as = "dosage",
  ncores = 1
)

test <- snp_attach(rds)

RPinfo = test$map[,c(1,2,4,5,6)]
RPinfo = as.data.frame(RPinfo)
names(RPinfo)[2:3] = c("rsID","BP")  # c("chromosome","rsID","BP","allele1","allele2") 
RPinfo$rsID <- paste0("chr",chr,"_bp",RPinfo$BP,"_rsid",RPinfo$rsID)

chr <- RPinfo$chromosome[1]
genotypes = test$genotypes[1:dim(test$genotypes)[1],1:dim(test$genotypes)[2]]
colnames(genotypes) = RPinfo$rsID
genotypes[1:20,1:5]


# keep only traits that will be fine-mapped
tmp1 <- latent_gwas[fm_traits_latent] 
tmp2 <- obsGWAS[fm_traits_obs]

# align obsGWAS and latent_gwas_25 to alleles coded for LD
sig_gwas_latent <- lapply(tmp1,alignGWAS,RPinfo=RPinfo)
sig_gwas_obs <- lapply(tmp2,alignGWAS,RPinfo=RPinfo)

#ksnp <- sig_gwas_obs[[1]]$rsID[which(sig_gwas_obs[[1]]$MAF>.005)]
#isnp <- intersect(ksnp, colnames(genotypes) )

isnp <- intersect(sig_gwas_obs[[1]]$rsID, colnames(genotypes) )
X <- genotypes[,isnp]


Xqc <- apply(X,2,LDqc,theta=0.2,BestGuess=TRUE) # best-guess genotype matrix

Nprop <- apply(Xqc,2,function(x) sum(!is.na(x)))/nrow(Xqc)
keep <- which(Nprop >= 0.80)

Xqc <- Xqc[,keep]

for(i in 1:length(sig_gwas_obs)) sig_gwas_obs[[i]] <- sig_gwas_obs[[i]][keep,]
for(i in 1:length(sig_gwas_latent)) sig_gwas_latent [[i]] <- sig_gwas_latent[[i]][keep,] 

ldout <- bigcor(Xqc, size=min(2000, ncol(Xqc)), verbose=FALSE, fun= "cor", use="p" )
corX <- as.matrix(ldout[1:dim(ldout)[1],1:dim(ldout)[1]])
rownames(corX) <- colnames(corX) <- colnames(Xqc)

rafqc <- apply(Xqc,2,mean,na.rm=T)/2
isnp <- intersect(names(rafqc), sig_gwas_obs[[1]]$rsID)


for(i in 1:length(sig_gwas_obs)) sig_gwas_obs[[i]] <- sig_gwas_obs[[i]][isnp,]
for(i in 1:length(sig_gwas_latent)) sig_gwas_latent [[i]] <- sig_gwas_latent[[i]][isnp,] 



### fine-mapping
library(flashfmZero)

save.path="tmpDIR"

# for one latent factor or one blood cell trait
if(length(sig_gwas_latent)==1){
  FM_single_trait_latent <- JAMdwithGroups(sig_gwas_latent[[1]], N=N, corX, save.path, cred = 0.99,
                                                   jam.nM.iter =5, maxcv = 1, maxcv_stop = 20, 
                                                   min.mppi = 0.01, r2.minmerge = 0.8)
}

 #  for one blood cell trait
if(length(sig_gwas_obs)==1){
  FM_single_trait_raw <- JAMdwithGroups(sig_gwas_obs[[1]], N=N, corX, cred = 0.99,
                                                    jam.nM.iter =5, maxcv = 1, maxcv_stop = 20, 
                                                    min.mppi = 0.01, r2.minmerge = 0.8)
}


# # JAM on multiple blood cell traits 
 if(length(sig_gwas_obs)>1){
   FM_single_trait_raw <- multiJAMd(sig_gwas_obs,  corX, N=N, save.path=save.path,
                                    maxcv = 1, maxcv_stop = 20, jam.nM.iter =5, r2.minmerge=0.8, minsnpmppi = 0.01,
                                    #NCORES=length(gwas.list.interval.raw99) #if running on the hpc
                                    NCORES=1) 
   FM_single_trait_raw_CS99  <- multiJAMdCS(FM_single_trait_raw, cred = 0.99)
 }


# multiple latent factors flashfmZero - returns single and multi-trait results
if(length(sig_gwas_latent)>1){
  FM_flashfmzero_latent <- FLASHFMZEROwithJAMd(sig_gwas_latent, 
                                                     corX, 
                                                     N = N, 
                                                     save.path=save.path, 
                                                     TOdds = 1,
                                                     cpp = 0.99, 
                                                     #NCORES=length(gwas.list.interval.fa25), #if running on the hpc
                                                     NCORES = 1,
                                                     maxcv = 1, 
                                                     maxcv_stop = 20, 
                                                     jam.nM.iter = 5, 
                                                     r2.minmerge = 0.8, 
                                                     minsnpmppi = 0.01)
  FM_flashfmzero_latent_CS99 <- allcredsetsPP(FM_flashfmzero_latent$mpp.pp,cred=.99)
}


# save all output pbjects here in RData
