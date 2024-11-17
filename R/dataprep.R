

#' @title Flip GWAS variants to align with SNP correlation matrix (LD) from reference panel or in-sample 
#' @param gwas GWAS data.frame with the following columns (any order, could have others but require these): "rsID" (variant ID), "EA" (effect allele), "NEA" (non-effect allele),
#' "beta" (effect size), "EAF" (effect allele frequency)
#' @param RPinfo reference panel (or in-sample) details data.frame with the following columns (any order): "rsID" (variant ID), "allele1", "allele2" 
#' (this function will align the gwas to allele1, but either allele1 or allele2 may be used -  need consistency for correlation signs)
#' @param details default value is FALSE to return only the aligned GWAS data.frame; TRUE will also return names of excluded snps
#' @return if details=FALSE, return only data.frame of gwas aligned to reference panel (RP) or in-sample LD
#' if details=TRUE: a list with three components: 
#' gwasA = gwas aligned to reference panel (RP) or in-sample LD;
#' excluded = gwas variants that are excluded (not in RP, incompatible alleles (even if flip));
#' ind_excl = indices of removed rows from imput gwas
#' @author Jenn Asimit
#' @export
alignGWAS <- function(gwas,RPinfo,details=FALSE) {

 snpkeep <- intersect(gwas$rsID,RPinfo$rsID)
 if(length(snpkeep)==0) error("There is no overlap between the GWAS and reference panel SNP names. Check input.")
 rownames(gwas) <- gwas$rsID
 rownames(RPinfo) <- RPinfo$rsID
 gwas <- gwas[snpkeep,]
 RPinfo <- RPinfo[snpkeep,]
  
 
 flip1 <- which(gwas$EA != RPinfo$allele1 | gwas$NEA != RPinfo$allele2) # check for discrepancies
 check <- excluded <- indrm <- c()
 
 if(length(flip1)>0) check <- which( gwas$EA[flip1] == RPinfo$allele2[flip1] & gwas$NEA[flip1] == RPinfo$allele1[flip1] )  # snps to flip

 if(length(check)>0){
   gwas$NEA[flip1[check]] <- RPinfo$allele2[flip1[check]]
   gwas$EA[flip1[check]] <- RPinfo$allele1[flip1[check]]
   gwas$beta[flip1[check]] <- -gwas$beta[flip1[check]]
   gwas$EAF[flip1[check]] <- 1-gwas$EAF[flip1[check]]

   indrm <-  setdiff(flip1,flip1[check]) # rm these snps from both datasets since allele codings still don't agree if flipped 
   if(length(indrm)>0) {
    excluded <- gwas[indrm,] 
    gwas <- gwas[-indrm,]      
   }
   
 }

if(details) {
out <- list(gwasA=gwas,excluded=excluded,ind_exc=indrm)
} else { out <- gwas}

return(out)

}



LDqc <- function(g,theta=0.1,BestGuess=TRUE) {
  ind0 <- which(g<=theta)
  ind1 <- which(g>=1-theta & g<=1+theta)
  ind2 <- which(g>=2-theta)
  bg <- rep(NA,length(g))
  if(BestGuess) {
    if(length(ind0)>0) bg[ind0] <- 0
    if(length(ind1)>0) bg[ind1] <- 1
    if(length(ind2)>0) bg[ind2] <- 2 
  } else{
    if(length(ind0)>0) bg[ind0] <- g[ind0]
    if(length(ind1)>0) bg[ind1] <- g[ind1]
    if(length(ind2)>0) bg[ind2] <- g[ind2]   
  }  
  return(bg)
}

