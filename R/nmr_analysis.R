
library(tidyverse)
library(psych)

OUTPUT_DIR="/home/wja24/data/jenn_flashfmm_output/"
R_DIR="/home/wja24/repos/analysis_jenn_flashfmm.git/R"
source(sprintf("%s/latentGWAS.R", R_DIR))

lipid_phenos<-read_csv("/home/wja24/rds/rds-jmmh2-post_qc_data/interval/phenotype/nightingale_metabolomics/gwasqc/nmr_qcgwas.csv")

Y_start<-lipid_phenos[,-1]
N <- nrow(Y_start)

#check the number of missing values
count_non_missing_values_start <- unlist(as.data.frame(colSums(!is.na(Y_start))))
missing_data_fraction <- 1-count_non_missing_values_start/N


# trait correlation
corY_start <- cor(Y_start, use="pairwise.complete.obs", method='pearson')


# Build the adjacency matrix
adj_mat <- abs(corY_start) > 0.99
diag(adj_mat) <- FALSE  # Remove self-loops


# Initialize the list of variables to keep
variables_to_trim <- NULL

# Iteratively remove variables
while (any(adj_mat)) {
  # Compute the degree (number of low-correlation connections) for each variable
  deg <- rowSums(adj_mat)
  
  # Compute the score: higher degree and more missing data increase the score
  score <- deg * (1 + missing_data_fraction)
  
  remove_index<-which.max(score)
  variables_to_trim<-c(variables_to_trim, remove_index)

  # Update the adjacency matrix by removing the corresponding row and column
  adj_mat[remove_index, ] = FALSE
  adj_mat[,remove_index] = FALSE
 }

# Subsetdata frame to include only the variables to keep
Y<-Y_start[,! 1:dim(Y_start)[2]  %in% variables_to_trim]
count_non_missing_values <- unlist(as.data.frame(colSums(!is.na(Y))))
P <- ncol(Y)


#check pairwise missing values
count_non_missing_values_pairwise = matrix(0,P,P)
colnames(count_non_missing_values_pairwise) <- rownames(count_non_missing_values_pairwise) <-  colnames(Y)
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


median_nonmissing2 <- median(count_non_missing_values_pairwise)
#fa.parallel(covY, fm="ml", n.obs=median_nonmissing2,fa="fa",show.legend=F,main="") # use median non-missing pairwise


## set nfactors to the number of latent factors suggested above
nfactors <- 21 

pdf(sprintf("%s/lipid/factor_anal/scree.pdf",OUTPUT_DIR)) # jpeg("scree.jpg") tiff("scree.tiff") 
fa.parallel(covY, fm="ml", n.obs=median_nonmissing2,fa="fa",show.legend=F,main="")
abline(v=nfactors,col="red")
dev.off()

write.csv(covY, file=sprintf("%s/lipid/factor_anal/Trait_cov.csv", OUTPUT_DIR))


# using number of latent factors from above, run factor analysis using covariance matrix
faY <- fa(r=covY, nfactors=nfactors, fm="ml", rotate="varimax")
faYloading <- as.matrix(faY$loadings)


####rescaled factors to get contributions 
latent_terms <- factor_contributions(faYloading) # output is factor loadings and re-scaled factor loadings (contributions) matrices with observed traits ordered by maximum contributing latent factor
rescaledFAloading <- latent_terms[[1]]
FAloading <- latent_terms[[2]]

write.csv(rescaledFAloading, file=sprintf("%s/Re-scaled-Factor_loading_matrix.csv",OUTPUT_DIR))
write.csv(FAloading, file=sprintf("%s/Factor_loading_matrix.csv",OUTPUT_DIR))
write.csv(covY[rownames(FAloading),rownames(FAloading)], file=sprintf("%s/Trait_cov.csv",OUTPUT_DIR))




