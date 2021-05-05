
## RUN: Rscript /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/R_Scripts/prior_prob_finemap.R


args <- commandArgs(TRUE)
K <- as.numeric(args[1]) ## No of maximum causal SNPs
m <- as.numeric(args[2]) ## No of SNPs in this region
file1 <- args[3] ## Output file for prior probability, required in FINEMAP


####################### PRIOR PROBABILITY CALCULATION ###############################

prior <- choose( m, seq_len( K ) ) * ( 1 / m ) ^ seq_len( K ) * ( ( m - 1 ) / m ) ^ ( m - seq_len( K ) )
write.table(t(prior),file=file1,row.names=F,quote=F,col.names=F,append=F, sep=" ") ## NOTE the prior is transposed to get into a single line