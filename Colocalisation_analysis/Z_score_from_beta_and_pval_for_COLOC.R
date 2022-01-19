## RUN: Rscript /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/R_Scripts/Z_score_from_beta_and_pval_for_COLOC.R < INPUTFILE > < OUTPUTFILE >


args <- commandArgs(TRUE)
file1 <- args[1] ## Input file
file2 <- args[2] ## Output file for z scores and se; required for gwas-pw (For Colocalization)



############## DATASETS ##############################################################
df1.use <- read.table(file1, header = FALSE, stringsAsFactors=F)
colnames(df1.use) <- c("SNP", "pval", "beta")
df1.use$pval<- as.numeric(df1.use$pval)

df1.use$z.score <- sign(df1.use$beta) * qnorm(df1.use$pval/2 , mean=0, sd=1, lower=FALSE)
df1.use$se <- (df1.use$beta/df1.use$z.score)

out.file<-unique(df1.use[,c("SNP", "pval","beta","se","z.score")])

write.table(out.file,file=file2,row.names=F,quote=F,col.names=F,append=F,sep="\t")

