## RUN: Rscript /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/R_Scripts/Z_score_from_beta_and_pval.R < INPUTFILE > < sample size > < OUTPUTFILE >


args <- commandArgs(TRUE)
file1 <- args[1] ## Input file
file2 <- args[2] ## Output file for z scores; required in FINEMAP and CAVIARBF



############## DATASETS ##############################################################
df1.use <- read.table(file1, header = TRUE,stringsAsFactors=F)


df1.use$z.score <- sign(df1.use$beta) * qnorm( df1.use$pval/2 , mean=0, sd=1, lower=FALSE)


out.file<-unique(df1.use[,c("RSid","alt","ref","alt.AF","beta","StdErr","pval","N","z.score","StdID","GeneID","chr","Pos")])

write.table(out.file,file=file2,row.names=F,quote=F,col.names=T,append=F,sep="\t")
