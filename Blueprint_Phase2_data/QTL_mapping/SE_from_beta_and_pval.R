## RUN: Rscript ~/Projects/Scripts/Blueprint/Re_analysis/Limix/SE_from_beta_and_pval.R < INPUTFILE > < OUTPUTFILE >


args <- commandArgs(TRUE)
file1 <- args[1] ## Input file
file2 <- args[2] ##  Output file



############## DATASETS ##############################################################

df1.use <- read.table(file1, header = TRUE, stringsAsFactors=F)


## For methylation - if it's separated into several chunks
# df1.use <- read.table(file1, header = FALSE, stringsAsFactors=F)
# colnames(df1.use) <- c("chr:pos_ref_alt", "rsid", "phenotypeID", "p.value", "beta", "Bonferroni.p", "lFDR", "alt.AF")

z.score <- sign(df1.use$beta) * qnorm( df1.use$p.value/2 , mean=0, sd=1, lower=FALSE)
df1.use$se <- (df1.use$beta/z.score)


write.table(df1.use,file=file2,row.names=F,quote=F,col.names=T,append=F)

