#!/usr/bin/bash Rscript

library(qvalue)
library(graphics)


args <- commandArgs(TRUE)



file=args[1]
name=args[2]


inQTL<-read.csv(file,sep="\t",header=TRUE)
Eigen <-inQTL$BF
Eigen[Eigen>1]<-1
gFDR <- qvalue(p=Eigen)
outQTL<-data.frame(inQTL,gFDR$qvalues)

# headers<-c("chr:pos_ref_alt","rsid","phenotypeID","p.value","beta","Bonferroni.p","lFDR", "alt.AF","gFDR")
# oname=paste("gFDR.added.",name,sep="")
# colnames(outQTL)<-headers

write.table(outQTL,file=name,quote=FALSE,row.names=FALSE,col.names=TRUE, sep = "\t")



## plotting gFDR vs pvalue
#tab<-read.csv(file,header=FALSE,sep=" ")
# headers<-c("pos","rsid","phenotypeID","pv","beta","Bonferroni","lFDR", "AF","qv_all")
# colnames(outQTL)<-headers

tab<-outQTL
pv<-tab$p.value[tab$gFDR.qvalues<=0.05]
max.pv<-max(pv)
max.pv

#linear fit
y<-as.vector(-log10(tab$gFDR.qvalues))
x<-as.vector(-log10(tab$gFDR.qvalues))

mod<-lm(y~x)
mod
mod.pv<-(-log10(0.05)-as.numeric(mod$coefficients[1]))/as.numeric(mod$coefficients[2])
predict.pv<-10^(-mod.pv)
#as.numeric(mod$coefficients[1])
#mod$coefficients[2]

#plot
filepdf<-paste(name,".pdf",sep="")
pdf(filepdf)
plot(-log10(tab$p.value),-log10(tab$gFDR.qvalues),main=paste("-log10: max pv data",round(-log10(max.pv),2),"predicted max pv",round(mod.pv,2),sep=" "),xlab="-log10 Pval",ylab="-log10 gFDR")
abline(h=-log10(0.05),col="red")
abline(v=-log10(max.pv))
abline(v=mod.pv,col="blue")
abline(lm(y~x),col="blue")
legend(40,20,"5% gFDR",text.col="red",bty="n")
legend(40,15,"lm",text.col="blue",bty="n")
dev.off()
