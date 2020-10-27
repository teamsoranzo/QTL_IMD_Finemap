
suppressMessages(library("data.table", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("optparse", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("splitstackshape", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
library("desiR", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/")
library("plyr", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/")


opt = NULL



Rank_ATAC_Seq_signal = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
 
  #### READ and transform subset ----
  
  subset = read.table(opt$subset, sep="\t", stringsAsFactors = F, header = T)
  
  colnames(subset)[which(colnames(subset) == "MarkerName")]<-"VAR"
  subset$VAR<-paste('chr',subset$VAR,sep='')
  subset$VAR<-gsub(':','_',subset$VAR)
  colnames(subset)[which(colnames(subset) == "Celltype")]<-"phenotype"
  
  cat("-------------------->phenotypes\n")
  cat(sprintf(as.character(levels(factor(subset$phenotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(factor(subset$phenotype)))))
  cat("\n")
  
  cat("SUBSET_\n")
  cat(str(subset))
  cat("\n")
  
  cat("SUBSET_\n")
  cat(str(subset))
  cat("\n")
  
  
  
  #  quit(status=1)
  
  ### Read the ATAC_RESULT file ----
  
 
  ATAC_original = read.table(file=opt$ATAC_RESULT, sep="\t", stringsAsFactors = F, header = T)
  
  colnames(ATAC_original)[which(colnames(ATAC_original) == "Factor4")]<-"Lineage"
  
  cat("ATAC_original\n")
  str(ATAC_original)
  cat("\n")
  
  ATAC_original$Lineage<-gsub("lineage.+","lineage",ATAC_original$Lineage,perl=T)
  
  ATAC_original<-unique(ATAC_original)
 
  if(length(ATAC_original$VAR) > 0)
  {
  
  #### read parameters ----
  
  # parameters<-c(3,18,0.05,0.8,1,3)
  
  parameters<-unlist(strsplit(opt$parameters, split=","))
  
  cat("parameters\n")
  cat(sprintf(as.character(parameters)))
  cat("\n")
  
  a<-as.numeric(parameters[1])
  b<-as.numeric(parameters[2])
  
  cat("check.a\n")
  cat(sprintf(as.character(a)))
  cat("\n")

  cat("check.b\n")
  cat(sprintf(as.character(b)))
  cat("\n")
  
  #### READ and transform out ----
  
  Freq_Threshold = as.numeric(opt$Freq_Threshold)
  
  cat("Freq_Threshold_\n")
  cat(sprintf(as.character(Freq_Threshold)))
  cat("\n")
  
  
  #### MAX_ATAC_original ----
  
  ATAC_original.dt<-data.table(ATAC_original, key="VAR")
  
  cat("ATAC_original.dt\n")
  cat(str(ATAC_original.dt))
  cat("\n")
  
  MAX_ATAC_original<-as.data.frame(ATAC_original.dt[, .SD[which.max(ATAC_value)], by=VAR])
  
  cat("MAX_ATAC_original\n")
  cat(str(MAX_ATAC_original))
  cat("\n")
  
  #### Frequencies ----
  
  ATAC_original_Thresholded<-ATAC_original[which(ATAC_original$ATAC_value >= Freq_Threshold),]
  
  cat("ATAC_original_Thresholded\n")
  cat(str(ATAC_original_Thresholded))
  cat("\n")
  
  #### IF condition for chunks ----
  
  if(length(ATAC_original_Thresholded$VAR) > 0)
  {
    Frequencies<-as.data.frame(table(ATAC_original_Thresholded$VAR))
    colnames(Frequencies)<-c("VAR","Freq")
    
    Frequencies$VAR<-as.character(Frequencies$VAR)
    
    cat("------------------->Frequencies_POST\n")
    str(Frequencies)
    cat("\n")
    
    
    #### desiR on frequency ----
    setwd(out)
    
    pdf(file="desIR_ATAC_Frequencies.pdf")
    hist(Frequencies$Freq, breaks=50, col="grey", border="white", main="",
         xlab="ATAC-Seq in a Relevant Cell Type")
    des.line(Frequencies$Freq, "d.high", des.args=c(cut1=a, cut2=b, scale=0.5))
    dev.off()
    
    Frequencies$Frequency_weight <- d.high(Frequencies$Freq, cut1=a, cut2=b, scale=0.5)
    
    #### desiR on MAX_ATAC_original ----
    
    a<-as.numeric(parameters[3])
    b<-as.numeric(parameters[4])
    
    cat("check.a\n")
    cat(sprintf(as.character(a)))
    cat("\n")
    
    cat("check.b\n")
    cat(sprintf(as.character(b)))
    cat("\n")
    
    pdf(file="desIR_ATAC_value.pdf")
    hist(MAX_ATAC_original$ATAC_value, breaks=50, col="grey", border="white", main="",
         xlab=paste("ATAC_value"))
    des.line(MAX_ATAC_original$ATAC_value, "d.high", des.args=c(cut1=a, cut2=b, scale=0.5))
    dev.off()
    
    
    
    MAX_ATAC_original$ATAC_MAX_weight<-d.high(MAX_ATAC_original$ATAC_value, cut1=a, cut2=b, scale=0.5)
    
    
    #### MERGE ----
    
    DEF<-merge(Frequencies,
               MAX_ATAC_original,
               by="VAR",
               all=T)
    
    
    cat("------------------->DEF_1\n")
    str(DEF)
    cat("\n")
    
    DEF$Frequency_weight[is.na(DEF$Frequency_weight)]<-0.1
    
    #quit(status = 1)
    
    #### Overall desirability ---- 
    
    a<-as.numeric(parameters[5])
    b<-as.numeric(parameters[6])
    
    
    cat("check.a\n")
    cat(sprintf(as.character(a)))
    cat("\n")
    
    cat("check.b\n")
    cat(sprintf(as.character(b)))
    cat("\n")
    
    DEF$Overall_weight <- d.overall(DEF$Frequency_weight, DEF$ATAC_MAX_weight, 
                                    weights=c(a,b))
    
    par(las=1)
    pdf(file= "desIR_overall.pdf")
    plot(rev(sort(DEF$Overall_weight)), type="l", xlab="Rank", ylab="Overall Desirability")
    dev.off()
    
    #### Discretize ----
    
    
    # A<-c(summary(DEF$Overall_weight)[1],
    #      summary(DEF$Overall_weight)[2],
    #      summary(DEF$Overall_weight)[3],
    #      summary(DEF$Overall_weight)[5],
    #      (summary(DEF$Overall_weight)[6] +0.1))
    # 
    # 
    # DEF$Rnk_discretized<- cut(DEF$Overall_weight,breaks = A,right = FALSE)
    # 
    # 
    # 
    # cat("Rnk_discretized\n")
    # cat(sprintf(as.character(levels(as.factor(DEF$Rnk_discretized)))))
    # cat("\n")
    # cat(sprintf(as.character(summary(as.factor(DEF$Rnk_discretized)))))
    # cat("\n")
    
    
    DEF_subset<-DEF[,c(which(colnames(DEF)== "Overall_weight"),
                       which(colnames(DEF)== "VAR"))]
    
    
    DEF_subset<-DEF_subset[order(DEF_subset$Overall_weight, decreasing = F),]
    
    
    #ordered_levels_ATAC<-unique(DEF_subset$Rnk_discretized)
    
    MERGE_1<-merge(DEF_subset,
                   ATAC_original,
                   by="VAR",
                   all=T)
    
    MERGE_1_subset<-unique(MERGE_1[,c(which(colnames(MERGE_1)== "VAR"),
                                      which(colnames(MERGE_1)=="Overall_weight"),
                                      which(colnames(MERGE_1)== "Lineage"))])
    
    #### SAVE FILE ----
    
    setwd(out)
    
    
    
    filename18<-paste("ATAC_RESULT_Ranked_",type,".txt", sep='')
    
    write.table(MERGE_1_subset,file=filename18,sep="\t",row.names = F,quote = F)
    
  }# length(ATAC_original_Thresholded$VAR)
}# length(ATAC_original$VAR)
  
}


printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}

#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--subset"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ATAC_RESULT"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--parameters"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Freq_Threshold"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
    
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --ATAC_original FILE.txt
                        --GeneEXP FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  
   
   Rank_ATAC_Seq_signal(opt)
  
  
}


###########################################################################

system.time( main() )
