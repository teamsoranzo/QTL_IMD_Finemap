

suppressMessages(library("data.table", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("optparse", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))



opt = NULL

LOOP_ATACSeq = function(option_list)
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
  
  subset = read.table(opt$subset, sep="\t", header = T, stringsAsFactors = F)  
  
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
  
  #quit(status = 1)
  
  #### Read ATAC.scaled ----
  
  #filename_5<-paste("ATAC_scaled_relevant_tissue_info_",type,".txt", sep='')
  
  ATACSeq.scaled.m.no.NA = as.data.frame(fread(opt$ATAC_Scaled, sep="\t", stringsAsFactors = F, header = T))
   
  #filename_6<-paste("ATAC_scaled_Lineage_hierarchy",type,".txt", sep='')
  
  Lineage.hierarchy =  read.table(opt$Lineage_hierarchy, sep="\t", stringsAsFactors = F, header = T)
  
  #filename_7<-paste("ATAC_scaled_trait_table_",type,".txt", sep='')
  
  trait.table =  read.table(opt$trait_table, sep="\t", stringsAsFactors = F, header = T)
  
  cat("Check READS_\n")
  str(ATACSeq.scaled.m.no.NA)
  str(Lineage.hierarchy)
  str(trait.table)
  cat("\n")
  
  #### Read ATAC.RAW ----
  
  #filename_5<-paste("ATAC_RAW_relevant_tissue_info_",type,".txt", sep='')
  
  ATACSeq.RAW.m.no.NA = as.data.frame(fread(opt$ATAC_Raw, sep="\t", stringsAsFactors = F, header = T))
  
  cat("Check ATACSeq.RAW.m.no.NA\n")
  str(ATACSeq.RAW.m.no.NA)
  cat("\n")
  
  #### Merge kolf_ATAC with dB ----
  
  selection<-subset
  
  cat("Check selection\n")
  str(selection)
  cat("\n")

  
  #quit(status =1)
  
  ###### LOOP SUPER.GATHER through lineages and clasiff VJ ---------
  
  
  VARS<-unique(as.character(subset$VAR))
  
  
  
  SUPER.Gather<-data.frame(matrix(vector(), 0, 5,
                                  dimnames=list(c(), 
                                                c("VAR","ATAC_Cell_Type",
                                                  "ATAC_value","ATAC_value_RAW",
                                                  "Lineage"
                                                ))),
                           stringsAsFactors=F)
  
  
  
  
  for(i in 1:length(VARS))
  {
    variant<-VARS[i]
    
    cat("-------------------->variant\n")
    cat(sprintf(as.character(variant)))
    cat("\n")
    
    selection_variant<-selection[(selection$VAR==variant),]
    ATACSeq_variant<-ATACSeq.scaled.m.no.NA[which(ATACSeq.scaled.m.no.NA$VAR == variant),]
    ATACSeq_RAW_variant<-ATACSeq.RAW.m.no.NA[which(ATACSeq.RAW.m.no.NA$VAR == variant),]
    
    
    if(length(ATACSeq_variant$VAR) > 0)
    {
      
      # BIAS!!! Select ATACSeq only in the phenotypes that are relevant for the variant
      
      phenotypes<-unique(as.character(selection_variant$phenotype))
      
      cat("-------------------->phenotypes\n")
      cat(sprintf(as.character(phenotypes)))
      cat("\n")
      
     
      lineages_relevant_table<-trait.table[which(trait.table$Trait%in%phenotypes),]
      
      cat("lineages_relevant_table\n")
      str(lineages_relevant_table)
      cat("\n")
      
      # if("tcel"%in%phenotypes)
      # {
      #   quit(status = 1)
      # }
      
      lineages_relevant<-unique(lineages_relevant_table$Factor4)
      cat("-------------------->lineages_relevant\n")
      cat(sprintf(as.character(lineages_relevant)))
      cat("\n")
      
      Lineage.hierarchy_sel<-unique(Lineage.hierarchy[which(Lineage.hierarchy$Lineage%in%lineages_relevant),
                                                      -3])
      
      cat("Lineage.hierarchy_sel\n")
      str(Lineage.hierarchy_sel)
      cat("\n")
      
     
      
      
      ATACSeq_variant_subset<-unique(ATACSeq_variant[,c(which(colnames(ATACSeq_variant) == "VAR"),
                                                        which(colnames(ATACSeq_variant) == "variable"),
                                                        which(colnames(ATACSeq_variant) == "value"))])
      
      cat("ATACSeq_variant_subset\n")
      str(ATACSeq_variant_subset)
      cat("\n")
      
      ATACSeq_RAW_variant_subset<-unique(ATACSeq_RAW_variant[,c(which(colnames(ATACSeq_RAW_variant) == "VAR"),
                                                        which(colnames(ATACSeq_RAW_variant) == "variable"),
                                                        which(colnames(ATACSeq_RAW_variant) == "value"))])
      
      cat("ATACSeq_RAW_variant_subset\n")
      str(ATACSeq_RAW_variant_subset)
      cat("\n")
      
      
      ATAC_merge<-merge(ATACSeq_variant_subset,
                        ATACSeq_RAW_variant_subset,
                        by=c("VAR","variable"),
                        all=T)
      
      cat("ATAC_merge_PRE\n")
      str(ATAC_merge)
      cat("\n")
      
      
      colnames(ATAC_merge)<-c("VAR","ATAC_Cell_Type","ATAC_value","ATAC_value_RAW")
      
      cat("ATAC_merge_POST_1\n")
      str(ATAC_merge)
      cat("\n")
      
      cat("-------------------->ATAC_Cell_Type\n")
      cat(sprintf(as.character(levels(factor(ATAC_merge$ATAC_Cell_Type)))))
      cat("\n")
      cat(sprintf(as.character(summary(factor(ATAC_merge$ATAC_Cell_Type)))))
      cat("\n")
      
      colnames(Lineage.hierarchy_sel)<-c("ATAC_Cell_Type","Lineage")
      
      DEF_ATACSeq_variant<-merge(ATAC_merge,Lineage.hierarchy_sel,
                                 by="ATAC_Cell_Type")
      
      cat("DEF_ATACSeq_variant_POST_2\n")
      str(DEF_ATACSeq_variant)
      cat("\n")
      
      cat("-------------------->ATAC_Cell_Type\n")
      cat(sprintf(as.character(levels(factor(DEF_ATACSeq_variant$ATAC_Cell_Type)))))
      cat("\n")
      cat(sprintf(as.character(summary(factor(DEF_ATACSeq_variant$ATAC_Cell_Type)))))
      cat("\n")
      
      cat("-------------------->Lineage\n")
      cat(sprintf(as.character(levels(factor(DEF_ATACSeq_variant$Lineage)))))
      cat("\n")
      cat(sprintf(as.character(summary(factor(DEF_ATACSeq_variant$Lineage)))))
      cat("\n")
      
      
      DEF_ATACSeq_variant_ordered<-unique(DEF_ATACSeq_variant[,c(2,1,3,4,5)])
      
      SUPER.Gather<-rbind(SUPER.Gather,DEF_ATACSeq_variant_ordered)
      
    }else{
      
      # A<-as.data.frame(cbind(variant,"NA",0,"NA"))
      # colnames(A)<-colnames(SUPER.Gather)
      # SUPER.Gather<-rbind(SUPER.Gather,A)
    }
    
  }
  
  
  
  
  #### SAVE FILE ----
  
  setwd(out)
  
  filename_8<-paste("ATAC_RESULT_",type,".txt", sep='')
  
  cat("filename_8_9_\n")
  str(filename_8)
  cat("\n")
  
  write.table(SUPER.Gather,
              file=filename_8, sep="\t", quote=F, row.names = F)
  
   
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
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ATAC_Scaled"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ATAC_Raw"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Lineage_hierarchy"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--trait_table"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
    
  )
  
  
  
  parser = OptionParser(usage = "133__Rscript_edition_105.R
                        --dB FILE.txt
                        --subset FILE.txt
                        --ATACSeq FILE.txt
                        --Threshold_ATAC type
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  LOOP_ATACSeq(opt)
    
}


###########################################################################

system.time( main() )
