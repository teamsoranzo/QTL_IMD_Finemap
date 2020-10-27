

suppressMessages(library("data.table", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("Biostrings", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("optparse", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("splitstackshape", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("data.table", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))


opt = NULL



scaled_VJ_per_column = function(option_list)
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
  
  #### READ ATACSeq_1 ----
  
  ATACSeq_1 = as.data.frame(fread(opt$ATACSeq_1, sep="\t", stringsAsFactors = F, header = T))
  
  indx.deplete<-which(colnames(ATACSeq_1) == "rsid")
 
  ATACSeq_1<-ATACSeq_1[,-indx.deplete] 
  
  cat("ATACSeq_1_fields_\n")
  str(ATACSeq_1)
  cat("\n")
  
  #### READ ATACSeq_2 ----
  
  ATACSeq_2 = as.data.frame(fread(opt$ATACSeq_2, sep="\t", stringsAsFactors = F, header = T))
  
  
  cat("ATACSeq_2_fields_\n")
  str(ATACSeq_2)
  cat("\n")
  
  #### combine ATAC files ----
  
  ATACSeq = rbind(ATACSeq_1,ATACSeq_2)
  
  
  cat("ATACSeq_fields_\n")
  str(ATACSeq)
  cat("\n")
  
  
  #### VJ cells to lineages ----
  
  erythroid <- c("mchc","ret_p","ret","mrv","mscv","mch","mcv","rdw_cv","MicroR",
                 "hgb","rbc","hct","irf","hlr_p","hlr")
  erythroid_lineage <- c("HSC", "MPP", "CMP", "MEP", "Ery")
  
  mega <- c("plt","pct","mpv","pdw","H???IPF","IPF%","IPF","P???LCR")
  mega_lineage <- c("HSC", "MPP", "CMP", "MEP", "Mega")
  
  gran_mono <- c("neut","mono_p","mono","neut_p","wbc","baso","baso_p","eo_p","eo")
  gran_mono_lineage <- c("HSC", "MPP","CMP", "GMP.A", "Mono")
  
  lymph <- c("lymph_p","lymph")
  lymph_lineage.NK <- c("HSC", "MPP", "LMPP", "CLP","NK")
  
  lymph <- c("lymph_p","lymph")
  lymph_lineage.B <- c("HSC", "MPP", "LMPP", "CLP","B")
  
  lymph <- c("lymph_p","lymph")
  lymph_lineage.CD8 <- c("HSC", "MPP", "LMPP", "CLP", "CD8")
  
  lymph <- c("lymph_p","lymph")
  lymph_lineage.CD4 <- c("HSC", "MPP", "LMPP", "CLP", "CD4")
  
  
  ############## Normalize the ATACSeq data ----------
  # Maximum of openness per cell type, that's a 1 per cell type. Then
  # scale down the values in the same cell type
  
  ATACSeq.dt<-data.table(ATACSeq, key=c("MarkerName","seqnames","start","end","width","strand"))
  
  cat("ATACSeq_fields_2_\n")
  str(ATACSeq.dt)
  cat("\n")
  
  #quit(status = 1)
  
  ATACSeq.scaled<-ATACSeq.dt[,.(B/max(ATACSeq.dt$B),
                                                CD4/max(ATACSeq.dt$CD4),CD8/max(ATACSeq.dt$CD8),CLP/max(ATACSeq.dt$CLP),
                                                CMP/max(ATACSeq.dt$CMP),Ery/max(ATACSeq.dt$Ery),GMP.A/max(ATACSeq.dt$GMP.A),
                                                GMP.B/max(ATACSeq.dt$GMP.B),GMP.C/max(ATACSeq.dt$GMP.C),HSC/max(ATACSeq.dt$HSC),
                                                LMPP/max(ATACSeq.dt$LMPP),mDC/max(ATACSeq.dt$mDC),Mega/max(ATACSeq.dt$Mega),
                                                MEP/max(ATACSeq.dt$MEP),Mono/max(ATACSeq.dt$Mono),MPP/max(ATACSeq.dt$MPP),
                                                NK/max(ATACSeq.dt$NK),pDC/max(ATACSeq.dt$pDC)),.(MarkerName)]
  
  cat("ATACSeq_fields_3_\n")
  str(ATACSeq.scaled)
  cat("\n")
  
  ATACSeq.scaled<-as.data.frame(ATACSeq.scaled)
  
  colnames(ATACSeq.scaled)[2:length(colnames(ATACSeq.scaled))]<-c(
    "B",
    "CD4","CD8","CLP",
    "CMP","Ery","GMP.A",
    "GMP.B","GMP.C","HSC",
    "LMPP","mDC","Mega",
    "MEP","Mono","MPP",
    "NK","pDC" 
  )
  
  cat("ATACSeq_fields_4_\n")
  str(ATACSeq.scaled)
  cat("\n")
  
  
 
  ########## Lineage.hierarchy -----------
  
  Lineage.hierarchy<- data.frame(matrix(vector(), 0, 3,
                                        dimnames=list(c(), 
                                                      c("Lineage","CellType","Factor5"
                                                      ))),
                                 stringsAsFactors=T)
  
  
  
  
  ###
  
  
  ATACSeq.scaled.m<-melt(ATACSeq.scaled)
  
  ATACSeq.scaled.m$Factor4<-rep("NA",length(ATACSeq.scaled.m$var))
  
  cat("ATACSeq_fields_5_\n")
  str(ATACSeq.scaled.m)
  cat("\n")
  
  toMatch1 <- mega_lineage
  length(toMatch1)
  ATACSeq.scaled.m$Factor4[which(ATACSeq.scaled.m$variable%in%toMatch1)]<-"mega_lineage"
  
  temp.table<-data.frame(cbind(toMatch1,rep("mega_lineage",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","CMP","MEP","Mega"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
    
  toMatch1 <- erythroid_lineage
  length(toMatch1)
  ATACSeq.scaled.m$Factor4[which(ATACSeq.scaled.m$variable%in%toMatch1)]<-"erythroid_lineage"
  
  temp.table<-data.frame(cbind(toMatch1,rep("erythroid_lineage",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","CMP","MEP","Ery"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  toMatch1 <- gran_mono_lineage
  length(toMatch1)
  ATACSeq.scaled.m$Factor4[which(ATACSeq.scaled.m$variable%in%toMatch1)]<-"gran_mono_lineage"
  
  temp.table<-data.frame(cbind(toMatch1,rep("gran_mono_lineage",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","CMP","GMP.A","Mono"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  ###
  
  toMatch1 <- lymph_lineage.CD4
  length(toMatch1)
  ATACSeq.scaled.m$Factor4[which(ATACSeq.scaled.m$variable%in%toMatch1)]<-"lymph_lineage"
  
  
  temp.table<-data.frame(cbind(toMatch1,rep("lymph_lineage.CD4",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","LMPP","CLP","CD4"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  ###
  
  toMatch1 <- lymph_lineage.CD8
  length(toMatch1)
  ATACSeq.scaled.m$Factor4[which(ATACSeq.scaled.m$variable%in%toMatch1)]<-"lymph_lineage"
  
  
  temp.table<-data.frame(cbind(toMatch1,rep("lymph_lineage.CD8",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","LMPP","CLP","CD8"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  ###
  
  toMatch1 <- lymph_lineage.B
  length(toMatch1)
  ATACSeq.scaled.m$Factor4[which(ATACSeq.scaled.m$variable%in%toMatch1)]<-"lymph_lineage"
  
  
  temp.table<-data.frame(cbind(toMatch1,rep("lymph_lineage.B",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","LMPP","CLP","B"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  ###
  
  toMatch1 <- lymph_lineage.NK
  length(toMatch1)
  ATACSeq.scaled.m$Factor4[which(ATACSeq.scaled.m$variable%in%toMatch1)]<-"lymph_lineage"
  
  
  temp.table<-data.frame(cbind(toMatch1,rep("lymph_lineage.NK",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","LMPP","CLP","NK"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  
  
  # na cell types; 
  
  na.values<-ATACSeq.scaled.m[which(ATACSeq.scaled.m$Factor4 == "NA"),]
  
  summary(droplevels(as.factor(na.values$variable))) # checked GMP.B GMP.C   mDC   pDC
  
  
  
  
  
  #### BIAS!!! ----
  
    
  # Not assigned "GMP.B" "GMP.C" "mDC"   "pDC"
  
  ATACSeq.scaled.m.no.NA<-ATACSeq.scaled.m[which(ATACSeq.scaled.m$Factor4 != "NA"),]
  
  cat("ATACSeq.scaled.m.no.NA_8_\n")
  str(ATACSeq.scaled.m.no.NA)
  cat("\n")
  
  #### Transform MarkerName -> VAR ----
  
  colnames(ATACSeq.scaled.m.no.NA)[which(colnames(ATACSeq.scaled.m.no.NA) == 'MarkerName')]<-"VAR"
  
  cat("ATACSeq.scaled.m.no.NA_9_\n")
  str(ATACSeq.scaled.m.no.NA)
  cat("\n")
  #quit(status = 1)
  
  #### SAVING ATACSeq ----
  
  setwd(out)
  
  filename_5<-paste("ATAC_scaled_relevant_tissue_info_",type,".txt", sep='')
  
  cat("filename_5_9_\n")
  str(filename_5)
  cat("\n")
  
  write.table(ATACSeq.scaled.m.no.NA,
              file=filename_5, sep="\t", quote=F, row.names = F)
  
  filename_6<-paste("ATAC_scaled_Lineage_hierarchy",type,".txt", sep='')

  write.table(Lineage.hierarchy,
              file=filename_6, sep="\t", quote=F, row.names = F)


  
  
  
}

RAW_VJ_per_column = function(option_list)
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
  
  #### READ ATACSeq_1 ----
  
  ATACSeq_1 = as.data.frame(fread(opt$ATACSeq_1, sep="\t", stringsAsFactors = F, header = T))
  
  indx.deplete<-which(colnames(ATACSeq_1) == "rsid")
  
  ATACSeq_1<-ATACSeq_1[,-indx.deplete] 
  
  cat("ATACSeq_1_fields_\n")
  str(ATACSeq_1)
  cat("\n")
  
  #### READ ATACSeq_2 ----
  
  ATACSeq_2 = as.data.frame(fread(opt$ATACSeq_2, sep="\t", stringsAsFactors = F, header = T))
  
  
  cat("ATACSeq_2_fields_\n")
  str(ATACSeq_2)
  cat("\n")
  
  #### combine ATAC files ----
  
  ATACSeq = rbind(ATACSeq_1,ATACSeq_2)
  
  
  cat("ATACSeq_fields_\n")
  str(ATACSeq)
  cat("\n")
  
  
  
  #### VJ cells to lineages ----
  
  erythroid <- c("mchc","ret_p","ret","mrv","mscv","mch","mcv","rdw_cv","MicroR",
                 "hgb","rbc","hct","irf","hlr_p","hlr")
  erythroid_lineage <- c("HSC", "MPP", "CMP", "MEP", "Ery")
  
  mega <- c("plt","pct","mpv","pdw","H???IPF","IPF%","IPF","P???LCR")
  mega_lineage <- c("HSC", "MPP", "CMP", "MEP", "Mega")
  
  gran_mono <- c("neut","mono_p","mono","neut_p","wbc","baso","baso_p","eo_p","eo")
  gran_mono_lineage <- c("HSC", "MPP","CMP", "GMP.A", "Mono")
  
  lymph <- c("lymph_p","lymph")
  lymph_lineage.NK <- c("HSC", "MPP", "LMPP", "CLP","NK")
  
  lymph <- c("lymph_p","lymph")
  lymph_lineage.B <- c("HSC", "MPP", "LMPP", "CLP","B")
  
  lymph <- c("lymph_p","lymph")
  lymph_lineage.CD8 <- c("HSC", "MPP", "LMPP", "CLP", "CD8")
  
  lymph <- c("lymph_p","lymph")
  lymph_lineage.CD4 <- c("HSC", "MPP", "LMPP", "CLP", "CD4")
  
  
  ############## Normalize the ATACSeq data ----------
  # Maximum of openness per cell type, that's a 1 per cell type. Then
  # scale down the values in the same cell type
  
  ATACSeq.dt<-data.table(ATACSeq, key=c("MarkerName","seqnames","start","end","width","strand"))
  
  cat("ATACSeq_fields.dt\n")
  str(ATACSeq.dt)
  cat("\n")
  
  #quit(status = 1)
  
  ATACSeq.RAW<-ATACSeq.dt[,.(B,CD4,CD8,CLP,
                                CMP,Ery,GMP.A,
                                GMP.B,GMP.C,HSC,
                                LMPP,mDC,Mega,
                                MEP,Mono,MPP,
                                NK,pDC), .(MarkerName)]
  
  cat("ATACSeq_fields_3_\n")
  str(ATACSeq.RAW)
  cat("\n")
  
  ATACSeq.RAW<-as.data.frame(ATACSeq.RAW)
  
  colnames(ATACSeq.RAW)[2:length(colnames(ATACSeq.RAW))]<-c(
    "B",
    "CD4","CD8","CLP",
    "CMP","Ery","GMP.A",
    "GMP.B","GMP.C","HSC",
    "LMPP","mDC","Mega",
    "MEP","Mono","MPP",
    "NK","pDC" 
  )
  
  cat("ATACSeq_fields_4_\n")
  str(ATACSeq.RAW)
  cat("\n")
  
  
  
  ########## Lineage.hierarchy -----------
  
  Lineage.hierarchy<- data.frame(matrix(vector(), 0, 3,
                                        dimnames=list(c(), 
                                                      c("Lineage","CellType","Factor5"
                                                      ))),
                                 stringsAsFactors=T)
  
  
  
  
  ###
  
  
  ATACSeq.RAW.m<-melt(ATACSeq.RAW)
  
  ATACSeq.RAW.m$Factor4<-rep("NA",length(ATACSeq.RAW.m$var))
  
  cat("ATACSeq_fields_5_\n")
  str(ATACSeq.RAW.m)
  cat("\n")
  
  toMatch1 <- mega_lineage
  length(toMatch1)
  ATACSeq.RAW.m$Factor4[which(ATACSeq.RAW.m$variable%in%toMatch1)]<-"mega_lineage"
  
  temp.table<-data.frame(cbind(toMatch1,rep("mega_lineage",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","CMP","MEP","Mega"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  
  toMatch1 <- erythroid_lineage
  length(toMatch1)
  ATACSeq.RAW.m$Factor4[which(ATACSeq.RAW.m$variable%in%toMatch1)]<-"erythroid_lineage"
  
  temp.table<-data.frame(cbind(toMatch1,rep("erythroid_lineage",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","CMP","MEP","Ery"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  toMatch1 <- gran_mono_lineage
  length(toMatch1)
  ATACSeq.RAW.m$Factor4[which(ATACSeq.RAW.m$variable%in%toMatch1)]<-"gran_mono_lineage"
  
  temp.table<-data.frame(cbind(toMatch1,rep("gran_mono_lineage",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","CMP","GMP.A","Mono"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  ###
  
  toMatch1 <- lymph_lineage.CD4
  length(toMatch1)
  ATACSeq.RAW.m$Factor4[which(ATACSeq.RAW.m$variable%in%toMatch1)]<-"lymph_lineage"
  
  
  temp.table<-data.frame(cbind(toMatch1,rep("lymph_lineage.CD4",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","LMPP","CLP","CD4"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  ###
  
  toMatch1 <- lymph_lineage.CD8
  length(toMatch1)
  ATACSeq.RAW.m$Factor4[which(ATACSeq.RAW.m$variable%in%toMatch1)]<-"lymph_lineage"
  
  
  temp.table<-data.frame(cbind(toMatch1,rep("lymph_lineage.CD8",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","LMPP","CLP","CD8"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  ###
  
  toMatch1 <- lymph_lineage.B
  length(toMatch1)
  ATACSeq.RAW.m$Factor4[which(ATACSeq.RAW.m$variable%in%toMatch1)]<-"lymph_lineage"
  
  
  temp.table<-data.frame(cbind(toMatch1,rep("lymph_lineage.B",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","LMPP","CLP","B"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  ###
  
  toMatch1 <- lymph_lineage.NK
  length(toMatch1)
  ATACSeq.RAW.m$Factor4[which(ATACSeq.RAW.m$variable%in%toMatch1)]<-"lymph_lineage"
  
  
  temp.table<-data.frame(cbind(toMatch1,rep("lymph_lineage.NK",length(toMatch1))), stringsAsFactors=T)
  colnames(temp.table)<-c("CellType","Lineage")                  
  temp.table$Factor5<-factor(temp.table$CellType,
                             levels=c("HSC","MPP","LMPP","CLP","NK"),
                             ordered =T)             
  
  
  Lineage.hierarchy<-rbind(Lineage.hierarchy,temp.table)
  
  
  
  # na cell types; 
  
  na.values<-ATACSeq.RAW.m[which(ATACSeq.RAW.m$Factor4 == "NA"),]
  
  summary(droplevels(as.factor(na.values$variable))) # checked GMP.B GMP.C   mDC   pDC
  
  
 
  cat("Lineage.hierarchy_6_\n")
  str(Lineage.hierarchy)
  cat("\n")
 
  
  
  
  #### BIAS!!! ----
  
  
  # Not assigned "GMP.B" "GMP.C" "mDC"   "pDC"
  
  ATACSeq.RAW.m.no.NA<-ATACSeq.RAW.m[which(ATACSeq.RAW.m$Factor4 != "NA"),]
  
  cat("ATACSeq.RAW.m.no.NA_8_\n")
  str(ATACSeq.RAW.m.no.NA)
  cat("\n")
  
  #### Transform MarkerName -> VAR ----
  
  colnames(ATACSeq.RAW.m.no.NA)[which(colnames(ATACSeq.RAW.m.no.NA) == 'MarkerName')]<-"VAR"
  
  cat("ATACSeq.RAW.m.no.NA_9_\n")
  str(ATACSeq.RAW.m.no.NA)
  cat("\n")
  #quit(status = 1)
  
  #### SAVING ATACSeq ----
  
  setwd(out)
  
  filename_5<-paste("ATAC_RAW_",type,".txt", sep='')
  
  cat("filename_5_9_\n")
  str(filename_5)
  cat("\n")
  
  write.table(ATACSeq.RAW.m.no.NA,
              file=filename_5, sep="\t", quote=F, row.names = F)
  
  
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
    make_option(c("--ATACSeq_1"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ATACSeq_2"), type="character", default=NULL, 
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
                        --subset FILE.txt
                        --Equiv_table FILE.txt
                        --PCHIC FILE.txt 
                        --Relev_Tissue_EXP FILE.txt
                        --ENSEMBL_bedfile FILE.txt
                        --PCHIC_Threshold FILE.txt
                        --ATACSeq FILE.txt
                        --chrst_states FILE.txt
                        --chrst_legend FILE.txt
                        --chrst_prob_threshold FILE.txt
                        --cell_type_trait_table FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  
  scaled_VJ_per_column(opt)
  RAW_VJ_per_column(opt)
  
  
}


###########################################################################

system.time( main() )
