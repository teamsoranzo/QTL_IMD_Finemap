#!/bin/bash -l
# Conditional analysis PipeLine
# Conditional analysis using GCTA



######################################## PATHS ##########################################
#																						#
#  			Paths and Env. variable for all required softwares and files				#
#########################################################################################

# BCFTOOL="/nfs/team151/software/bcftools_2017/bcftools-1.4/bcftools";
BCFTOOL="/nfs/team151/software/Anaconda2/bin/bcftools"; 
# PLINK19="/software/hgi/pkglocal/plink-1.90b3w/bin/plink";
PLINK2="/nfs/team151/software/plink2_18_April_2015/plink";
# RSCRIPT="/software/R-3.2.2/bin/Rscript";
# SQLITE3="/usr/bin/sqlite3";
RSCRIPT="/software/hgi/installs/anaconda3/envs/hgi_base/bin/Rscript"  ## Farm5 (hgi_base env)
SQLITE3="/software/hgi/installs/anaconda3/envs/hgi_base/bin/sqlite3"  ## Farm5 (hgi_base env)
GCTA="/nfs/team151/software/GCTA/gcta_1.26.0/gcta64";
QCTOOL="/nfs/team151/software/qctool/qctool_v1.4-linux-x86_64/qctool";

DataPATH="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/VCF_QC/Processed_VCF/Final_VCFs";
GeneAnno="/lustre/scratch114/projects/hematopoiesis/Blueprint/limix/data_sep24/final_pheno/RNA-seq/gene_20151028_anno_sort.txt";


########################## INPUT DATA & PARAMETERS SET ##################################
#																						#
#  					All the parameters need to be set here			  					#
#########################################################################################

ID=$1;
infile=$2;

CellType=$(cat $infile|awk '{if($1 == "'$ID'") print $2}');	## TAG ID, e.g., ENSG00000100321.10 (for gene_exp), 10:12307226:12307485 (for H3K27ac) etc.
QTL=$(cat $infile|awk '{if($1 == "'$ID'") print $3}');
TAG=$(cat $infile|awk '{if($1 == "'$ID'") print $4}');

Database_PATH=/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Data/Phase2_data/QTL_summary_stat/QTL_summary_DB/$CellType\_$QTL.db;
leadSNP_PATH=$(ls -l /lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Data/Phase2_data/QTL_summary_stat/$CellType\_$QTL\_*10_summary.Beta_changed.SE.Eigen.pval.txt| awk '{print $9}');


SUFFIX="BPWP10_23_05_17"; ## Phase 2



## Parameter for conditional analysis (GCTA software)
MAF=0.05 ;

P_val_Threshold_COJO="1e-5"; ## P value threshold set for the step-wise conditional analysis

##########################################################################################


SNP_Info=$(cat $leadSNP_PATH | awk '{if($3 == "'$TAG'") print $1"\t"$2}');


TEMPDIR="$CellType"_"$QTL/$CellType.$QTL.$TAG.cond";
mkdir -p $TEMPDIR

if [ $CellType == "mono" ]; then

	SAMPLE_SIZE=194;

elif [ $CellType == "neut" ]; then

	SAMPLE_SIZE=196;

else

	SAMPLE_SIZE=169;

fi


SNPid=$(echo $SNP_Info | awk '{print $1}');
RSid=$(echo $SNP_Info | awk '{print $2}');
Chr=$(echo $SNPid | awk -F: '{print $1}');

if [[ ! $RSid =~ ^"rs" ]]; then
	echo "RSid not given; considered this id ($SNPid) as a standard SNP ID"
fi


################################### Conditional analysis #################################
# 																						 #
#   									  												 #
##########################################################################################

$SQLITE3 $Database_PATH  "select * from qtl WHERE tag = '$TAG'" | \
awk -F"|" '{OFS="\t"; print $2,$3,$1,$4,$5,$6,$7,$8,$9}' |sort -u > $TEMPDIR/Dup_SNP_list.txt; ## Duplicate entry deleted; however, different value for same SNPids are still there (this is observed in phase1; hasn't been checked for phase2 - used just for the safer side).


No_Dup=$(cat $TEMPDIR/Dup_SNP_list.txt |awk '{print $1}'|sort|uniq -c | awk '{if($1 > 1) print $2}'| wc -l);

if [[ $No_Dup -eq  0 ]]; then
	echo "Checked in";
	cp $TEMPDIR/Dup_SNP_list.txt $TEMPDIR/SNP_list.txt

else
	echo "Checked out";
	cat $TEMPDIR/Dup_SNP_list.txt |awk '{print $1}'|sort|uniq -c | awk '{if($1 > 1) print $2}' > $TEMPDIR/Dup.txt;

	while IFS= read -r name; do

		grep $name $TEMPDIR/Dup_SNP_list.txt | sort -gk3| sed 1d >> $TEMPDIR/Dup_info.txt;

	done < $TEMPDIR/Dup.txt;

	grep -vf $TEMPDIR/Dup_info.txt $TEMPDIR/Dup_SNP_list.txt > $TEMPDIR/SNP_list.txt;
	rm $TEMPDIR/Dup_info.txt ;
fi


######################### Extract variants from big VCF files ############################
# 																						 #
#   									  												 #
##########################################################################################

# Pos_start=$(cat $TEMPDIR/SNP_list.txt | cut -f3| tr "_" "\t" | tr ":" "\t" | cut -f2| sort -n | head -1);
# Pos_end=$(cat $TEMPDIR/SNP_list.txt | cut -f3| tr "_" "\t" | tr ":" "\t" | cut -f2| sort -n | tail -1);
#
# $BCFTOOL view -r $Chr\:$Pos_start\-$Pos_end $DataPATH/$Chr\.$SUFFIX\.vcf.gz

# cat  $TEMPDIR/SNP_list.txt | cut -f3 > $TEMPDIR/extract.txt
#
# $PLINK19 --noweb \
# --bfile  $DataPATH/Final_PLINKs/$Chr\.$SUFFIX \
# --extract $TEMPDIR/extract.txt \
# --make-bed \
# --keep-allele-order \
# --out $TEMPDIR/$Chr\.$SUFFIX\.extracted \
# --memory $PlinkMem \
# --threads $NoOfThreads
#
# rm $TEMPDIR/extract.txt

################################### Conditional analysis #################################
# 																						 #
#   									  												 #
##########################################################################################

echo -e "RSid\talt\tref\talt.AF\tbeta\tStdErr\tpval\tN" > $TEMPDIR/summary_stat_gcta_input.ma;

cat $TEMPDIR/SNP_list.txt | cut -f3| tr "_" "\t" | awk '{print $1"_"$2"_"$3"\t"$2"\t"$3}' > $TEMPDIR/SNP_list_SNPid.txt;

awk 'NR==FNR {h[$1] = $0; next} {print $0,h[$3]}' $TEMPDIR/SNP_list_SNPid.txt $TEMPDIR/SNP_list.txt | \
awk '{OFS="\t"; if($3 == $10) print $3,$12,$11,$8,$5,$9,$4,"'$SAMPLE_SIZE'"}' >> $TEMPDIR/summary_stat_gcta_input.ma;



echo "$SNPid" > $TEMPDIR/cond.snplist;

echo -e "SNPid\tAlt.AF\tbeta\tStdErr\tpval\tc_beta\tc_pval" > $TEMPDIR/cond.snplist.details;
cat $TEMPDIR/summary_stat_gcta_input.ma | grep -w "$SNPid" | awk '{OFS="\t"; print $1,$4,$5,$6,$7,"-","-"}' >> $TEMPDIR/cond.snplist.details;

Check="0";

while [ $Check -lt "1" ] ; do

	$GCTA --bfile $DataPATH/Final_PLINKs/$Chr\.$SUFFIX \
	--chr $Chr \
	--maf $MAF \
	--cojo-file $TEMPDIR/summary_stat_gcta_input.ma \
	--cojo-cond $TEMPDIR/cond.snplist \
	--out  $TEMPDIR/SNPcond

# 	--cojo-collinear 0.98 \
# 	--cojo-actual-geno \

	count=$(cat $TEMPDIR/SNPcond.cma.cojo | sort -gk13 | awk '{if($13 != "NA" && $13 < '''$P_val_Threshold_COJO''') print $0}'| wc -l);
# 	echo "INFO: $count $P_val_Threshold_COJO $Check";
	if [ $count -eq "0" ]; then

		Check=1;

	else

		cat $TEMPDIR/SNPcond.cma.cojo | sort -gk13 | awk '{if($13 != "NA" && $13 < '''$P_val_Threshold_COJO''') print $2}'| head -1 >>  $TEMPDIR/cond.snplist ;
		cat $TEMPDIR/SNPcond.cma.cojo | sort -gk13 | awk '{ OFS="\t"; if($13 != "NA" && $13 < '''$P_val_Threshold_COJO''') print $2,$5,$6,$7,$8,$11,$13}'| head -1 >>  $TEMPDIR/cond.snplist.details ;

		## Sometimes GCTA throws some error, explaining the collinearity problem, In order to tackle that, the
		## option "--cojo-collinear" in GCTA needs to be increased. However, if you don't want to increase it,
		## this code is required.
		dupcondsnps=$(cat $TEMPDIR/cond.snplist|sort|uniq -c|sort -nrk1| head -1|awk '{print $1}');
		if [[ $dupcondsnps -gt 1 ]]; then
			sort -u $TEMPDIR/cond.snplist > $TEMPDIR/temp.txt;
			mv $TEMPDIR/temp.txt $TEMPDIR/cond.snplist;
			Check=1;
		fi
	fi

done

# rm $TEMPDIR/Dup_SNP_list.txt $TEMPDIR/SNPcond.cma.cojo $TEMPDIR/SNPcond.given.cojo $TEMPDIR/SNP_list_SNPid.txt $TEMPDIR/SNP_list.txt
# rm $TEMPDIR/summary_stat_gcta_input.ma
