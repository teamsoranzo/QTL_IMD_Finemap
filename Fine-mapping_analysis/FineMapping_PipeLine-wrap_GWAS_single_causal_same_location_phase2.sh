#!/bin/bash -l
# FineMapping PipeLine
# This is a fine mapping pipeline using CAVIARBF and FINEMAP
# AUTHOR: Kousik Kundu


######################################## PATHS ##########################################
#																						#
#  			Paths and Env. variable for all required softwares and files				#
#########################################################################################

BCFTOOL="/nfs/team151/software/bcftools_2015/bcftools-1.2/bcftools";
PLINK19="/software/hgi/pkglocal/plink-1.90b3w/bin/plink";
FINEMAP="/nfs/team151/software/Finemap/v1.1/finemap";
FINEMAP_Beta="/nfs/team151/software/Finemap/v1.1.Beta/finemap";
CAVIARBF="/nfs/team151/software/CAVIARBF/v0.1.4.1/caviarbf";
CAVIARBF_MS="/nfs/team151/software/CAVIARBF/v0.1.4.1/model_search";
RSCRIPT="/software/R-3.2.2/bin/Rscript";
SQLITE3="/usr/bin/sqlite3";
GCTA="/nfs/team151/software/GCTA/gcta_1.25.2/gcta64";
QCTOOL="/nfs/team151/software/qctool/qctool_v1.4-linux-x86_64/qctool";

DataPATH="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/VCF_QC/Processed_VCF/Final_VCFs";
Ref_data="/lustre/scratch116/vr/ref/human/GRCh37/imputation/uk10k+1000g-phase3/vcf";


# FINEMAP: This library is needed ! Don't insert this path to your .bashrc,
# as it might create problems.
export LD_LIBRARY_PATH=/software/hgi/pkglocal/gcc-4.9.1/lib64:/software/hgi/pkglocal/gcc-4.9.1/lib:$LD_LIBRARY_PATH

########################## INPUT DATA & PARAMETERS SET ##################################
#																						#
#  					All the parameters need to be set here			  					#
#########################################################################################

ID=$1;
infile=$2;
CellType=$3;

GWAS_file=$4;
QTL=$5;


MAF=0.05;



echo "$ID $infile";

# SNP=$(cat $infile|awk '{if($1 == "'$ID'") print $2}'); ## SNPID (either rs9611170 OR 22:39784845:G:C)
TAG=$(cat $infile|awk '{if($1 == "'$ID'") print $2}');	## TAG ID, e.g., ENSG00000100321.10 (for gene_exp), 10:12307226:12307485 (for H3K27ac) etc.

Database_PATH=/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Data/Phase2_data/QTL_summary_stat/QTL_summary_DB/$CellType\_$QTL.db;
leadSNP_PATH=$(ls -l /lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Data/Phase2_data/QTL_summary_stat/$CellType\_$QTL\_*10_summary.Beta_changed.SE.Eigen.pval.txt| awk '{print $9}');

Dis=$(basename ${GWAS_file%.txt.*});
SNP_Info=$(cat $leadSNP_PATH | awk '{if($3 == "'$TAG'") print $1"\t"$2}');

SAMPLE_SIZE=$(cat /lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/Phenotype_after_PEER_sample_size.txt | grep -w $CellType | grep -w $QTL | cut -f3);
## NOTE: GWAS sample size might NOT be very accurate.
GWAS_SAMPLE_SIZE=$(cat /nfs/users/nfs_k/kk8/Projects/GWAS_data/GWAS_sample_size.txt | grep -w $Dis | cut -f2);


# if [ $CellType == "mono" ]; then
#
# 	SAMPLE_SIZE=194;
#
# elif [ $CellType == "neut" ]; then
#
# 	SAMPLE_SIZE=196;
#
# else
#
# 	SAMPLE_SIZE=169;
#
# fi

echo -e "$Dis\t$CellType\t$QTL\t$SAMPLE_SIZE";
#################################### PRE-PROCESSING ######################################
# 																						 #
#   									  												 #
##########################################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Standard RSid conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

SNPid=$(echo $SNP_Info | awk '{print $1}');
RSid=$(echo $SNP_Info | awk '{print $2}');
Chr=$(echo $SNPid | awk -F: '{print $1}');

if [[ ! $RSid =~ ^"rs" ]]; then
	echo "RSid not given; considered this id ($SNPid) as a standard SNP ID"
fi




TEMPDIR=$Dis/$CellType\_$QTL/$SNPid.$RSid.$TAG.$CellType.temp;
mkdir -p $TEMPDIR

echo "SNP: $SNPid";
echo "TAG: $TAG";

echo "$ID $infile $Database_PATH $leadSNP_PATH $SAMPLE_SIZE $GWAS_SAMPLE_SIZE $GWAS_file";


########################### eQTLs fetching from SQLite database ##########################
# 											and 										 #
#  								  Duplicate entry deleting 			 					 #
##########################################################################################


# sort -t'|' -gk4 | head -1 \
## NOTE - if($4 >= 0.05 && $4<= 0.95) condition in the below code is added later in the analysis (after noticing that the QTL data has MAF <= 1%).
$SQLITE3 $Database_PATH  "select * from qtl WHERE tag = '$TAG'" | \
awk -F"|" '{OFS="\t"; if($8 >= 0.05 && $8<= 0.95) print $2,$3,$1,$4,$5,$6,$7,$8,$9}' | sort -u > $TEMPDIR/Dup_SNP_list.txt; ## Duplicate entry deleted; however, different value for same SNPids are still there.


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

########################### Z score calculation for the LD SNPs ##########################
# 										AND												 #
#   						    Prior probability set									 #
##########################################################################################



echo "RSid alt ref alt.AF beta StdErr pval N StdID GeneID chr Pos"|tr " " "\t" > $TEMPDIR/SNP_list_full_info.txt; ### Heading ###


while IFS2= read -r name; do line=`echo $name | awk '{print $3}'`;

# RSid	alt	ref	alt.AF	beta	StdErr	pval	N	z.score	StdID	GeneID	chr	Pos 1:168944348_G_A
# chr:pos_ref_alt	rsid	phenotypeID	p.value	beta	Bonferroni.p	FDR	alt.AF	std.error
s_chr=$(echo $line | awk -F\: '{print $1}');
s_pos=$(echo $line | awk -F\: '{print $2}'| awk -F\_ '{print $1}');
s_ref=$(echo $line | awk -F\_ '{print $2}');
s_alt=$(echo $line | awk -F\_ '{print $3}');
# z_score=$(echo $name | awk '{if($9 == "Inf") print "0"; else printf ("%.15f\n", $5/$9)}');
# s_stdid=$(echo $line|tr "_" ":");

echo -e "$name\t$s_chr\t$s_pos\t$s_ref\t$s_alt\t$line";
# echo -e "$SNPfetch\t$name";

done < $TEMPDIR/SNP_list.txt | awk '{OFS="\t"; print $1,$13,$12,$8,$5,$9,$4,"'$SAMPLE_SIZE'",$14,$2,$10,$11}' | sort -nk12  >> $TEMPDIR/SNP_list_full_info.txt;

# (head -n 1 $TEMPDIR/SNP_list_full_info.txt && tail -n +2 $TEMPDIR/SNP_list_full_info.txt | sort -nk4) > newfile


$RSCRIPT /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/FineMapping/R_Scripts/Z_score_from_beta_and_pval.R \
$TEMPDIR/SNP_list_full_info.txt \
$TEMPDIR/summary_stat.txt \
#

cat $TEMPDIR/summary_stat.txt | sort -nk13 | awk '{print $1,$9}'|sed 1d > $TEMPDIR/data.z
## Sorting based on chromosome position (-nk13) is very important, as PLINK (e.g. --extract) does it automatically

# Overwrite "$TEMPDIR/LD\_$r2\.txt" as some SNPs are filtered out in the limix run
cat $TEMPDIR/summary_stat.txt |sort -nk13 | awk '{print $10}' |sed 1d > $TEMPDIR/SNPid_list.txt;
# cat $TEMPDIR/summary_stat.txt | sort -nk13 | awk '{print $1}' |sed 1d > $TEMPDIR/SNPid_list_RSid.txt;
## Sorting based on chromosome position (-nk13) is very important, as PLINK (e.g. --extract) does it automatically

NoOfSNPs=$(cat $TEMPDIR/summary_stat.txt|wc -l| awk '{print $1-1}');


############# GWAS SNPs collections ##########
pos_start=$(cat $TEMPDIR/summary_stat.txt | sed '1d'| awk '{print $13}'| sort -n | head -1);
pos_end=$(cat $TEMPDIR/summary_stat.txt | sed '1d'| awk '{print $13}'| sort -n | tail -1);


if [ "$Dis" == "Advanced_AMD_unified" ] || [ "$Dis" == "CAD_C4D_cardiogram_Nikpay_2015_additive_unified" ] || [ "$Dis" == "GWAS3_CD" ] || [ "$Dis" == "GWAS3_IBD" ] || [ "$Dis" == "GWAS3_UC" ] || [ "$Dis" == "Allergy_Ferreira_2017_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.1.AS.EUR_Malik_2018_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.1.AS.TRANS_Malik_2018_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.2.AIS.EUR_Malik_2018_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.2.AIS.TRANS_Malik_2018_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.4.CES.TRANS_Malik_2018_NatGen_GWAS" ]; then

	zcat $GWAS_file | awk '{if($2=='''$Chr''' && $3 >= '''$pos_start''' && $3 <= '''$pos_end''') print $1,$2":"$3,$7,$12,$2,$3,$4,$5,$6,$8,$11,'''$GWAS_SAMPLE_SIZE'''}' | \
	sort -gk3 | awk '!seen[$2]++'| sort -gk2|  awk '{OFS="\t"; if($1=="NA" || $1=="." || $1 !~ /^rs/) print $2,$2,$5,$6,$3,$4,$7,$8,$9,$10,$11,$12; else print $1,$2,$5,$6,$3,$4,$7,$8,$9,$10,$11,$12}'> $TEMPDIR/summary_stat_gwas.txt;

else

	zcat $GWAS_file | awk '{if($2=='''$Chr''' && $3 >= '''$pos_start''' && $3 <= '''$pos_end''') print $1,$2":"$3,$7,$12,$2,$3,$4,$5,$6,$10,$11,'''$GWAS_SAMPLE_SIZE'''}' | \
	sort -gk3 | awk '!seen[$2]++'| sort -gk2|  awk '{OFS="\t"; if($1=="NA" || $1=="." || $1 !~ /^rs/) print $2,$2,$5,$6,$3,$4,$7,$8,$9,$10,$11,$12; else print $1,$2,$5,$6,$3,$4,$7,$8,$9,$10,$11,$12}'> $TEMPDIR/summary_stat_gwas.txt;

fi

####################################### MAF threshold ####################################
# 																						 #
##   If you want to put MAF threshold, un-command the below mentioned code 				 #
#   									  												 #
##########################################################################################

mv $TEMPDIR/summary_stat_gwas.txt $TEMPDIR/summary_stat_gwas_all.txt;
MAF_evi=$(cat $TEMPDIR/summary_stat_gwas_all.txt | awk '{if($9 == "NA") print $0}' | wc -l);

if [[ $MAF_evi -gt 0 ]]; then
	cat $TEMPDIR/summary_stat_gwas_all.txt | awk '{print $2}'| awk -F":" '{print $1"\t"$2"\t"$2}' > $TEMPDIR/summary_stat_gwas_bedfile.txt;
	$BCFTOOL view -R $TEMPDIR/summary_stat_gwas_bedfile.txt $Ref_data/$Chr.sites.bcf | grep -v ^"#" | awk '{print $1,$2,$4,$5,$8}' | awk -F"AF=" '{print $1,$2}' | awk '{print $1,$2,$3,$4,$NF}' | awk -F";" '{print $1}'| awk '{print $1":"$2"\t"$3"\t"$4"\t"$5}' > $TEMPDIR/summary_stat_gwas_AF.txt;

	awk 'NR==FNR {h[$1] = $0; next} {print $0,h[$2]}' $TEMPDIR/summary_stat_gwas_AF.txt $TEMPDIR/summary_stat_gwas_all.txt | tr " " "\t" | \
	awk '{if($2 == $13) print $0}'| awk '{OFS="\t"; if($16>=0.5 && $9 =="NA") print $1,$2,$3,$4,$5,$6,$7,$8,1-$16,$10,$11,$12; \
	else if($16<0.5 && $9 =="NA") print $1,$2,$3,$4,$5,$6,$7,$8,$16,$10,$11,$12; else if($9 != "NA") print $0}'> $TEMPDIR/summary_stat_gwas_all_new.txt;
	mv $TEMPDIR/summary_stat_gwas_all_new.txt $TEMPDIR/summary_stat_gwas_all.txt;
	rm $TEMPDIR/summary_stat_gwas_bedfile.txt $TEMPDIR/summary_stat_gwas_AF.txt;
fi

################### MAF cut-off #####################
cat $TEMPDIR/summary_stat_gwas_all.txt | awk -F"\t" '{OFS="\t"; if($9>0.5) print $1,$2,$3,$4,$5,$6,$7,$8,1-$9,$10,$11,$12; else print $0 }' | awk -F"\t" '{if($9 >= '''$MAF''') print $0}' > $TEMPDIR/summary_stat_gwas.txt;
rm $TEMPDIR/summary_stat_gwas_all.txt;

#########################################################################################


### Adding a header in Summary_stat_gwas file
(echo -e 'RSid\tSNPid\tChr\tPos\tpval\tZ_score\tEA\tOA\tMAF\tbeta\tse\tSample_size' && cat $TEMPDIR/summary_stat_gwas.txt) > $TEMPDIR/summary_stat_gwas1.txt;
mv $TEMPDIR/summary_stat_gwas1.txt $TEMPDIR/summary_stat_gwas.txt
# cat $TEMPDIR/summary_stat_gwas1.txt | grep -v rs10562650 | grep -v rs13029144 | grep -v rs6731125 | grep -v rs4667282 | grep -v rs4667283 | grep -v rs7573465 | grep -v rs1449263 | grep -v rs2124440  > $TEMPDIR/summary_stat_gwas.txt

cat $TEMPDIR/summary_stat_gwas.txt| sed '1d' | awk -F"\t" '{if($6 != "NA" || $6 != "." || $6 != "") print $1,$6}' | grep -vw "NA" | awk '!seen[$1]++' > $TEMPDIR/data_gwas.z; ## Sometime SNPid is different but RSid is same - hence get rid of those. E.g. AMD disease GWAS - rs71689664/7:99662511 and rs71689664/7:99662512.

NoOfSNPs_gwas=$(cat $TEMPDIR/data_gwas.z|wc -l);

########################### Correlation matrix for FineMapping ###########################
# 																						 #
#   									  												 #
##########################################################################################

## Since this script is only assuming one causal SNP in a locus, the ld matrix is practically
## not needed, and therefore here I am creating a dummy ld matrix for CAVIARBF and FINEMAP input requirements.

$RSCRIPT /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/FineMapping/R_Scripts/r_calculation_dummay.R \
$NoOfSNPs \
$TEMPDIR/data.ld

$RSCRIPT /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/FineMapping/R_Scripts/r_calculation_dummay.R \
$NoOfSNPs_gwas \
$TEMPDIR/data_gwas.ld


NoOfCausalSNP=1;

echo "NUMBER OF CAUSAL SNPs: $NoOfCausalSNP";
echo "Number of SNPs in WP10: $NoOfSNPs";
echo "WP10 Sample Size: $SAMPLE_SIZE"

echo "Number of SNPs in GWAS: $NoOfSNPs_gwas";
echo "GWAS Sample Size: $GWAS_SAMPLE_SIZE";
############################# FineMapping with FINEMAP program ###########################
# 																						 #
#   									  												 #
##########################################################################################

# data.k (prior probability for FINEMAP) - This parameter is not needed for v1.1 onwards (still can be used)
$RSCRIPT  /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/FineMapping/R_Scripts/prior_prob_finemap.R \
$NoOfCausalSNP \
$NoOfSNPs \
$TEMPDIR/data.k

echo "z;ld;snp;config;k;log;n-ind" > $TEMPDIR/data.master
echo "$TEMPDIR/data.z;$TEMPDIR/data.ld;$TEMPDIR/data.snp;$TEMPDIR/data.config;$TEMPDIR/data.k;$TEMPDIR/data.log;$SAMPLE_SIZE" >> $TEMPDIR/data.master

echo "z;ld;snp;config;k;log;n-ind" > $TEMPDIR/data_gwas.master
echo "$TEMPDIR/data_gwas.z;$TEMPDIR/data_gwas.ld;$TEMPDIR/data_gwas.snp;$TEMPDIR/data_gwas.config;$TEMPDIR/data.k;$TEMPDIR/data_gwas.log;$GWAS_SAMPLE_SIZE" >> $TEMPDIR/data_gwas.master


$FINEMAP --sss --in-files $TEMPDIR/data.master --corr-config 0.8 --n-causal-max $NoOfCausalSNP --prior-std 0.3 --log --n-ind $SAMPLE_SIZE
$FINEMAP --sss --in-files $TEMPDIR/data_gwas.master --corr-config 0.8 --n-causal-max $NoOfCausalSNP --prior-std 0.1281429 --log --n-ind $GWAS_SAMPLE_SIZE


## Modified BF after pruning- "--corr-pruning" option is important
# $FINEMAP_Beta --prune --in-files $TEMPDIR/data.master --corr-pruning 0.8 --corr-config 0.8 --n-causal-max $NoOfCausalSNP --prior-std 0.1281429 --log --n-ind $SAMPLE_SIZE


######################################### NOTE ###########################################
# The default --prior-std id 0.05, but when comparing with CAVIARBF, one needs to use the
# same --prior-std for both.
#
# --prior-std 0.1281429 looks more convincing with eQTL data than --prior-std 0.05 (default)
## Higher value indicates signal is more spread and lower value indicates otherwise.



############################# FineMapping with CAVIARBF program ##########################
# 																						 #
#   									  												 #
##########################################################################################


$CAVIARBF -z $TEMPDIR/data.z -r $TEMPDIR/data.ld -t 0 -a 0.3 -n $SAMPLE_SIZE -c $NoOfCausalSNP -o $TEMPDIR/data.bf

# $CAVIARBF_MS -i $TEMPDIR/data.bf -o $TEMPDIR/data.prior0 -m $NoOfSNPs -p 0 -sex > $TEMPDIR/temp.console
$CAVIARBF_MS -i $TEMPDIR/data.bf -o $TEMPDIR/data.prior0 -m $NoOfSNPs -p 0  > $TEMPDIR/temp.console


$CAVIARBF -z $TEMPDIR/data_gwas.z -r $TEMPDIR/data_gwas.ld -t 0 -a 0.1281429 -n $GWAS_SAMPLE_SIZE -c $NoOfCausalSNP -o $TEMPDIR/data_gwas.bf
$CAVIARBF_MS -i $TEMPDIR/data_gwas.bf -o $TEMPDIR/data_gwas.prior0 -m $NoOfSNPs_gwas -p 0  > $TEMPDIR/temp_gwas.console

############ POST processing of CAVIARBF output (only *.prior0.marginal) #################

echo -e "snp\tsnp_prob\tindex" > $TEMPDIR/data.prior0.marginal.mod
while CBF= read -r line; do
lineno=$(echo $line | awk '{print $1}');
RSid_CBF=$(cat $TEMPDIR/data.z| head -$lineno| tail -1 | awk '{print $1}');

echo $line | awk '{print "'$RSid_CBF'",$2,$1}';

done < $TEMPDIR/data.prior0.marginal  >> $TEMPDIR/data.prior0.marginal.mod

sleep 5;
rm $TEMPDIR/temp.console $TEMPDIR/data.bf $TEMPDIR/data.prior0.statistics

##############################

echo -e "snp\tsnp_prob\tindex" > $TEMPDIR/data_gwas.prior0.marginal.mod
while CBF= read -r line; do
lineno=$(echo $line | awk '{print $1}');
RSid_CBF=$(cat $TEMPDIR/data_gwas.z| head -$lineno| tail -1 | awk '{print $1}');

echo $line | awk '{print "'$RSid_CBF'",$2,$1}';

done < $TEMPDIR/data_gwas.prior0.marginal  >> $TEMPDIR/data_gwas.prior0.marginal.mod

sleep 5;
rm $TEMPDIR/temp_gwas.console $TEMPDIR/data_gwas.bf $TEMPDIR/data_gwas.prior0.statistics


##########################################################################################

sleep 5;
rm $TEMPDIR/data.ld $TEMPDIR/data.k $TEMPDIR/data.master
rm $TEMPDIR/data_gwas.ld $TEMPDIR/data_gwas.master
rm $TEMPDIR/Dup* $TEMPDIR/SNP*;


# rm $TEMPDIR/summary_stat_gwas_all.txt;

# Other summary*
# rm $TEMPDIR/summary* (Don't delete it; this is required for after processing)







