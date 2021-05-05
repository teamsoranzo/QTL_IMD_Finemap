#!/bin/bash -l
# Conditional analysis PipeLine
# Conditional analysis using GCTA


######################################## PATHS ##########################################
#																						#
#  			Paths and Env. variable for all required softwares and files				#
#########################################################################################

# BCFTOOL="/nfs/team151/software/bcftools_2017/bcftools-1.4/bcftools";
BCFTOOL="/nfs/team151/software/Anaconda2/bin/bcftools"; 
PLINK19="/software/hgi/pkglocal/plink-1.90b3w/bin/plink";
# PLINK2="/software/hgi/installs/plink/2.00-alpha/plink2"; ## Farm5
PLINK2="/nfs/team151/software/plink2_18_April_2015/plink";
RSCRIPT="/software/R-3.2.2/bin/Rscript";
SQLITE3="/usr/bin/sqlite3";
GCTA="/nfs/team151/software/GCTA/gcta_1.26.0/gcta64";
QCTOOL="/nfs/team151/software/qctool/qctool_v1.4-linux-x86_64/qctool";

Ref_data="/lustre/scratch116/vr/ref/human/GRCh37/imputation/uk10k+1000g-phase3/vcf";
DataPATH="/lustre/scratch119/realdata/mdt2/teams/soranzo/users/kk8/uk10k+1000G-phase3_PLINK";

GWAS_data="/nfs/team151_data03/PublicData/GWAS_summary_stats/GWAS_unified_generic_formatting/NEW";

########################## INPUT DATA & PARAMETERS SET ##################################
#																						#
#  					All the parameters need to be set here			  					#
#########################################################################################

ID=$1;
infile=$2;

Dis=$(cat $infile|awk '{if($1 == "'$ID'") print $2}');
SNPid=$(cat $infile|awk '{if($1 == "'$ID'") print $3}');

SUFFIX="23_08_17"; ## uk10k+1000g-phase3 PLINK file

PlinkMem=5000
NoOfThreads=8

## Parameter for conditional analysis (GCTA software)
MAF=0.05 ;

# P_val_Threshold_COJO="1e-5"; ## P value threshold set for the step-wise conditional analysis
P_val_Threshold_COJO="5e-8"; ## P value threshold set for the step-wise conditional analysis

##########################################################################################

TEMPDIR="$Dis/$SNPid.cond";
mkdir -p $TEMPDIR



Chr=$(echo $SNPid | awk -F: '{print $1}');
Pos=$(echo $SNPid | awk -F: '{print $2}');

loci_start=$(echo $Pos| awk '{print $0-1000000}');
loci_end=$(echo $Pos| awk '{print $0+1000000}');

if [[ $loci_start -lt 0 ]]; then
	loci_start=1
fi

echo -e "$Dis\t$SNPid";
echo -e "$Chr\t$Pos\t$loci_start\t$loci_end";

################################### Data preparation #####################################
# 																						 #
#   									  												 #
##########################################################################################

GWAS_SAMPLE_SIZE=$(cat /nfs/users/nfs_k/kk8/Projects/GWAS_data/GWAS_sample_size.txt | grep -w $Dis | cut -f2);
echo -e "$Dis\t$GWAS_SAMPLE_SIZE";

echo -e "RSid\tEA\tOA\tEAF\tbeta\tStdErr\tpval\tN" > $TEMPDIR/summary_stat_gcta_input.ma;


if [ "$Dis" == "Advanced_AMD_unified" ] || [ "$Dis" == "CAD_C4D_cardiogram_Nikpay_2015_additive_unified" ] || [ "$Dis" == "GWAS3_CD" ] || [ "$Dis" == "GWAS3_IBD" ] || [ "$Dis" == "GWAS3_UC" ] || [ "$Dis" == "Allergy_Ferreira_2017_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.1.AS.EUR_Malik_2018_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.1.AS.TRANS_Malik_2018_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.2.AIS.EUR_Malik_2018_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.2.AIS.TRANS_Malik_2018_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.4.CES.TRANS_Malik_2018_NatGen_GWAS" ]; then
	zcat $GWAS_data/$Dis\.txt.gz | awk '{OFS="\t"; if($2=='''$Chr''' && $3 >= '''$loci_start''' && $3 <= '''$loci_end''') print $1,$2":"$3,$3,$4,$5,$8,$11,$7,'''$GWAS_SAMPLE_SIZE'''}' | \
	sort -gk8 | awk '!seen[$2]++'| sort -gk3|  awk '{OFS="\t"; if($1=="NA" || $1==".") print $2,$2,$3,$4,$5,$6,$7,$8,$9; else print $0}'| \
	grep -vw "NA" | awk '!seen[$1]++' > $TEMPDIR/summary_stat_gwas.txt;
else
	zcat $GWAS_data/$Dis\.txt.gz | awk '{OFS="\t"; if($2=='''$Chr''' && $3 >= '''$loci_start''' && $3 <= '''$loci_end''') print $1,$2":"$3,$3,$4,$5,$10,$11,$7,'''$GWAS_SAMPLE_SIZE'''}' | \
	sort -gk8 | awk '!seen[$2]++'| sort -gk3|  awk '{OFS="\t"; if($1=="NA" || $1==".") print $2,$2,$3,$4,$5,$6,$7,$8,$9; else print $0}'| \
	grep -vw "NA" | awk '!seen[$1]++' > $TEMPDIR/summary_stat_gwas.txt;
fi

cat $TEMPDIR/summary_stat_gwas.txt | awk '{print $2}'| awk -F":" '{print $1"\t"$2"\t"$2}' > $TEMPDIR/summary_stat_gwas_bedfile.txt;


$BCFTOOL view -R $TEMPDIR/summary_stat_gwas_bedfile.txt $Ref_data/$Chr.sites.bcf | grep -v ^"#" | awk '{print $1,$2,$4,$5,$8}' | awk -F"AF=" '{print $1,$2}' | awk '{print $1,$2,$3,$4,$NF}' | awk -F";" '{print $1}'| awk '{print $1":"$2"\t"$3"\t"$4"\t"$5}' > $TEMPDIR/summary_stat_gwas_AF.txt


# awk 'NR==FNR {h[$2] = $0; next} {print $0,h[$1]}' $TEMPDIR/summary_stat_gwas.txt $TEMPDIR/summary_stat_gwas_AF.txt | awk '{if($1 == $6) print $0}' | \
# awk '{OFS="\t"; if($3 != $8) print $1"_"$2"_"$3,$8,$9,1-$4,$10,$11,$12,$13; else print $1"_"$2"_"$3,$8,$9,$4,$10,$11,$12,$13}' >> $TEMPDIR/summary_stat_gcta_input.ma;
### This is to check whether the effect allele is Alt allele (IF) or the minor allele (ELSE)
if [ "$Dis" == "Advanced_AMD_unified" ] || [ "$Dis" == "GWAS3_CD" ] || [ "$Dis" == "GWAS3_IBD" ] || [ "$Dis" == "GWAS3_UC" ]; then
	awk 'NR==FNR {h[$1] = $0; next} {print $0,h[$2]}' $TEMPDIR/summary_stat_gwas_AF.txt $TEMPDIR/summary_stat_gwas.txt | tr " " "\t" | \
	awk '{if($2 == $10) print $0}'| awk '{OFS="\t"; if($11 == $5 && $12 == $4) print $0; \
	else if ($12 != $4 && $11 != $5) print $1,$2,$3,$12,$11,$6,$7,$8,$9,$10,$11,$12,$13; }'| \
	awk '{OFS="\t"; print $10"_"$11"_"$12,$4,$5,$13,$6,$7,$8,$9}' >> $TEMPDIR/summary_stat_gcta_input.ma;

elif [ "$Dis" == "Allergy_Ferreira_2017_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.1.AS.EUR_Malik_2018_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.1.AS.TRANS_Malik_2018_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.2.AIS.EUR_Malik_2018_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.2.AIS.TRANS_Malik_2018_NatGen_GWAS" ] || [ "$Dis" == "MEGASTROKE.4.CES.TRANS_Malik_2018_NatGen_GWAS" ]  || [ "$Dis" == "Inflammatory_bowel_disease_CD_Liu_2015_NatGen_GWAS_Immunochip" ] || [ "$Dis" == "Inflammatory_bowel_disease_Liu_2015_NatGen_GWAS_Immunochip" ] || [ "$Dis" == "Inflammatory_bowel_disease_UC_Liu_2015_NatGen_GWAS_Immunochip" ] || [ "$Dis" == "Inflammatory_bowel_disease_CD_Liu_2015_NatGen_GWAS" ] || [ "$Dis" == "Inflammatory_bowel_disease_Liu_2015_NatGen_GWAS" ] || [ "$Dis" == "Inflammatory_bowel_disease_UC_Liu_2015_NatGen_GWAS" ] || [ "$Dis" == "Inflammatory_bowel_disease_CD_DeLange_2017_NatGen_GWAS_meta" ] || [ "$Dis" == "Inflammatory_bowel_disease_DeLange_2017_NatGen_GWAS_meta" ] || [ "$Dis" == "Inflammatory_bowel_disease_UC_DeLange_2017_NatGen_GWAS_meta" ] || [ "$Dis" == "Multiple_sclerosis_IMSGC-Patsopoulos_2017_BioRxiv_GWAS_meta" ]; then
	awk 'NR==FNR {h[$1] = $0; next} {print $0,h[$2]}' $TEMPDIR/summary_stat_gwas_AF.txt $TEMPDIR/summary_stat_gwas.txt | tr " " "\t" | \
	awk '{if($2 == $10) print $0}'| awk '{OFS="\t"; if($11 == toupper($5) && $12 == toupper($4)) print $10"_"$11"_"$12,$4,$5,$13,$6,$7,$8,$9; \
	else if ($11 == toupper($4) && $12 == toupper($5)) print $10"_"$11"_"$12,$4,$5,1-$13,$6,$7,$8,$9; }' >> $TEMPDIR/summary_stat_gcta_input.ma;

# else
# 	awk 'NR==FNR {h[$1] = $0; next} {print $0,h[$2]}' $TEMPDIR/summary_stat_gwas_AF.txt $TEMPDIR/summary_stat_gwas.txt | tr " " "\t" | \
# 	awk '{if($2 == $10) print $0}'| awk '{OFS="\t"; if($13>=0.5 && $11 == $4 && $12 == $5) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,1-$13; \
# 	else if ($13<0.5 && $12 == $4 && $11 == $5) print $0; else if($13>=0.5 && $11 != $4 && $12 != $5) print $1,$2,$3,$11,$12,$6,$7,$8,$9,$10,$11,$12,1-$13; \
# 	else if ($13<0.5 && $12 != $4 && $11 != $5) print $1,$2,$3,$12,$11,$6,$7,$8,$9,$10,$11,$12,$13; }'| \
# 	awk '{OFS="\t"; print $10"_"$11"_"$12,$4,$5,$13,$6,$7,$8,$9}' >> $TEMPDIR/summary_stat_gcta_input.ma;
else
	awk 'NR==FNR {h[$1] = $0; next} {print $0,h[$2]}' $TEMPDIR/summary_stat_gwas_AF.txt $TEMPDIR/summary_stat_gwas.txt | tr " " "\t" | \
	awk '{if($2 == $10) print $0}'| awk '{OFS="\t"; if($13>=0.5) print $10"_"$11"_"$12,$4,$5,1-$13,$6,$7,$8,$9; \
	else print $10"_"$11"_"$12,$4,$5,$13,$6,$7,$8,$9;}' >> $TEMPDIR/summary_stat_gcta_input.ma;

fi



cat $TEMPDIR/summary_stat_gcta_input.ma | sed '1d' | awk -F"\t" '{OFS="\t"; if($4 > 0.5) print $1,$2,$3,1-$4,$5,$6,$7,$8; else print $0 }' | awk '{if($4>='''$MAF''') print $1 }' > $TEMPDIR/summary_stat_gcta_input.ma.MAF.txt;
awk 'NR==FNR {h[$1] = $0; next} {print $0,h[$1]}' $TEMPDIR/summary_stat_gcta_input.ma.MAF.txt  $TEMPDIR/summary_stat_gcta_input.ma | awk '{OFS="\t"; if($1 == $9) print $1,$2,$3,$4,$5,$6,$7,$8}' > $TEMPDIR/summary_stat_gcta_input.ma.final.txt;
mv $TEMPDIR/summary_stat_gcta_input.ma.final.txt $TEMPDIR/summary_stat_gcta_input.ma;

# SNP=$(cat $TEMPDIR/summary_stat_gcta_input.ma | grep ^$SNPid\_ | cut -f1);

################################################################################################################
###### In some GWAS, there are some issues with the strand; hence some error might (will) occur during the run -
# Error: none of the given SNPs can be matched to the genotype and summary data.
# For fixing this issue, a non-bobust step can be applied; NO NEED to use this, if you don't get such error

cat $TEMPDIR/summary_stat_gcta_input.ma| tr "_" "\t"| awk '{OFS ="\t"; if($2 != $4 && $2 != $5 && $4 == "G" && $5 == "T") print $1"_"$2"_"$3,"C","A",$6,$7,$8,$9,$10; else if ($2 != $4 && $2 != $5 && $4 == "G" && $5 == "A") print $1"_"$2"_"$3,"C","T",$6,$7,$8,$9,$10; else if ($2 != $4 && $2 != $5 && $4 == "T" && $5 == "G") print $1"_"$2"_"$3,"A","C",$6,$7,$8,$9,$10; else if ($2 != $4 && $2 != $5 && $4 == "T" && $5 == "C") print $1"_"$2"_"$3,"A","G",$6,$7,$8,$9,$10; else if ($2 != $4 && $2 != $5 && $4 == "A" && $5 == "G") print $1"_"$2"_"$3,"T","C",$6,$7,$8,$9,$10; else if ($2 != $4 && $2 != $5 && $4 == "A" && $5 == "C") print $1"_"$2"_"$3,"T","G",$6,$7,$8,$9,$10; else if ($2 != $4 && $2 != $5 && $4 == "C" && $5 == "T") print $1"_"$2"_"$3,"G","A",$6,$7,$8,$9,$10; else if ($2 != $4 && $2 != $5 && $4 == "C" && $5 == "A") print $1"_"$2"_"$3,"G","T",$6,$7,$8,$9,$10; else print $1"_"$2"_"$3,$4,$5,$6,$7,$8,$9,$10;}' > $TEMPDIR/summary_stat_gcta_input.ma.stand_fixed.txt;
mv $TEMPDIR/summary_stat_gcta_input.ma.stand_fixed.txt  $TEMPDIR/summary_stat_gcta_input.ma;

################################################################################################################

SNP=$(cat $TEMPDIR/summary_stat_gcta_input.ma | sed '1d' | awk '{if($7 != 0) print $0}'| sort -gk7| head -1 | cut -f1);


SNPid=$SNP;

######################### Extract variants from big PLINK files ##########################
# 																						 #
#   									  												 #
##########################################################################################

echo -e "$Chr\t$loci_start\t$loci_end\tSet1" > $TEMPDIR/extract_range.txt


$PLINK2 --noweb \
--bfile $DataPATH/$Chr\.$SUFFIX \
--extract range $TEMPDIR/extract_range.txt \
--make-bed \
--keep-allele-order \
--out $TEMPDIR/$Chr\.$SUFFIX\.extracted \
--memory $PlinkMem \
--threads $NoOfThreads

rm $TEMPDIR/extract_range.txt

# cp $DataPATH/$Chr\.$SUFFIX.fam $TEMPDIR/$Chr\.$SUFFIX\.extracted.fam

################################### Conditional analysis #################################
# 																						 #
#   									  												 #
##########################################################################################


echo "$SNPid" > $TEMPDIR/cond.snplist;

echo -e "SNPid\tEAF\tbeta\tStdErr\tpval\tc_beta\tc_pval" > $TEMPDIR/cond.snplist.details;
cat $TEMPDIR/summary_stat_gcta_input.ma | grep -w "$SNPid" | awk '{OFS="\t"; print $1,$4,$5,$6,$7,"-","-"}' >> $TEMPDIR/cond.snplist.details;

Check="0";

while [ $Check -lt "1" ] ; do

# 	$GCTA --bfile $DataPATH/$Chr\.$SUFFIX \
	$GCTA --bfile $TEMPDIR/$Chr\.$SUFFIX\.extracted \
	--chr $Chr \
	--maf $MAF \
	--cojo-file $TEMPDIR/summary_stat_gcta_input.ma \
	--cojo-cond $TEMPDIR/cond.snplist \
	--out  $TEMPDIR/SNPcond

# 	--cojo-collinear 0.98 \

	count=$(cat $TEMPDIR/SNPcond.cma.cojo | sort -gk13 | awk '{if($13 != "NA" && $13 < '''$P_val_Threshold_COJO''') print $0}'| wc -l);
# 	echo "INFO: $count $P_val_Threshold_COJO $Check";
	if [ $count -eq "0" ]; then

		Check=1;

	else

		cat $TEMPDIR/SNPcond.cma.cojo | sort -gk13 | awk '{if($13 != "NA" && $13 != 0 && $13 < '''$P_val_Threshold_COJO''' && $5 > '''$MAF''') print $2}'| head -1 >>  $TEMPDIR/cond.snplist ;
		cat $TEMPDIR/SNPcond.cma.cojo | sort -gk13 | awk '{ OFS="\t"; if($13 != "NA" && $13 != 0 && $13 < '''$P_val_Threshold_COJO''' && $5 > '''$MAF''') print $2,$5,$6,$7,$8,$11,$13}'| head -1 >>  $TEMPDIR/cond.snplist.details ;

		## Sometimes GCTA throws some error, explaining the collinearity problem, In order to tackle that, the
		## option "--cojo-collinear" in GCTA need to be increased. However, if you don't want to increase it,
		## this code is required.
		dupcondsnps=$(cat $TEMPDIR/cond.snplist|sort|uniq -c|sort -nrk1| head -1|awk '{print $1}');
		if [[ $dupcondsnps -gt 1 ]]; then
			sort -u $TEMPDIR/cond.snplist > $TEMPDIR/temp.txt;
			mv $TEMPDIR/temp.txt $TEMPDIR/cond.snplist;
			Check=1;
		fi
	fi

	Check=1;

done


rm $TEMPDIR/summary_stat_gwas* $TEMPDIR/SNPcond.cma.cojo $TEMPDIR/SNPcond.given.cojo $TEMPDIR/*extracted* $TEMPDIR/extract_range.txt
rm $TEMPDIR/summary_stat_gcta_input.*





