#!/bin/bash -l
# Colocalization PipeLine
# This is a Colocalization pipeline using gwas-pw
# AUTHOR: Kousik Kundu


######################################## PATHS ##########################################
#																						#
#  			Paths and Env. variable for all required softwares and files				#
#########################################################################################


RSCRIPT="/software/R-3.2.2/bin/Rscript";
SQLITE3="/usr/bin/sqlite3";
GWASPW="/nfs/team151/software/gwas-pw-21/gwas-pw-0.21/src/gwas-pw";

GWAS_data="/nfs/team151_data03/PublicData/GWAS_summary_stats/GWAS_unified_generic_formatting/NEW";

# Lead_SNP="/lustre/scratch114/projects/hematopoiesis/Blueprint/limix/data_sep24/eqtl_results/summary_touse/2016_QTL_summary";
Lead_SNP="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Data/Phase2_data/QTL_summary_stat"

########################## INPUT DATA & PARAMETERS SET ##################################
#																						#
#  					All the parameters need to be set here			  					#
#########################################################################################

ID=$1;
infile=$2;
CellType=$3;
QTL=$4;
outdir=$5
GWAS_file=$6;

Disease=$(basename ${GWAS_file%.txt.gz});

Dis=$(grep -w ^$Disease $GWAS_data/Short.names| awk '{print $2}');

TAG=$(cat $infile|awk '{if($1 == "'$ID'") print $2}');	## TAG ID, e.g., ENSG00000100321.10 (for gene_exp), 10:12307226:12307485 (for H3K27ac) etc.

# leadSNP_PATH=$(ls -l $Lead_SNP/$CellType\_$QTL\_*summary.txt.20052016.simple.txt| awk '{print $9}');
leadSNP_PATH=$(ls -l $Lead_SNP/$CellType\_$QTL*_10_summary.Beta_changed.SE.Eigen.pval.txt| awk '{print $9}');

output=$outdir/$CellType\_$QTL\_${TAG%.*}\_$Dis\_coloc;
mkdir $output;

########################################################################################

# cat $GWAS_file| sed '1d' | sort -gk7 | \
# awk '{if($7 != "NA" || $11 != "NA" || $12 != "NA") print $2":"$3"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$11"\t"$12}'| awk '!seen[$1]++' > $output/GWAS.txt;



awk 'NR==FNR {h[$1] = $2"\t"$3"\t"$4"\t"$6"\t"$7; next} {OFS="\t"; print $0,h[$1]}' $outdir/$TAG $Disease/GWAS.txt | \
awk -F"\t" '{if($9!="" && $10 != "" && $11 != ""  && $12 != ""  && $13 != "") print $0}' | \
awk -F"\t" '{OFS="\t"; if ($4 != "A" && $4 != "T" && $4 != "G" && $4 != "C")  print $1,$2,$3,$6,$7,$8,$11,$12,$13; \
else if($4 != $10) print $1,$2,$3,$6,$7,$8*-1,$11,$12,$13; else print $1,$2,$3,$6,$7,$8,$11,$12,$13}' | \
awk '{if($4 == 0 || $5 == 0 || $6 == 0 || $7 == 0 || $8 == 0 || $9 == 0) next; else print $0 }'> $outdir/${TAG%.*}\_master.txt;

## For above command: For the INDELS, the sign will not be changed, since I only checked whether the GWAS SNPs ($4) in either of the nucleotide.

(echo -e "SNPID\tCHR\tPOS\tZ_$Disease\tV_$Disease\tZ_${TAG%.*}\tV_${TAG%.*}" && cat $outdir/${TAG%.*}\_master.txt | \
awk '{OFS="\t"; print $1,$2,$3,$6,$5,$9,$8}' | grep -wv "NA" | grep -wv "Inf" | sort -nk3) > $outdir/${TAG%.*}\_coloc.input;


########################### Colocalization with gwas-pw program ##########################
# 																						 #
#   									  												 #
##########################################################################################
Num=$(cat $outdir/${TAG%.*}\_coloc.input | sed '1d'|wc -l);

$GWASPW -i $outdir/${TAG%.*}\_coloc.input \
-phenos $Disease ${TAG%.*} \
-o $output/$CellType\_$QTL\_${TAG%.*}\_$Dis \
-k $Num > $output/$CellType\_$QTL\_${TAG%.*}\_$Dis.log


##########################################################################################

cat $outdir/${TAG%.*}\_master.txt | awk '{OFS="\t"; print $2,$3,$4,$7}' > $output/$CellType\_$QTL\_${TAG%.*}\_$Dis.pval ;

zcat $output/*.bfs.gz | sed '1d'| awk '{print $2"\t"$3"\t"$15}' > $output/$CellType\_$QTL\_${TAG%.*}\_$Dis.pp;


#################################### Plotting with R #####################################
# 																						 #
#   									  												 #
##########################################################################################


WP10_lead_pos=$(cat $leadSNP_PATH | grep -w ${TAG%.*} | awk -F_ '{print $1}'| awk -F":" '{print $2}');
GWAS_lead_pos=$(cat $output/$CellType\_$QTL\_${TAG%.*}\_$Dis.pval | sort -gk3 | awk '{print $2}'| head -1);

echo $WP10_lead_pos;
echo $GWAS_lead_pos;


$RSCRIPT /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/R_Scripts/Locus_Zoom_plot_for_Coloc.R \
$output/$CellType\_$QTL\_${TAG%.*}\_$Dis.pval \
$output/$CellType\_$QTL\_${TAG%.*}\_$Dis.pp \
$WP10_lead_pos \
$GWAS_lead_pos \
$CellType\_$QTL\_${TAG%.*}\_$Dis \
$output/$CellType\_$QTL\_${TAG%.*}\_$Dis.pdf;


################################# Number of SNPs count ###################################
# 																						 #
#   									  												 #
##########################################################################################

echo -e "GWAS_SNPs\tWP10_SNPs\tColoc_input" >  $output/NumberOfSNPs.txt
snp_start=$(cat $outdir/$TAG| cut -f1| awk -F":" '{print $2}'| sort -n | head -1);
snp_end=$(cat $outdir/$TAG| cut -f1| awk -F":" '{print $2}'| sort -n | tail -1);
snp_chr=$(cat $outdir/$TAG| head -1| awk -F":" '{print $1}');
cat $Disease/GWAS.txt| awk '{if($2 == '$snp_chr' && $3 >= '$snp_start' && $3<= '$snp_end') print $0}' |wc -l | tr "\n" "\t"  >> $output/NumberOfSNPs.txt
# cat $outdir/${TAG%.*}\_master.txt | wc -l|tr "\n" "\t" >> $output/NumberOfSNPs.txt ;
cat $outdir/$TAG  |wc -l | tr "\n" "\t" >> $output/NumberOfSNPs.txt;
cat $outdir/${TAG%.*}\_coloc.input | wc -l >> $output/NumberOfSNPs.txt;




rm $outdir/$TAG;
rm $outdir/${TAG%.*}\_master.txt;
rm $outdir/${TAG%.*}\_coloc.input;
rm $output/*.log;

















