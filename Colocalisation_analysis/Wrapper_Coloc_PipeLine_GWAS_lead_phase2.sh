#!/bin/bash -l
# Colocalization PipeLine Wrapper Script (Wrapper_Coloc_PipeLine.sh)
# This is a LSF job wrapper script to run Colocalization_GWAS_WP10.sh


# i) 	Handles job arrays
# ii)	Creating array index file
# iii)	Checking the input data path
# iv)	Setting required parameters for Colocalization


######################################## PATHS ##########################################
#																						#
#  			Paths and Env. variable for all required softwares and files				#
#########################################################################################


# RSCRIPT="/software/R-3.2.2/bin/Rscript";
# SQLITE3="/usr/bin/sqlite3";

RSCRIPT="/software/hgi/installs/anaconda3/envs/hgi_base/bin/Rscript"; ## Farm5 (hgi_base env)
SQLITE3="/software/hgi/installs/anaconda3/envs/hgi_base/bin/sqlite3"; ## Farm5 (hgi_base env)

WP10_GWAS_Overlap="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/GWAS_WP10_QTLs_overlap/GWAS_WP10_overlap_with_in_house_summery_stat_WP10_phase2_GWAS_lead_uk10k_1kgp_LD";
# WP10_GWAS_Overlap="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/GWAS_WP10_QTLs_overlap/GWAS_loci_LDetect_Pickrell/GWAS3_benchmarking/GWAS_WP10_overlap";

GWAS_data="/nfs/team151_data03/PublicData/GWAS_summary_stats/GWAS_unified_generic_formatting/NEW";

###### For Barrett data ######
# GWAS_data="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/FineMapping/Data/Barrett_Data/GWAS3";
###### For AMD and CAD ######
## I will delete this zip after analysis
# GWAS_data="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/GWAS_WP10_QTLs_overlap/GWAS_WP10_overlap_with_in_house_summery_stat_WP10_phase2_GWAS_lead_uk10k_1kgp_LD/COLOC"
##############################


Database="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Data/Phase2_data/QTL_summary_stat/QTL_summary_DB";
Lead_SNP="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Data/Phase2_data/QTL_summary_stat";


########################## INPUT DATA & PARAMETERS SET ##################################
#																						#
#  					All the parameters need to be set here			  					#
#########################################################################################

inputfile=$1;

#########################################################################################

cat $inputfile | \
while read line; do
	CellType=$(echo $line| awk '{print $1}');
	QTL=$(echo $line| awk '{print $2}');
	Disease=$(echo $line| awk '{print $3}');

	infile=$WP10_GWAS_Overlap/$Disease\_loci/$CellType\_$QTL\_WP10_GWAS_overlap_SNPs.txt;
# 	infile=$WP10_GWAS_Overlap/$Disease\_LDetect_loci/$CellType\_$QTL\_WP10_GWAS_overlap_SNPs.txt;

	GWAS_file=$GWAS_data/$Disease.txt.gz;

	echo $CellType\_$QTL $GWAS_file $infile;

	dir=$(basename $GWAS_file);
	file=$(basename $infile);
	ArrayIndex=${dir%.txt.*}/${file%.*}-ArrayIndex.txt;

	database_PATH=$Database/$CellType\_$QTL.db;
# 	database_PATH=$Lead_SNP/$CellType\_$QTL\*_10_all_summary.Beta_changed.SE.Eigen.pval.txt;
	leadSNP_PATH=$(ls -l $Lead_SNP/$CellType\_$QTL*_10_summary.Beta_changed.SE.Eigen.pval.txt| awk '{print $9}');


	outdir=${dir%.txt.*}/$CellType\_$QTL;
	mkdir -p $outdir


	cat $infile| awk '{if($3 == 1) print $4}'|sort -u| grep -f - $leadSNP_PATH | awk '{if($7<=0.05) print $3}'|sort -u | \
	while read tag; do
		$SQLITE3 $database_PATH  "select * from qtl WHERE tag = '$tag'" | awk -F"|" '{if($4 < 1) print $1,$4,$5}'> $outdir/$tag.txt;
		$RSCRIPT ~/Projects/Scripts/Blueprint/FineMapping/R_Scripts/Z_score_from_beta_and_pval_for_COLOC.R $outdir/$tag.txt $outdir/$tag.out ;
		cat $outdir/$tag.out | sed 's/_/\t/g' > $outdir/$tag.txt;
		rm $outdir/$tag.out;

	done

	ls -l $outdir | sed '1d'| cat -n| awk '{print $1"\t"$NF}' > $ArrayIndex ;

	array_size=$(cat $ArrayIndex | wc -l);

	##### RUN LSF JOB
	JOB=$outdir/JOB;
	mkdir -p $JOB;

	echo -e "Tested QTL loci - $outdir\t$array_size";

	file_G="$Disease/GWAS.txt";

	if [ ! -f "$file_G" ]; then
		echo "GWAS file creating.....";
		zcat $GWAS_file| sed '1d' | sort -gk7 | \
		awk '{if($7 != "NA" && $11 != "NA" && $12 != "NA") print $2":"$3"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$11"\t"$12}'| awk '!seen[$1]++' > $Disease/GWAS.txt;

	fi
	if [ ! -s $Disease/GWAS.txt ]; then
		echo "Error: GWAS file has missing data .. (either Pval, Z_score, or SE)";
		rm -rf ${dir%.txt.*};
		continue;
	fi;



start=1;
end=1000;

loop=$(echo "($array_size/1000)+1"|bc);

	for i in $(seq 1 $loop); do

		if [ "$end" -gt "$array_size" ]; then

			end=$array_size;
		fi


# 		echo 'bash /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/FineMapping/Colocalization_GWAS_WP10.sh ${LSB_JOBINDEX} '$ArrayIndex' '$CellType' '$QTL' '$outdir' '$GWAS_file'' | \
# 		bsub -G hematopoiesis -J "JobFetch[1-$array_size]" -o $JOB/%J.%I.${file%.*}.o -q small -R "select[mem>=1000] rusage[mem=1000] span[hosts=1]" -M1000 ;

		echo 'bash /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/FineMapping/Colocalization_GWAS_WP10.sh ${LSB_JOBINDEX} '$ArrayIndex' '$CellType' '$QTL' '$outdir' '$GWAS_file'' | \
		bsub -G hematopoiesis -J "JobFetch[$start-$end]" -o $JOB/%J.%I.${file%.*}.o -q small -R "select[mem>=5000] rusage[mem=5000] span[hosts=1]" -M5000 ;


		start=$(echo "($i*1000)+1"|bc);
		end=$(echo "($i*1000)+1000"|bc);


	done


done














































