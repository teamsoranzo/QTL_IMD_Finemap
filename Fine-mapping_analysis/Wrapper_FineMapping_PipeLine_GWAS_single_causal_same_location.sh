#!/bin/bash -l
# FineMapping PipeLine Wrapper Script (Wrapper_FineMapping_PipeLine.sh)
# This is a LSF job wrapper script to run FineMapping_PipeLine.sh


# i) 	Handles job arrays
# ii)	Creating array index file
# iii)	Checking the input data path
# iv)	Setting required parameters for FineMapping

# bash ~/Projects/Scripts/Blueprint/FineMapping/Wrapper_RUN_FineMapping_PipeLine_GWAS_single_causal_same_location.sh

inputfile=$1;
Dis=$2 ## Inflammatory_bowel_disease_CD_Liu_2015_NatGen_GWAS, Rheumatoid_Arthritis_Okada_2014_Nature_GWAS_meta
CellType=$3; ## mono, neut, tcel
QTL=$4; ## gene, K27AC, K4ME1, meth, psi

ArrayIndex=$Dis\_$CellType\_$QTL\-ArrayIndex.txt;


# if [[ $QTL == "psi" ]]; then
	ShortName=$(grep -w $Dis /nfs/team151_data03/PublicData/GWAS_summary_stats/GWAS_unified_generic_formatting/NEW/Short.names | awk '{print $2}');
	cat $inputfile|grep "$Dis\/" | awk '{if($2>=0.98) print $1}' | grep $CellType\_$QTL | \
	tr "/" "\t"|awk '{print $3}'| sed "s/^$CellType\_$QTL\_//" | sed "s/\_$ShortName\_coloc$//" |sort -u | cat -n | sed "s/^[ ]*//" > $ArrayIndex;
# else
# 	cat $inputfile|grep "$Dis\/" | awk '{if($2>=0.99) print $1}' | grep $CellType\_$QTL | \
# 	tr "/" "\t"|awk '{print $3}'| awk -F_ '{print $3}'| sort -u| cat -n | sed "s/^[ ]*//" > $ArrayIndex;
# fi

array_size=$(cat $ArrayIndex | wc -l);
echo $array_size;

##################################################  SET THE PARAMETERS  #########################################################

#
# if [[ ! $ArrayIndex =~ $CellType ]]; then
# 	echo "Wrong files input..";
# 	exit;
# fi
#################################################################################################################################

#
# if [ ! -f $data_PATH ]; then
# 	echo "$data_PATH does not exist !!";
# 	echo "Program terminated..... !!";
# 	exit ;
# fi

if [ ! -s $ArrayIndex ]; then
	echo "Nothing to Fine-map; Program terminated.....";
	exit;
fi


GWAS_PATH=/nfs/team151_data03/PublicData/GWAS_summary_stats/GWAS_unified_generic_formatting/NEW/$Dis\.txt\.gz;

## For Barrett GWAS3 data use below GWAS_PATH
# GWAS_PATH=/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/FineMapping/Data/Barrett_Data/GWAS3/$Dis\.txt\.gz;

## For CAD or anyother GWAS (temporary place) data use below GWAS_PATH
# GWAS_PATH=/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/GWAS_WP10_QTLs_overlap/GWAS_WP10_overlap_with_in_house_summery_stat_WP10_phase2_GWAS_lead_uk10k_1kgp_LD/COLOC/Fine-mapping/$Dis\.txt\.gz;



##### RUN LSF JOB
mkdir -p $Dis/JOB;


#################### Decide memory to make use FARM efficiently ##########################

while IFS= read -r name
do
	line=$(echo $name | awk '{print $2}');
	chr=$(cat /lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Data/Phase2_data/QTL_summary_stat/$CellType\_$QTL\_*10_summary.Beta_changed.SE.Eigen.pval.txt | \
	grep -w $line| awk -F":" '{print $1}');
	if [[ $chr -eq 6 ]]; then
		Mem=10000;
	else
		Mem=2000;

	fi

echo -e "$name\t$Mem";

done < "$ArrayIndex"

# Mem=20000;
##########################################################################################

start=1;
end=500;

loop=$(echo "($array_size/500)+1"|bc);

for i in $(seq 1 $loop); do

	if [ "$end" -gt "$array_size" ]; then

		end=$array_size;

	fi

## Phase 1
# 	echo 'bash /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/FineMapping/FineMapping_PipeLine-wrap_GWAS_single_causal_same_location.sh ${LSB_JOBINDEX} '$ArrayIndex' '$CellType' '$Database_PATH' '$leadSNP_PATH' '$GWAS_PATH' '$GWAS_SAMPLE_SIZE' '$QTL'' | \
# 	bsub -G hematopoiesis -J "JobFetch[$start-$end]" -o $Dis/JOB/%J.%I.$Dis.$QTL.o -q normal -R "select[mem>=$Mem] rusage[mem=$Mem] span[hosts=1]" -M$Mem ;

## Phase 2
	echo 'bash /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/FineMapping/FineMapping_PipeLine-wrap_GWAS_single_causal_same_location_phase2.sh ${LSB_JOBINDEX} '$ArrayIndex' '$CellType' '$GWAS_PATH' '$QTL'' | \
	bsub -G hematopoiesis -J "JobFetch[$start-$end]" -o $Dis/JOB/%J.%I.$Dis.$QTL.o -q normal -R "select[mem>=$Mem] rusage[mem=$Mem] span[hosts=1]" -M$Mem ;

## Phase 2 - Prior-std comparison
# 	echo 'bash /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/FineMapping/FineMapping_PipeLine-wrap_GWAS_single_causal_same_location_phase2_prior-std_comp.sh ${LSB_JOBINDEX} '$ArrayIndex' '$CellType' '$GWAS_PATH' '$QTL'' | \
# 	bsub -G hematopoiesis -J "JobFetch[$start-$end]" -o $Dis/JOB/%J.%I.$Dis.$QTL.o -q long -R "select[mem>=$Mem] rusage[mem=$Mem] span[hosts=1]" -M$Mem ;



	start=$(echo "($i*500)+1"|bc);
	end=$(echo "($i*500)+500"|bc);


done






































