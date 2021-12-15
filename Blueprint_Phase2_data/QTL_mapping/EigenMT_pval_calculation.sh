
#!/bin/bash -l


###################################  INPUT DATA	 ########################################
#																						#
#				   All the input arguments are being received 	 	  					#
#									  													#
#########################################################################################

### Other parameters ###


Chr=$1;
CellType=$2;
Tag=$3;

echo "Chr = $Chr";
echo "CellType = $CellType";
echo "Feature = $Tag";

TEMPDIR=$CellType\_$Tag\_TEMP
mkdir -p $TEMPDIR;

######################################## PATHS ##########################################
#																						#
#  						Paths for all required softwares and files		  				#
#########################################################################################

PYTHON="/usr/bin/python";
Rscript="/software/R-3.2.2/bin/Rscript";
PY_EigenMT="/nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/Re_analysis/Limix/eigenMT.py";

GenoMatrixScript="/nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/Re_analysis/Limix/hd5_to_txt_genotype_file.R";


QTL_path="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/qtl_result/$CellType\_$Tag*/runs/split_folder/$CellType\_$Tag_*_all_summary.Beta_changed.SE.txt";

mono_gene_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/sorted_gene_20151028_anno_sort.txt";
neut_gene_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/sorted_gene_20151028_anno_sort.txt";
tcel_gene_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/sorted_gene_20151028_anno_sort.txt";

mono_K27AC_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/anno.matrix.01092016.log2rpm.MON.K27AC.BEDMAP.ALL.K27AC.bed";
neut_K27AC_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/anno.matrix.01092016.log2rpm.NEU.K27AC.BEDMAP.ALL.K27AC.bed";
tcel_K27AC_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/anno.matrix.01092016.log2rpm.TCELL.K27AC.BEDMAP.ALL.K27AC.bed";

mono_K4ME1_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/anno.matrix.01092016.log2rpm.MON.K4ME1.BEDMAP.ALL.K4ME1.bed";
neut_K4ME1_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/anno.matrix.01092016.log2rpm.NEU.K4ME1.BEDMAP.ALL.K4ME1.bed";
tcel_K4ME1_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/anno.matrix.01092016.log2rpm.TCELL.K4ME1.BEDMAP.ALL.K4ME1.bed";

mono_meth_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/meth_20151028_anno_sort.txt";
neut_meth_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/meth_20151028_anno_sort.txt";
tcel_meth_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/meth_20151028_anno_sort.txt";

mono_psi_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/sorted_psi_20151028_anno_sort3.txt";
neut_psi_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/sorted_psi_20151028_anno_sort3.txt";
tcel_psi_pheno="/nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/sorted_psi_20151028_anno_sort3.txt";



Geno_path=/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/$CellType\_$Tag/chrom$Chr.h5;
Stdid_rsid_map="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/Beta_change_info/All_chromosome_beta_change.alt.AF.txt"


Pheno_path=$CellType\_$Tag\_pheno;

#########################################################################################

# echo $QTL_path;
# echo ${!Pheno_path};
# echo $Geno_path;

echo "QTL processing....."
(echo -e "SNP\tphenotypeID\tp.value" && cat $QTL_path | sed '1d'| grep ^$Chr\: | awk '{OFS="\t"; print $1,$3,$4}') > $TEMPDIR/$Chr\_cis.eqtls.txt;

echo "Gen position processing....."
(echo -e "snp\tchr_snp\tpos" && cat $TEMPDIR/$Chr\_cis.eqtls.txt |sed '1d'|awk '{print $1}'|awk -F"_" '{OFS="\t"; print $0,$1}'| awk -F":" '{print $1,$2,$3}'| awk '{print $1":"$2"\t"$3"\t"$4}') > $TEMPDIR/$Chr\_gene.pos.txt;




echo "Pheno position processing....."
# if [ $Tag == "psi" ]; then
#
# 	(echo -e "phenotypeID\tchrom_probe\ts1\ts2" && cat ${!Pheno_path} | sed '1d'| awk '{OFS="\t"; if($2=='"$Chr"') print $22,$2,$3,$4}') > $TEMPDIR/$Chr\_phe.pos.txt;
# else

	(echo -e "phenotypeID\tchrom_probe\ts1\ts2" && cat ${!Pheno_path} | sed '1d'| awk '{OFS="\t"; if($2=='"$Chr"') print $1,$2,$3,$4}') > $TEMPDIR/$Chr\_phe.pos.txt;

# fi


echo "Gen matrix processing....."
$Rscript $GenoMatrixScript $Geno_path "/genotypes/col_headers/rs" "/genotypes/matrix" $TEMPDIR/$Chr\_geno_matrix_pre0.txt;

## Mapping rsID to StdID
awk 'NR==FNR {h[$2] = $1; next} {print h[$1],$0}' $Stdid_rsid_map $TEMPDIR/$Chr\_geno_matrix_pre0.txt | cut -d" " -f1,3- > $TEMPDIR/$Chr\_geno_matrix_pre1.txt
rm $TEMPDIR/$Chr\_geno_matrix_pre0.txt;

## Making the genotype matrix for variants that have been tested (or present in the *cis.eqtls.txt) in Limix. Otherwise EigenMT would not work.
cat $TEMPDIR/$Chr\_cis.eqtls.txt |  sed '1d'| cut -f1| awk '!seen[$0]++' > $TEMPDIR/$Chr\_var.txt;
awk 'NR==FNR {h[$1] = $0; next} {print $0,h[$1]}' $TEMPDIR/$Chr\_geno_matrix_pre1.txt $TEMPDIR/$Chr\_var.txt | cut -d" " -f2- > $TEMPDIR/$Chr\_geno_matrix_pre.txt
rm $TEMPDIR/$Chr\_geno_matrix_pre1.txt $TEMPDIR/$Chr\_var.txt;



Ncol=$(awk '{print NF-1; exit}' $TEMPDIR/$Chr\_geno_matrix_pre.txt);
(echo "ID" && for i in $(seq 1 $Ncol); do echo "SAM$i"; done) | tr "\n" " "| sed '$s/ $/\n/' > $TEMPDIR/$Chr\_header.txt;

cat $TEMPDIR/$Chr\_header.txt $TEMPDIR/$Chr\_geno_matrix_pre.txt > $TEMPDIR/$Chr\_geno_matrix.txt
rm $TEMPDIR/$Chr\_header.txt $TEMPDIR/$Chr\_geno_matrix_pre.txt;



echo "Running EigenMT for Chr $Chr....."

$PYTHON $PY_EigenMT \
--CHROM $Chr \
--QTL $TEMPDIR/$Chr\_cis.eqtls.txt \
--GEN $TEMPDIR/$Chr\_geno_matrix.txt \
--GENPOS $TEMPDIR/$Chr\_gene.pos.txt \
--PHEPOS $TEMPDIR/$Chr\_phe.pos.txt \
--cis_dist 1e6 \
--OUT $TEMPDIR/$Chr\_EigenMT.txt

## Cleanning temp files
rm $TEMPDIR/$Chr\_cis.eqtls.txt  $TEMPDIR/$Chr\_geno_matrix.txt $TEMPDIR/$Chr\_gene.pos.txt $TEMPDIR/$Chr\_phe.pos.txt ;































