
ID=$1;
infile=$2;
beta=$3;

echo "$ID $infile $beta";

FILE=$(cat $infile|awk '{if($1 == "'$ID'") print $2}');

# For all summary stat files
# infile -
# 1	qtl_result/mono_gene_nor_combat_31052017.txt_sort_peer_10/runs/split_folder/mono_gene_nor_combat_31052017.txt_sort_peer_10_all_summary.txt
# 2	qtl_result/mono_K27AC_log2rpm_31052017_bedmap_peer_10/runs/split_folder/mono_K27AC_log2rpm_31052017_bedmap_peer_10_all_summary.txt
# 3	qtl_result/mono_K4ME1_log2rpm_31052017_bedmap_peer_10/runs/split_folder/mono_K4ME1_log2rpm_31052017_bedmap_peer_10_all_summary.txt
# 4	qtl_result/mono_meth_M_31052017.txt_sort_peer_10/runs/split_folder/mono_meth_M_31052017.txt_sort_peer_10_all_summary.txt
# 5	qtl_result/mono_psi_nor_31052017.txt_sort_peer_10/runs/split_folder/mono_psi_nor_31052017.txt_sort_peer_10_all_summary.txt
# 6	qtl_result/neut_gene_nor_combat_31052017.txt_sort_peer_10/runs/split_folder/neut_gene_nor_combat_31052017.txt_sort_peer_10_all_summary.txt
# 7	qtl_result/neut_K27AC_log2rpm_31052017_bedmap_peer_10/runs/split_folder/neut_K27AC_log2rpm_31052017_bedmap_peer_10_all_summary.txt
# 8	qtl_result/neut_K4ME1_log2rpm_31052017_bedmap_peer_10/runs/split_folder/neut_K4ME1_log2rpm_31052017_bedmap_peer_10_all_summary.txt
# 9	qtl_result/neut_meth_M_31052017.txt_sort_peer_10/runs/split_folder/neut_meth_M_31052017.txt_sort_peer_10_all_summary.txt
# 10	qtl_result/neut_psi_nor_31052017.txt_sort_peer_10/runs/split_folder/neut_psi_nor_31052017.txt_sort_peer_10_all_summary.txt
# 11	qtl_result/tcel_gene_nor_combat_31052017.txt_sort_peer_10/runs/split_folder/tcel_gene_nor_combat_31052017.txt_sort_peer_10_all_summary.txt
# 12	qtl_result/tcel_K27AC_log2rpm_31052017_bedmap_peer_10/runs/split_folder/tcel_K27AC_log2rpm_31052017_bedmap_peer_10_all_summary.txt
# 13	qtl_result/tcel_K4ME1_log2rpm_31052017_bedmap_peer_10/runs/split_folder/tcel_K4ME1_log2rpm_31052017_bedmap_peer_10_all_summary.txt
# 14	qtl_result/tcel_meth_M_31052017.txt_sort_peer_10/runs/split_folder/tcel_meth_M_31052017.txt_sort_peer_10_all_summary.txt
# 15	qtl_result/tcel_psi_nor_31052017.txt_sort_peer_10/runs/split_folder/tcel_psi_nor_31052017.txt_sort_peer_10_all_summary.txt
####################################################################################################################################################

# echo -e "chr:pos_ref_alt rsid phenotypeID p.value beta Bonferroni.p lFDR alt.AF" >  ${FILE%.*}.Beta_changed.txt;
# awk 'NR==FNR {h[$2] = $0; next} {print $0,h[$2]}' $beta $FILE | awk '{if($2==$9) print $8,$2,$1,$4,$5*$11,$6,$7,$12}' >> ${FILE%.*}.Beta_changed.txt;

####################################################################################################################################################
# For simple summary stat file -
# infile -
# 1	qtl_result/mono_gene_nor_combat_31052017.txt_sort_peer_10/runs/split_folder/mono_gene_nor_combat_31052017.txt_sort_peer_10_summary.txt
# 2	qtl_result/mono_K27AC_log2rpm_31052017_bedmap_peer_10/runs/split_folder/mono_K27AC_log2rpm_31052017_bedmap_peer_10_summary.txt
# 3	qtl_result/mono_K4ME1_log2rpm_31052017_bedmap_peer_10/runs/split_folder/mono_K4ME1_log2rpm_31052017_bedmap_peer_10_summary.txt
# 4	qtl_result/mono_meth_M_31052017.txt_sort_peer_10/runs/split_folder/mono_meth_M_31052017.txt_sort_peer_10_summary.txt
# 5	qtl_result/mono_psi_nor_31052017.txt_sort_peer_10/runs/split_folder/mono_psi_nor_31052017.txt_sort_peer_10_summary.txt
# 6	qtl_result/neut_gene_nor_combat_31052017.txt_sort_peer_10/runs/split_folder/neut_gene_nor_combat_31052017.txt_sort_peer_10_summary.txt
# 7	qtl_result/neut_K27AC_log2rpm_31052017_bedmap_peer_10/runs/split_folder/neut_K27AC_log2rpm_31052017_bedmap_peer_10_summary.txt
# 8	qtl_result/neut_K4ME1_log2rpm_31052017_bedmap_peer_10/runs/split_folder/neut_K4ME1_log2rpm_31052017_bedmap_peer_10_summary.txt
# 9	qtl_result/neut_meth_M_31052017.txt_sort_peer_10/runs/split_folder/neut_meth_M_31052017.txt_sort_peer_10_summary.txt
# 10	qtl_result/neut_psi_nor_31052017.txt_sort_peer_10/runs/split_folder/neut_psi_nor_31052017.txt_sort_peer_10_summary.txt
# 11	qtl_result/tcel_gene_nor_combat_31052017.txt_sort_peer_10/runs/split_folder/tcel_gene_nor_combat_31052017.txt_sort_peer_10_summary.txt
# 12	qtl_result/tcel_K27AC_log2rpm_31052017_bedmap_peer_10/runs/split_folder/tcel_K27AC_log2rpm_31052017_bedmap_peer_10_summary.txt
# 13	qtl_result/tcel_K4ME1_log2rpm_31052017_bedmap_peer_10/runs/split_folder/tcel_K4ME1_log2rpm_31052017_bedmap_peer_10_summary.txt
# 14	qtl_result/tcel_meth_M_31052017.txt_sort_peer_10/runs/split_folder/tcel_meth_M_31052017.txt_sort_peer_10_summary.txt
# 15	qtl_result/tcel_psi_nor_31052017.txt_sort_peer_10/runs/split_folder/tcel_psi_nor_31052017.txt_sort_peer_10_summary.txt
####################################################################################################################################################
## Note that too get 15 decimal floating value, I have to use "CONVFMT=%.15g"

echo -e "chr:pos_ref_alt rsid phenotypeID p.value beta Bonferroni.p lFDR alt.AF gFDR" >  ${FILE%.*}.Beta_changed.txt;
awk 'NR==FNR {h[$2] = $0; next} {print $0,h[$1]}' $beta $FILE | awk '{if($1==$13) print $12,$1,$2,$8,$7,$9,$10,$16,$11,$15}'| awk -v CONVFMT=%.15g '{gsub($5, $5*$10, $5)}; {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' >> ${FILE%.*}.Beta_changed.txt;

####################################################################################################################################################


## Phenotype id change after beta sign correction ONLY for PSI
# (echo -e "chr:pos_ref_alt rsid phenotypeID p.value beta Bonferroni.p lFDR alt.AF gFDR" && awk 'NR==FNR {h[$1] = $0; next} {print $0,h[$3]}' /nfs/team151_data03/WP10_release/release_V1_Oct2016/PHENOTYPE/sorted_psi_20151028_anno_sort3.txt  qtl_result/mono_psi_nor_31052017.txt_sort_peer_10/runs/split_folder/mono_psi_nor_31052017.txt_sort_peer_10_summary.Beta_changed.txt | awk '{if($3==$10) print $1,$2,$31,$4,$5,$6,$7,$8,$9}') > qtl_result/mono_psi_nor_31052017.txt_sort_peer_10/runs/split_folder/mono_psi_nor_31052017.txt_sort_peer_10_summary.Beta_changed1.txt









