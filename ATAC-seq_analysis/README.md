### The ATAC-Seq ranking method starts with the per cell normalisation of the raw ATAC-Seq signal.

$/software/R-3.6.1/bin/Rscript ~/Scripts/R/133__ATAC_normalize_and_prepare_Kousik.R  --ATACSeq_1 /lustre/scratch119/humgen/teams/soranzo/users/mt19/web\
tool_atac_variants.txt --ATACSeq_2 /lustre/scratch119/humgen/teams/soranzo/users/mt19/webtool_atac_variants_remainder.txt  --type Kousik_overlap --out /lust\
re/scratch119/humgen/teams/soranzo/users/mt19/

### Afterwards we obtain the overlap for the IMD, QTLs and sahred QTL variants and restrict ourselves to the ATAC-Seq signals in the CD4 T-cells and the monocyte lineages

$/software/R-3.6.1/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/180_ATAC_I_RESULT_Kousik.R --subset Variant_set_final.txt --ATAC_Scaled ATAC_scaled_relevant_t\
issue_info_Kousik_overlap.txt --ATAC_Raw ATAC_RAW_Kousik_overlap.txt --Lineage_hierarchy /lustre/scratch115/teams/soranzo/projects/ALL_dB/ATAC_scaled_Lineag\
e_hierarchy_Kousik_QTL_generation.txt --trait_table  /lustre/scratch115/teams/soranzo/projects/ALL_dB/ATAC_scaled_trait_table_Kousik_QTL_generation.txt --ty\
pe Kousik_overlap --out /nfs/users/nfs_m/mt19/ATAC_Kousik/

### Finally we used the package desiR(Lazic ) to rank the signal using two criteria: (i) the maximum normalized signal in a relevant cell type and (ii) the number of relevant cell types with a normalized ATAC-Seq signal above 0.5. We gave three times more relative weight to the maximum normalised ATAC-Seq signal in a cell type.
$/software/R-3.6.1/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/182_ATAC_III_summary_measure_Kousik.R --subset Variant_set_final.txt --ATAC_RESU\
LT ATAC_RESULT_Kousik_overlap.txt --parameters 0.5,4,0.05,0.8,1,3 --Freq_Threshold 0.5 --type Kousik_overlap --out /nfs/users/nfs_m/mt19/ATAC_Kousik/
