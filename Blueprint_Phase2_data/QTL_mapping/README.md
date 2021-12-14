## QTL mapping by Limix pipeline

########### before you run the pipeline, please bear the following in mind ###################\
you should redirect the path to specific location to run each code.
you should create three main output folders before hand: eqtl_result  geno_file  phenotype_file.
b4_peer and after_peer two sub folders were created under phenotype_file
#############################################################################################

shortcuts \
`lpython = /lustre/scratch114/projects/hematopoiesis/Blueprint/software/python-2.7.9/bin/python`

#############################################################################################


---

#### PREPROCESSING STEPS



1. **Filter the count matrix of phenotype**


    filtering genes of the following criteria:
    1. on autosomal chromosome
    2. expressed (expr>THRESHOLD_EXPR) in at least half of the samples
       ((expr>THRESHOLD_EXPR).mean(0)>THRESHOLD_MEAN)


inputs: 
- count.bed
           read counts matrix in rows [features] X cols [donors] \
            `Ex. /lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/lv3/WP10X/BEDMAP.NEU.PU1/limix.FINAL.matrix.24112016.log2rpm.NEU.PU1.BEDMAP.NEU.PU1.bed`

- anno.bed
            annotation file of the features:   feature_name    chr    start    end

```
    cd pysrc/preprocess/
    python expr_stage1_csv2hdf5.py -h


    Usage: expStage1_cv2hdf5 [Options] --expr_file ..

Options:
  -h, --help            show this help message and exit
  -f EXPRESSION_FILENAME, --expr_file=EXPRESSION_FILENAME
                        your input expression file path
  -l LOC_FILENAME, --loc_file=LOC_FILENAME
                        your input gene location file path
  -o H5_FILENAME, --outfile=H5_FILENAME
                        your output h5 file path, please name it with sufix
                        expr.h5
  -t THRESHOLD_EXPR, --threshold_expr=THRESHOLD_EXPR
                        threshold for filtering expression read
  -m THRESHOLD_MEAN, --threshold_mean=THRESHOLD_MEAN
                        threshold for filtering mean expression read

```

`$lpython expr_stage1_csv2hdf5.py --expr_file <count.bed> --loc_file <anno.bed> --outfile <b4_peer/count.h5> --threshold_expr <0> --threshold_mean <0.5>`
    


**IMPORTANT:**

Gene: --threshold_expr 3.321928 (log2(10))--threshold_mean 0.50 \
PSI: --threshold_expr 0.1 --threshold_mean 0.50 \
K27AC: --threshold_expr 0 --threshold_mean 0.5 \
K4ME1: --threshold_expr 0 --threshold_mean 0.5 \
Meth: --threshold_expr -15 --threshold_mean 0.0 (we wanted to check all the methylation probe, since methylated and unmethylated both sites are important) - Used M value and NOT beta value (which is a percentage value) 

Others:
Exon: --threshold_expr 2.321928 (log2(5))--threshold_mean 0.50 \
SJ: --threshold_expr 2.321928 --threshold_mean 0.50 


**2. converts vcf geno file to plink chrom-by-chrom
    (takes chrom num as argument)**

```
FIXED arguments
out_dir = "$path/limix/data_sep24/WGS/filter/WGS_plink"
in_dir ="$path/limix/data_sep24/WGS/filter/vcf_files"
```

```
cd pysrc/preprocess/
geno_stage1_vcf2plink.py  —> EDIT
```

    
    
**3. convert plink geno file to hdf5 chrom-by-chrom
    (takes chrom num as argument)**

```
FIXED arguments
 temp_folder = '$path/limix/data_sep24/WGS/filter/WGS_hdf5/temp/%s' %group_name

    pysrc/preprocess/
geno_stage2_plink2hdf5_run.py  —> EDIT

depends on
FIXED arguments
    in_dir = "$path/limix/data_sep24/WGS/filter/WGS_plink"
    out_dir = "$path/limix/data_sep24/WGS/filter/WGS_hdf5"
```

```
cd pysrc/preprocess/
geno_stage2_plink2hdf5.py  —> EDIT
```

**NOTE:** \
**Effect allele based on the 1st allele of  donor 1**

```
call to function readPED is from line 24
geno_stage2_plink2hdf5.py:    RV = readPED(bfn, standardize = False)

readPED function is lines 70,77 from this script below.
$path/software/python-2.7.9/lib/python2.7/site-packages/limix-0.7.3-py2.7-linux-x86_64.egg/limix/io/plink.py
 
copied from plink.py – it is  vals[0,0]
 
        for i in xrange(snpsstr.shape[1]/2):
            snps[inan[:,2*i],i]=SP.nan
            vals=snpsstr[~inan[:,2*i],2*i:2*(i+1)]
            snps[~inan[:,2*i],i]+=(vals==vals[0,0]).sum(1)

```


**4. Filter the hfd5 genotype files**

Filtering steps for variants
- max of 5% missingness
- min of 4% MAF in the filtered donors
- sd > 0 across filtered donors
- impute missing Geno to mean of donors

```
FIXED arguments
maf_min = 0.04
in_dir = "$path/limix/data_sep24/WGS/filter/WGS_hdf5


inputs
- <dir_vcf_hfd5>
Ex $path/limix/data_sep24/WGS/filter/WGS_hdf5/
- <odir_filt_vcf>
Ex $path/wp10_ex_pu1_limix/geno_file/NEU_PU1/
    pysrc/preprocess/
```

```
$python geno_stage3_filter_run.py <dir_vcf_hfd5> <b4_peer/count.h5> <odir_filt_vcf>
$python geno_stage4_gather.py <odir_filt_vcf>
 ``` 
 
---
#### PEER

**5. run peer on count matrix in step 1**
```
cd pysrc/eqtl
 $python run_peer_combat_in_bulk <in_dir> <out_dir>
```
 
if the peer factor k need to be change, please modify line 18:
for example, if you want to have two sets of results, one set with 10 factors, another set has 20 factors, you could change the k_arr[]
K_arr=[10]--> k_arr=[10,20]
    
---

#### LIMIX

**6. run eqtl**

inputs:
- <indir> $path/wp10_ex_pu1_limix/phenotype_file/after_peer/
- <outdir> $path/wp10_ex_pu1_limix/eqtl_result/
- <window_size> 2500

```
cd pysrc/eqtl
$python run_eqtl_cis_run_in_bulk_wp10x.py <indir> <outdir> <window_size>   —> EDIT
```

Edit
filtered genotype hdf5 file
NEU_PU1="$path/wp10_ex_pu1_limix/geno_file/NEU_PU1/all_chrom.h5"


- To change the cis window, you could simply change the third argument from "run eqtl" command:\
for example: \
    `lpython run_eqtl_cis_run_in_bulk_wp10x.py $path/wp10_ex_pu1_limix/phenotype_file/after_peer/ $path/wp10_ex_pu1_limix/eqtl_result/ 2500`
- To change to 5k window, you should just do:
for example: \
    `lpython run_eqtl_cis_run_in_bulk_wp10x.py $path/wp10_ex_pu1_limix/phenotype_file/after_peer/ $path/wp10_ex_pu1_limix/eqtl_result/ 5000`

You should also check if the eqtl run successfully, so run the following to check if job has gone through:\
`python checkError.py <outdir>`
if one of the job failed, you should delete the output folder, and rerun the job. 
For example if the a job under NEU_PU1_peer_10(under eqtl_result) failed, you should delete the entire NEU_PU1_peer_10, and put the corresponding phenotype file in a folder to rerun this single job.

`eqtl_cis.py`


**7. run summary**
code is the same folder as eqtl($path/wp10extlimix/pysrc/eqtl)
lpython run_summary_wp10x.py $path/wp10_ex_pu1_limix/eqtl_result_2500win/ after_peer/

name the summary with trait/qtl type
lpython rename_all_summary.py $path/wp10_ex_pu1_limix/eqtl_result/


CHANGE:
in line 37 — eqtl_cis_summary_output_all.py / eqtl_cis_summary_new.py
CHANGE to your window size, here it is 2.5KB —>     data = QtlData(my_fg,my_fp,2500)


code is the same folder as eqtl($path/wp10extlimix/pysrc/eqtl)
lpython run_summary_wp10x.py $path/wp10_ex_pu1_limix/eqtl_result_2500win/ after_peer/

name the summary with trait/qtl type
lpython rename_all_summary.py $path/wp10_ex_pu1_limix/eqtl_result/

This script will generate All summary (in txt format) and simple summary (in hdf5 format) files.

to generate all summary table: /nfs/team151_data03/WP10_release/release_V1_Oct2016/LIMIX/pysrc/eqtl/eqtl_cis_summary_output_all.py

*OUTPUT table format
1          2    3     4          5      6            7
geneID rs pos raw_pv Beta Bonf_pv local_FDR


to generate simple table         : /nfs/team151_data03/WP10_release/release_V1_Oct2016/LIMIX/pysrc/eqtl/eqtl_cis_summary_new.py
                                                   the simple table generated here in hdf5 format will have a global FDR value with header 'qv_all'
The output simple table is in hdf5 format.

To produce a text file from hdf5 simple table: 
/nfs/team151_data03/WP10_release/release_V1_Oct2016/LIMIX/get.table.from.simple.hdf5.R <simple.hdf5> <output_filename>
This R script will produce a simple table with this format:
*OUTPUT table format
1   2          3                    4               5              6     7      8   9            10            11
rs geneID gene_chrom gene_start gene_end pos beta pv pv_bonf pv_storey qv_all
>> pv_storey is local FDR
>> qv_all is global FDR
---
8. Properly format summary statistics AND change sign of Beta

    chr:pos_ref_alt rsid phenotypeID p.value beta Bonferroni.p lFDR alt.AF std.error gFDR
    20:48969024_A_G rs6125980 ENSG00000000419.8 2.953e-05 -0.7903 0.120613960675 1.187e-01 0.0725 0.1892 0.263348426391099

RUN:
bash ~/Projects/Scripts/Blueprint/Re_analysis/Limix/summary_format_and_beta_change.sh <ID> <infile> <Beta_key_file>
## Change the Beta sign
#echo 'bash ~/Projects/Scripts/Blueprint/Re_analysis/Limix/summary_format_and_beta_change.sh ${LSB_JOBINDEX} infile.txt' | bsub -G hematopoiesis -J "Beta_change[1-22]" -o JOB/%J.%I.beta_change.o -q normal -R "select[mem>=2000] rusage[mem=2000] span[hosts=1]" -n1 -M2000;

Calculate SE:
Rscript ~/Projects/Scripts/Blueprint/Re_analysis/Limix/SE_from_beta_and_pval.R  <infile> <outfile>
Example:
## Calculate the SE
#Rscript ~/Projects/Scripts/Blueprint/Re_analysis/Limix/SE_from_beta_and_pval.R  qtl_result/mono_psi_nor_31052017.txt_sort_peer_10/runs/split_folder/mono_psi_nor_31052017.txt_sort_peer_10_summary.Beta_changed.txt qtl_result/mono_psi_nor_31052017.txt_sort_peer_10/runs/split_folder/mono_psi_nor_31052017.txt_sort_peer_10_summary.Beta_changed.SE.txt



9. EigenMT  https://github.com/joed3/eigenMT
/nfs/team151_data03/WP10_release/release_V1_Oct2016/LIMIX/eigenMT.py


*must install scikit
pip install --user -U scikit-learn

*For updated version of numpy
cd /nfs/team151/kk8/Py_LIB/numpy-1.11.2/
python setup.py build -j 4 install --prefix $HOME/.local

Before running - 
a)    EigenMT should be run on chromosome-wise 
b)    The variance present in the genotype matrix (genotype.txt) should also be present on the qtl file (cis.eqtls) 


python eigenMT.py \
        --CHROM 19 \
        --QTL example/cis.eqtls.txt \
        --GEN example/genotypes.txt \
        --GENPOS example/gen.positions.txt \
        --PHEPOS example/phe.positions.txt \
        --OUT example/exampleOut.txt
        —cis_dist 1e6



** cist_dist is from peak TSS

Used Default values:
--var_thresh Default is 99% variance explained.
--window Default is 200 SNPs.
 --external            indicates whether the provided genotype matrix is
                        different from the one used to call cis-eQTLs
                        initially (default = False)

- - QTL    
    simple format of  three columns: 1. the SNP ID, 2. the gene ID, and 3. the cis-eQTL p-value. This file must contain a header line with the p-value column indicated by p-value.

RUN:
bash EigenMT_pval_calculation.sh <chr> <mono> <gene>
RUN EigenMT -
#echo 'bash EigenMT_pval_calculation.sh ${LSB_JOBINDEX} tcel meth' | bsub -G hematopoiesis -J "eigenMT[1-22]" -o JOB/%J.%I.tcel.meth.o -q normal -R "select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -M4000;


10. Storey FDR collection on adjusted p-values from step 9 —> Calculate global FDR on EigenMT corrected P-value to be added to simple table
Use R function:
Rscript Rscript.Storey.Qval.R <input> <output>
Run manually for all Celltypes and features (total 15) - it was quick !!

Input file: EigenMT output file


11. Calculate the Eigen.pval 


## Eigen P.value merge (simple table)
# (echo "chr.pos_ref_alt rsid phenotypeID p.value beta Bonferroni.p lFDR alt.AF se Eigen.p gFDR" && awk 'NR==FNR {h[$2] = $0; next} {print $0,h[$3]}'  EigenMT/mono_K4ME1_all_EigenMT_gFDR.txt qtl_result/mono_K4ME1_log2rpm_31052017_bedmap_peer_10/runs/split_folder/mono_K4ME1_log2rpm_31052017_bedmap_peer_10_summary.Beta_changed.SE.txt |sed '1d'| awk '{if($3 ==$12) print $1,$2,$3,$4,$5,$6,$7,$8,$10,$14,$16}') > qtl_result/mono_K4ME1_log2rpm_31052017_bedmap_peer_10/runs/split_folder/mono_K4ME1_log2rpm_31052017_bedmap_peer_10_summary.Beta_changed.SE.Eigen.pval.txt
    
    









