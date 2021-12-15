__author__ = 'yy5'

import sys
sys.path.append('/nfs/team151_data03/WP10_release/release_V1_Oct2016/LIMIX/pysrc/')

import os



in_dir = sys.argv[1]

#main_out_dir="/lustre/scratch114/projects/hematopoiesis/Blueprint/limix/data_sep24/RNA-seq/outputs/eqtl_out/combat_nor_result/"
main_out_dir=sys.argv[2]

window_size = sys.argv[3]


####define genotype files"
####### wp10 re_analysis ########

mono_gene="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/mono_gene/all_chrom.h5"
mono_K27AC="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/mono_K27AC/all_chrom.h5"
mono_K4ME1="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/mono_K4ME1/all_chrom.h5"
mono_meth="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/mono_meth/all_chrom.h5"
mono_psi="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/mono_psi/all_chrom.h5"

neut_gene="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/neut_gene/all_chrom.h5"
neut_K27AC="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/neut_K27AC/all_chrom.h5"
neut_K4ME1="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/neut_K4ME1/all_chrom.h5"
neut_meth="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/neut_meth/all_chrom.h5"
neut_psi="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/neut_psi/all_chrom.h5"

tcel_gene="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/tcel_gene/all_chrom.h5"
tcel_K27AC="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/tcel_K27AC/all_chrom.h5"
tcel_K4ME1="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/tcel_K4ME1/all_chrom.h5"
tcel_meth="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/tcel_meth/all_chrom.h5"
tcel_psi="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/tcel_psi/all_chrom.h5"


#NEU_PU1="/lustre/scratch115/projects/hematopoiesis/wp10_ex_pu1_limix/geno_file/NEU_PU1/all_chrom.h5"
#NEU_CEBPB="/lustre/scratch115/projects/hematopoiesis/wp10_extension/geno_file/NEU_CEBPB/all_chrom.h5"
#NEU_K27ME3_on_K27ME3="/lustre/scratch115/projects/hematopoiesis/wp10_extension/geno_file/NEU_K27ME3_on_K27ME3/all_chrom.h5"
#NEU_K27ME3_on_K4ME3="/lustre/scratch115/projects/hematopoiesis/wp10_extension/geno_file/NEU_K27ME3_on_K4ME3/all_chrom.h5"
#NEU_K4ME3="/lustre/scratch115/projects/hematopoiesis/wp10_extension/geno_file/NEU_K4ME3/all_chrom.h5"
#TCELL="/lustre/scratch115/projects/hematopoiesis/wp10_main_tcell_limix/genotype/TCELL_K27AC_ALL_K27AC/all_chrom.h5"
#NEU_K27AC_CHOP="/lustre/scratch114/projects/hematopoiesis/Blueprint/new_eqtl_chip_seq/genotype/NEU_K27AC_CHOP500_K27AC_peer_10/all_chrom.h5"
#NEU_K27AC_ALL="/lustre/scratch114/projects/hematopoiesis/Blueprint/new_eqtl_chip_seq/genotype/NEU_K27AC_ALL_K27AC_peer_10/all_chrom.h5"
#MON_K27AC_CHOP="/lustre/scratch114/projects/hematopoiesis/Blueprint/new_eqtl_chip_seq/genotype/MON_K27AC_CHOP500_K27AC_peer_10/all_chrom.h5"
#MON_K27AC_ALL="/lustre/scratch114/projects/hematopoiesis/Blueprint/new_eqtl_chip_seq/genotype/MON_K27AC_ALL_K27AC_peer_10/all_chrom.h5"
#SC_neut_all_k27ac="/lustre/scratch114/projects/hematopoiesis/Blueprint/cross_sample_eqtl/genotype/SC_neut_all_k27ac/all_chrom.h5"
#SC_mon_all_k27ac="/lustre/scratch114/projects/hematopoiesis/Blueprint/cross_sample_eqtl/genotype/SC_mon_all_k27ac/all_chrom.h5"
#NJ_neut_all_k27ac="/lustre/scratch114/projects/hematopoiesis/Blueprint/cross_sample_eqtl/genotype/NJ_neut_all_k27ac/all_chrom.h5"
#NJ_mon_all_k27ac="/lustre/scratch114/projects/hematopoiesis/Blueprint/cross_sample_eqtl/genotype/NJ_mon_all_k27ac/all_chrom.h5"
#May_MON_K4ME1="/lustre/scratch114/projects/hematopoiesis/Blueprint/limix/data_sep24/hdf5/me1/May_MON_K4ME1/all_chrom.h5"
#May_NEU_K4ME1="/lustre/scratch114/projects/hematopoiesis/Blueprint/limix/data_sep24/hdf5/me1/May_NEU_K4ME1/all_chrom.h5"
#May_TCELL_K4ME1="/lustre/scratch114/projects/hematopoiesis/Blueprint/limix/data_sep24/hdf5/me1/May_TCELL_K4ME1/all_chrom.h5"
#May_TCELL_K27AC="/lustre/scratch114/projects/hematopoiesis/Blueprint/limix/data_sep24/hdf5/me1/May_TCELL_K27AC/all_chrom.h5"

#CTCF="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/yy5/CTCF/geno_file/all_chrom.h5"




my_python = "/lustre/scratch114/projects/hematopoiesis/Blueprint/software/python-2.7.9/bin/python"

import glob
peer_files= glob.glob(in_dir+"*.h5")



for my_in_file in peer_files:

    file_name = os.path.basename(my_in_file)

    folder_name = os.path.splitext(file_name)[0]

    my_out_dir = os.path.join(main_out_dir, folder_name)
    print my_out_dir


    if not os.path.exists(my_out_dir):
            os.makedirs(my_out_dir)

    my_cmd =''
    my_cmd += my_python
    my_cmd += ' /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/Re_analysis/Limix/eQTL/eqtl_cis_run.py --out_dir '+my_out_dir+'/'

    if "mono_gene" in file_name:
        my_cmd += ' --fg '+mono_gene
        my_cmd += ' --fp '+ my_in_file

    if "mono_K27AC" in file_name:
        my_cmd += ' --fg '+mono_K27AC
        my_cmd += ' --fp '+ my_in_file

    if "mono_K4ME1" in file_name:
        my_cmd += ' --fg '+mono_K4ME1
        my_cmd += ' --fp '+ my_in_file

    if "mono_meth" in file_name:
        my_cmd += ' --fg '+mono_meth
        my_cmd += ' --fp '+ my_in_file

    if "mono_psi" in file_name:
        my_cmd += ' --fg '+mono_psi
        my_cmd += ' --fp '+ my_in_file

    if "neut_gene" in file_name:
        my_cmd += ' --fg '+neut_gene
        my_cmd += ' --fp '+ my_in_file

    if "neut_K27AC" in file_name:
        my_cmd += ' --fg '+neut_K27AC
        my_cmd += ' --fp '+ my_in_file

    if "neut_K4ME1" in file_name:
        my_cmd += ' --fg '+neut_K4ME1
        my_cmd += ' --fp '+ my_in_file

    if "neut_meth" in file_name:
        my_cmd += ' --fg '+neut_meth
        my_cmd += ' --fp '+ my_in_file

    if "neut_psi" in file_name:
        my_cmd += ' --fg '+neut_psi
        my_cmd += ' --fp '+ my_in_file

    if "tcel_gene" in file_name:
        my_cmd += ' --fg '+tcel_gene
        my_cmd += ' --fp '+ my_in_file

    if "tcel_K27AC" in file_name:
        my_cmd += ' --fg '+tcel_K27AC
        my_cmd += ' --fp '+ my_in_file

    if "tcel_K4ME1" in file_name:
        my_cmd += ' --fg '+tcel_K4ME1
        my_cmd += ' --fp '+ my_in_file

    if "tcel_meth" in file_name:
        my_cmd += ' --fg '+tcel_meth
        my_cmd += ' --fp '+ my_in_file

    if "tcel_psi" in file_name:
        my_cmd += ' --fg '+tcel_psi
        my_cmd += ' --fp '+ my_in_file


    my_cmd += ' --my_window '+ window_size
    print my_cmd
    os.system(my_cmd)

