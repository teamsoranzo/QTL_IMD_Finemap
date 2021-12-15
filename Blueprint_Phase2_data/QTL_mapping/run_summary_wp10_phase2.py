__author__ = 'yy5'

__author__ = 'yy5'

import sys
sys.path.append('/nfs/team151_data03/WP10_release/release_V1_Oct2016/LIMIX/pysrc/')

import os



in_dir = sys.argv[1]
in_dir_phe=sys.argv[2]




####define genotype files"

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

my_python = "/lustre/scratch114/projects/hematopoiesis/Blueprint/software/python-2.7.9/bin/python"

sub_dirs = os.listdir(in_dir)

############ phenotype files#############

print sub_dirs



for my_dir in sub_dirs:

    #print my_dir

    sum_folder_name = os.path.join(in_dir, my_dir, "runs", "split_folder/")


    my_fp = os.path.join(in_dir_phe, my_dir+".h5" )
    out_o = os.path.join(sum_folder_name, my_dir+".o" )
    out_e = os.path.join(sum_folder_name, my_dir+".e" )
    #print my_fp
    if os.path.exists(sum_folder_name):


        my_cmd=''
        my_cmd += 'bsub -G hematopoiesis -o '+out_o+' -e '+ out_e
        my_cmd += ' -q long -M 2000 -R"select[mem>=2000] rusage[mem=2000]" -- "'
        my_cmd += my_python
        my_cmd += ' /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/Re_analysis/Limix/eQTL/eqtl_cis_summary_new.py --in_dir '+sum_folder_name



        if "mono_gene" in my_dir:
            my_cmd += ' --fg '+mono_gene
            my_cmd += ' --fp '+ my_fp

        if "mono_K27AC" in my_dir:
            my_cmd += ' --fg '+mono_K27AC
            my_cmd += ' --fp '+ my_fp

        if "mono_K4ME1" in my_dir:
            my_cmd += ' --fg '+mono_K4ME1
            my_cmd += ' --fp '+ my_fp

        if "mono_meth" in my_dir:
            my_cmd += ' --fg '+mono_meth
            my_cmd += ' --fp '+ my_fp


        if "mono_psi" in my_dir:
            my_cmd += ' --fg '+mono_psi
            my_cmd += ' --fp '+ my_fp


        if "neut_gene" in my_dir:
            my_cmd += ' --fg '+neut_gene
            my_cmd += ' --fp '+ my_fp

        if "neut_K27AC" in my_dir:
            my_cmd += ' --fg '+neut_K27AC
            my_cmd += ' --fp '+ my_fp

        if "neut_K4ME1" in my_dir:
            my_cmd += ' --fg '+neut_K4ME1
            my_cmd += ' --fp '+ my_fp

        if "neut_meth" in my_dir:
            my_cmd += ' --fg '+neut_meth
            my_cmd += ' --fp '+ my_fp

        if "neut_psi" in my_dir:
            my_cmd += ' --fg '+neut_psi
            my_cmd += ' --fp '+ my_fp


        if "tcel_gene" in my_dir:
            my_cmd += ' --fg '+tcel_gene
            my_cmd += ' --fp '+ my_fp

        if "tcel_K27AC" in my_dir:
            my_cmd += ' --fg '+tcel_K27AC
            my_cmd += ' --fp '+ my_fp

        if "tcel_K4ME1" in my_dir:
            my_cmd += ' --fg '+tcel_K4ME1
            my_cmd += ' --fp '+ my_fp

        if "tcel_meth" in my_dir:
            my_cmd += ' --fg '+tcel_meth
            my_cmd += ' --fp '+ my_fp

        if "tcel_psi" in my_dir:
            my_cmd += ' --fg '+tcel_psi
            my_cmd += ' --fp '+ my_fp


        my_cmd += '"'
        print my_cmd

        my_cmd_all = my_cmd.replace("/nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/Re_analysis/Limix/eQTL/eqtl_cis_summary_new.py", "/nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/Re_analysis/Limix/eQTL/eqtl_cis_summary_output_all.py")
        print my_cmd_all
        os.system(my_cmd_all)
        os.system(my_cmd)
        #if "peer_10" in my_dir:
            #os.system(my_cmd_all)

