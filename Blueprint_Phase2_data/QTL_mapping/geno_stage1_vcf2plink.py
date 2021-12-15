import sys
sys.path.append('/nfs/team151_data03/WP10_release/release_V1_Oct2016/LIMIX/pysrc/')
from CFG.settings import *
import h5py
import pdb
import os

if __name__ == "__main__":

    out_dir = "/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/WGS_plink"
    in_dir ="/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/VCF_QC/Processed_VCF/Final_VCFs"
    group_name = 'blueprint_vcf2plink'
    command = "bgadd -L %d /%s" % (22, group_name)
    print command
    os.system(command)

    #pdb.set_trace()

    chroms = range(1,23)
    for chrom_i in chroms:
        # VCF file
        gzvcf_file = '%d.BPWP10_23_05_17.vcf.gz' % chrom_i
        gzvcf_file = os.path.join(in_dir, gzvcf_file)

        # temp dir
        temp_dir = os.path.join(out_dir,'temp','chrom%d'%chrom_i)
        if not os.path.exists(temp_dir): os.makedirs(temp_dir)

        # sterr and stdout folder
        std_dir = os.path.join(out_dir, 'std')
        if not os.path.exists(std_dir): os.makedirs(std_dir)

        #submit
        command  = 'bsub -G hematopoiesis -R "select[mem>20000] rusage[mem=20000]" -M 20000 '#rusage[tmp=1024]
        command += '-o %s/chrom%d.o ' % (std_dir, chrom_i)
        command += '-e %s/chrom%d.e ' % (std_dir, chrom_i)
        command += '-g /%s /software/vcftools_0.1.12a/bin/vcftools --gzvcf %s ' % (group_name, gzvcf_file)
        command += '--plink --out %s/chrom%d' % (out_dir, chrom_i)
        #command += ' --temp %s' % temp_dir
        #command += ' --maf 0.01'
        print command
        os.system(command)

