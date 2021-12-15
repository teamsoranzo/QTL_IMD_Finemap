import sys
import pdb
import os
import scipy as SP

n_chroms = 22

if __name__ == "__main__":

    # CREATION OF THE GROUP
    group_name = 'blueprint_geno_filter'
    print ""
    command = "bgadd /%s" % group_name
    print command
    os.system(command)
    command = "bgmod -L %d /%s" % (n_chroms,group_name)
    print command
    os.system(command)
    print ""

    #Create temp dir
    geno_dir=sys.argv[1]
    e_file = sys.argv[2]
    out_dir= sys.argv[3]
    #in_file_e = "/lustre/scratch114/projects/hematopoiesis/Blueprint/limix/data_sep24/RNA-seq/outputs/cl10_hdf5/neut_gene.expr.h5"
    #temp_folder  = '/lustre/scratch114/projects/hematopoiesis/Blueprint/limix/data_sep24/WGS/filter/WGS_hdf5_filt/neut.geno/'+'temp/%s' % group_name
    temp_folder  = out_dir +'temp/%s' % group_name
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    #pdb.set_trace()

    # GO! GO! GO!
    chroms = range(1,n_chroms+1)
    for chrom in chroms:
        stdout_file = os.path.join(temp_folder,'stdout_%d.o'%chrom)
        stderr_file = os.path.join(temp_folder,'stderr_%d.e'%chrom)
        command  = "bsub -G hematopoiesis -q normal -R \"select[mem>3000] rusage[mem=3000]\" -M 3000 -g /%s " % group_name
        command += "-o %s " % stdout_file
        command += "-e %s " % stderr_file
        command += "/lustre/scratch114/projects/hematopoiesis/Blueprint/software/python-2.7.9/bin/python /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/Re_analysis/Limix/preprocess/geno_stage3_filter.py " + str(chrom) +" "+e_file+" "+out_dir
        print command
        os.system(command)

