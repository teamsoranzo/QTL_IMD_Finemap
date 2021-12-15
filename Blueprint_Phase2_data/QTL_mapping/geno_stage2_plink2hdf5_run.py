import sys
import pdb
import os
import scipy as SP

n_chroms = 22

if __name__ == "__main__":

    # CREATION OF THE GROUP
    group_name = 'blueprint_plink2hdf5_phase2'
    print ""
    command = "bgadd /%s" % group_name
    print command
    os.system(command)
    command = "bgmod -L %d /%s" % (n_chroms,group_name)
    print command
    os.system(command)
    print ""

    #Create temp dir
    temp_folder = '/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/WGS_hdf5/temp/%s' %group_name
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

#     pdb.set_trace()

    # GO! GO! GO!
    ## Note that for Chr 3, ~85g memory is required. Change accordingly.
    chroms = range(1,n_chroms+1)
    for chrom in chroms:
        stdout_file = os.path.join(temp_folder,'stdout_%d.o'%chrom)
        stderr_file = os.path.join(temp_folder,'stderr_%d.e'%chrom)
        command  = "bsub  -G hematopoiesis -q normal -R \"select[mem>90000] rusage[mem=90000]\" -M 90000 -g /%s " % group_name
        command += "-o %s " % stdout_file
        command += "-e %s " % stderr_file
        command += "/lustre/scratch114/projects/hematopoiesis/Blueprint/software/python-2.7.9/bin/python /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/Re_analysis/Limix/preprocess/geno_stage2_plink2hdf5.py %d" % chrom
        print command
        os.system(command)

