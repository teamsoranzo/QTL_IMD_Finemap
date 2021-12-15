import sys
sys.path.append('/nfs/team151_data03/WP10_release/release_V1_Oct2016/LIMIX/pysrc/')

import os



in_dir = sys.argv[1]
out_dir=sys.argv[2]



import glob
combat_files= glob.glob(in_dir+"*.h5")



K_arr=[10]
#K_arr=[5,10, 20, 30, 40]

for my_in_file in combat_files:
    print my_in_file

    file_name = os.path.basename(my_in_file)
    file_root = os.path.splitext(file_name)[0]


    for K in K_arr:
        out_o = os.path.join(out_dir,file_root+"_peer_"+str(K)+".o" )
        out_e = os.path.join(out_dir,file_root+"_peer_"+str(K)+".e" )

        cmd_R=''
        cmd_R += 'bsub -G hematopoiesis -o '+out_o+' -e '+ out_e
        cmd_R += ' -q normal -M 5500 -R"select[mem>=5500] rusage[mem=5500]"'
        cmd_R += ' -- "/software/R-3.1.0/bin/Rscript /nfs/users/nfs_k/kk8/Projects/Scripts/Blueprint/Re_analysis/Limix/eQTL/PEER_3.R ' + my_in_file +' ' +str(K)+'"'+' '+out_dir

        print cmd_R
        os.system(cmd_R)

