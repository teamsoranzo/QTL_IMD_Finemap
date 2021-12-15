__author__ = 'yy5'


__author__ = 'yy5'
import os,sys


in_dir=sys.argv[1]

for f in os.listdir(in_dir):
    my_in_dir_1 = os.path.join(in_dir,f)
    for f1 in os.listdir(my_in_dir_1):
        if "runs" in f1:
            my_in_dir_2= os.path.join(my_in_dir_1,f1)
            old_name=''
            for f2 in os.listdir(my_in_dir_2+"/split_folder"):

                 if "summary.hdf5" in f2 and "perm" in in_dir:
                     new_name = my_in_dir_2+"/split_folder/"+f+"_summary_perm.hdf5"
                     old_name=my_in_dir_2+"/split_folder/"+f2
                     os.rename(old_name, new_name)

                 if "summary.hdf5" in f2 and "perm" not in in_dir:
                     old_name=my_in_dir_2+"/split_folder/"+f2
                     new_name = my_in_dir_2+"/split_folder/"+f+"_summary.hdf5"
                     os.rename(old_name, new_name)

                 if "all_summary.txt" in f2 and "perm" in in_dir:
                     new_name = my_in_dir_2+"/split_folder/"+f+"_all_summary.txt"
                     old_name=my_in_dir_2+"/split_folder/"+f2
                     os.rename(old_name, new_name)

                 if "all_summary.txt" in f2 and "perm" not in in_dir:
                     old_name=my_in_dir_2+"/split_folder/"+f2
                     new_name = my_in_dir_2+"/split_folder/"+f+"_all_summary.txt"
                     os.rename(old_name, new_name)

    '''
    f1=f.replace(".","_")


    up_f = "/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/qtl_result"

    full_f = os.path.join(up_f, f)
    full_f_new = os.path.join(up_f, f1)



    if os.path.exists(full_f):
        os.rename(full_f, full_f_new)
    '''
