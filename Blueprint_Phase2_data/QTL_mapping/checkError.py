__author__ = 'yy5'

import os, sys, glob

in_dir = sys.argv[1]


sub_dirs = os.listdir(in_dir)

for my_dir in sub_dirs:

    my_error_array=[];
    log_folder_name = os.path.join(in_dir, my_dir, "log")


    my_log_files = glob.glob(log_folder_name+"/"+"*out*")


    su_count=0
    for my_log_file in my_log_files:
        cmd=''


        if os.path.isfile(my_log_file):
            cmd += 'tail -n 30 '+my_log_file+'| grep Successfully'


            val= os.popen(cmd).read()

            if val:
                su_count += 1
            else:
                my_error_array.append(my_log_file)

    if su_count == len(my_log_files):

        print my_dir+ " has been successfully run"
    else:

        print my_dir+ " has failed"
        print "please check the following files for error"
        print my_error_array
    print ""





