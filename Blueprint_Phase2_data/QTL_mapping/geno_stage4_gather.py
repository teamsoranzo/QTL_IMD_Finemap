import sys
sys.path.append('/nfs/team151_data03/WP10_release/release_V1_Oct2016/LIMIX/pysrc/')
from CFG.settings import *
import h5py
import pdb
import os
import scipy as sp
import limix
from include.utils import smartDumpDictHdf5

if __name__=='__main__':

    in_dir =sys.argv[1]
    #in_dir = "/lustre/scratch114/projects/hematopoiesis/Blueprint/limix/data_sep24/WGS/filter/WGS_hdf5_filt/neut.geno"
    out_file = 'all_chrom.h5'
    out_file = os.path.join(in_dir, out_file)
    fout = h5py.File(out_file, 'w')
    genoGroup = fout.create_group('genotypes')

    RV = {}

    n_chroms = 22
    for chrom_i in range(1, n_chroms+1):

        print '.. copying chromosome %d'%chrom_i

        # load genotype
        in_file = os.path.join(in_dir, 'chrom%d.h5' % chrom_i)
        fin = h5py.File(in_file, 'r')
        chromGroup = genoGroup.create_group('chrom%d' % chrom_i)
        chromGroup.create_dataset('matrix',data=fin['genotypes']['matrix'][:].T)
        chromGroup.create_dataset('RRM',data=fin['RRM'][:])
        h5py.h5o.copy(fin['genotypes'].id, 'col_headers', chromGroup.id, 'col_headers')
        h5py.h5o.copy(fin['genotypes'].id, 'row_headers', chromGroup.id, 'row_headers')

        # summing up Kpop
        if 'RRM' not in RV.keys():
            RV['RRM']  = fin['RRM'][:]
        else:
            RV['RRM'] += fin['RRM'][:]

        fin.close()

    smartDumpDictHdf5(RV, genoGroup)
    fout.close()

