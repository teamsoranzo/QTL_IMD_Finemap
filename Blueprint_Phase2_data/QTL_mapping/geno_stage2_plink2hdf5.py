import sys
sys.path.append('/nfs/team151_data03/WP10_release/release_V1_Oct2016/LIMIX/pysrc/')
from CFG.settings import *
import h5py
import pdb
import os
import limix
from limix.io.plink import readPED
from include.utils import smartDumpDictHdf5

if __name__=='__main__':

    chrom_i = int(sys.argv[1])

    if 'debug' in sys.argv:     pdb.set_trace()

    in_dir = "/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/WGS_plink"
    out_dir = "/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/Limix/geno_file/WGS_hdf5"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # importing bed file
    bfn = os.path.join(in_dir, 'chrom%d' % chrom_i)
    RV = readPED(bfn, standardize = False)

    out_file = os.path.join(out_dir, 'chrom%d.h5' % chrom_i)
    fout = h5py.File(out_file, 'w')
    smartDumpDictHdf5(RV, fout)
    fout.close()

