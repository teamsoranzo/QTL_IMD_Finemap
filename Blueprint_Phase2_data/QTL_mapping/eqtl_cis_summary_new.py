import sys
sys.path.append('/nfs/team151_data03/WP10_release/release_V1_Oct2016/LIMIX/pysrc/')
from CFG.settings import *
from include.utils import smartAppend
from include.utils import smartDumpDictHdf5
import limix.stats.fdr as FDR
from include.data import QtlData

import scipy as sp
import h5py
import os
import pdb
import glob
import cPickle
import numpy as np
from optparse import OptionParser

if __name__=='__main__':

    parser = OptionParser()
    parser.add_option("--peer", action="store_true", dest='peer', default=False)
    parser.add_option("--perm", action="store_true", dest='perm', default=False)
    parser.add_option("--recalc", action="store_true", dest='recalc', default=False)

    parser.add_option("--in_dir", dest='in_dir', type=str, default=False)
    parser.add_option("--fg", dest='fg', type=str, default=False)
    parser.add_option("--fp", dest='fp', type=str, default=False)

    (opt, args) = parser.parse_args()
    opt_dict = vars(opt)


    in_dir =opt.in_dir
    my_fg = opt.fg
    my_fp = opt.fp

    data = QtlData(my_fg,my_fp,1000000)



    fname = 'summary'
    if opt.peer:    fname += '_peer'
    if opt.perm:    fname += '_perm'
    fname+= '.hdf5'
    out_file = os.path.join(in_dir, fname)

    print "output file is: "+out_file
    if not os.path.exists(out_file) or opt.recalc:

        table = {}

        runs_folder = 'runs'
        if opt.peer:    runs_folder += '_peer'
        if opt.perm:    runs_folder += '_perm'
        f_name = os.path.join(in_dir, '*.hdf5')

        print "file name is "+f_name
        files = glob.glob(f_name)
        files = sp.sort(files)

        for file_i, file in enumerate(files):

            #print str(100 * file_i / float( len(files) )) + '%'
            print file
            try:
                f = h5py.File(file,'r')
            except:
                print 'file corrupted'
                continue
            geneIDs = f.keys()
            for gene_idx, geneID in enumerate(geneIDs):

                fgene = f[geneID]
                print "here0"
                try:
                    temp = {}
                    temp_all = {}
                    # gene info
                    temp['geneID'] = str(geneID)
                    temp['file']   = str(file)
                    # gene info
                    gene_info = data.getGeneInfo(geneID)
                    temp.update(gene_info)
                    # single trait
                    pv = fgene['st']['pv'][0, :]
                    nsnps = pv.shape[0]
                    idx = pv.argmin()
                    temp['pv'] = pv[idx]
                    temp['beta'] = fgene['st']['beta'][0, idx]
                    temp['pv_bonf'] = nsnps * temp['pv']
                    temp['pv_storey'] = fgene['st']['qv'][0, idx]
                    temp['pos'] = fgene['snp_info']['pos'][:][idx]
                    temp['rs'] = fgene['snp_info']['rs'][:][idx]

                    #for summary all
                    '''
                    print "here 1"
                    pv_all = fgene['st']['pv'][0, :]
                    nsnps_all = pv.shape[0]
                    qv_all = fgene['st']['qv'][0, :]
                    pv_bonf_all = nsnps * pv
                    beta_all = fgene['st']['beta'][0, :]
                    pos_all = fgene['snp_info']['pos'][:]
                    rs_all = fgene['snp_info']['rs'][:]
                    temp_all['pv'] = pv_all
                    temp_all['beta'] = beta_all
                    temp_all['pv_bonf'] = pv_bonf_all
                    temp_all['pv_storey'] = qv_all
                    temp_all['pos'] = pos_all
                    temp_all['rs'] = rs_all
                    temp_all['geneID'] = [str(geneID)]*len(pv)

                    DAT =  np.column_stack((temp_all['geneID'],temp_all['rs'],temp_all['pos'],temp_all['pv'],temp_all['beta'],temp_all['pv_bonf'],temp_all['pv_storey']))
                    out_all = os.path.join(in_dir + "/all_summary.txt")
                    print "my out file is : "+out_all
                    ftxt=open(out_all,'a')
                    np.savetxt(ftxt,DAT, delimiter=" ", fmt="%s")
                    ftxt.close()
                    '''

                except:
                    print geneID, 'failed'
                    continue
                #append the temp table into the big table
                for key in temp.keys():
                    smartAppend(table, key, temp[key])
            f.close()

        for key in table.keys():
            table[key] = sp.array(table[key])

        print '.. correct for multiple testing'
        table['pv_bonf'][table['pv_bonf'] > 1] = 1.
        table['qv_all'] = FDR.qvalues(table['pv_bonf'])
        print 'no eQTLs at FDR 0.10:', (table['qv_all']<0.10).sum()
        print 'no genes:', table['qv_all'].shape[0]

        fout = h5py.File(out_file, 'w')
        smartDumpDictHdf5(table,fout)
        fout.close()

    else:

        f = h5py.File(out_file, 'r')
        R = {}
        for key in f.keys():
            R[key] = f[key][:]
        f.close()




