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



    fname = in_dir.split("/")[-2]+'_summary'
    if opt.peer:    fname += '_peer'
    if opt.perm:    fname += '_perm'
    fname+= '.txt'
    out_file = os.path.join(in_dir, fname)

    print "output file is: "+out_file
    if not os.path.exists(out_file) or opt.recalc:

        table = {}

        runs_folder = 'runs'
        #if opt.peer:    runs_folder += '_peer'
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
                try:
                    temp = {}

                    gene_info = data.getGeneInfo(geneID)
                    temp.update(gene_info)
                    # single trait
                    pv = fgene['st']['pv'][0, :]
                    nsnps = pv.shape[0]
                    qv = fgene['st']['qv'][0, :]
                    pv_bonf = nsnps * pv
                    beta = fgene['st']['beta'][0, :]
                    pos = fgene['snp_info']['pos'][:]
                    rs = fgene['snp_info']['rs'][:]


                    temp['pv'] = ['{:.3e}'.format(float(elem)) for elem in pv ]

                    #new_pv =temp['pv']
                    #new_pv[:]=[x*nsnps for x in new_pv]

                    temp['beta'] = ['{:.3e}'.format(float(elem)) for elem in beta ]

                    temp['pv_bonf'] = pv*nsnps

                    temp['pv_bonf'][temp['pv_bonf'] > 1] = 1


                    temp['pv_storey'] = ['{:.3e}'.format(float(elem)) for elem in qv ]

                    temp['pos'] = pos

                    temp['rs'] = rs
                    temp['geneID'] = [str(geneID)]*len(pv)
                    temp['file']   = [str(file)]*len(pv)

                    DAT =  np.column_stack((temp['geneID'],temp['rs'],temp['pos'],temp['pv'],temp['beta'],temp['pv_bonf'],temp['pv_storey']))
                    #DAT =  np.column_stack((temp['pv_bonf'],temp['pv_storey']))


                    ftxt=open(in_dir + "/all_summary.txt",'a')


                    np.savetxt(ftxt,DAT, delimiter=" ", fmt="%s")
                    ftxt.close()
                except:
                    print geneID, 'failed'
                    continue

                #append the temp table into the big table
                #for key in temp.keys():
                #    smartAppend(table, key, temp[key])
            f.close()
#        for key in table.keys():
#            table[key] = sp.array(table[key])

#        print '.. correct for multiple testing'
        #table['pv_bonf'][table['pv_bonf'] > 1] = 1.
        #table['qv_all'] = FDR.qvalues(table['pv_bonf'])
        #print 'no eQTLs at FDR 0.10:', (table['qv_all']<0.10).sum()
        #print 'no genes:', table['qv_all'].shape[0]

#        fout = h5py.File(out_file, 'w')
#        smartDumpDictHdf5(table,fout)
#        fout.close()

#    else:

#        f = h5py.File(out_file, 'r')
#        R = {}
#        for key in f.keys():
#            R[key] = f[key][:]
#        f.close()



