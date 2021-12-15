import sys
sys.path.append('/nfs/team151_data03/WP10_release/release_V1_Oct2016/LIMIX/pysrc/')
#from CFG.settings import *
import h5py
import pdb
import os
import scipy as sp
from include.utils import smartDumpDictHdf5
import optparse;



usage = "usage: expStage1_cv2hdf5 [Options] --expr_file .."
parser = optparse.OptionParser(usage=usage)

parser.add_option("-f", "--expr_file", dest="expression_filename",
                  help="your input expression file path")

parser.add_option("-l", "--loc_file", dest="loc_filename",
                  help="your input gene location file path")

parser.add_option("-o", "--outfile", dest="h5_filename",
                  help="your output h5 file path, please name it with sufix expr.h5 ")

parser.add_option("-t", "--threshold_expr", dest="threshold_expr",
                  help="threshold for filtering expression read ")

parser.add_option("-m", "--threshold_mean", dest="threshold_mean",
                  help="threshold for filtering mean expression read ")


options, remainder = parser.parse_args()

if __name__=='__main__':


    #in_dir = os.path.join(CFG['data']['base'], 'RNA-seq')



    f_expr=options.expression_filename
    f_geneloc=options.loc_filename
    out_file = options.h5_filename

    cut1=float(options.threshold_expr)
    cut2=float(options.threshold_mean)

    print "input expression file is: " + f_expr
    print "input location file is: " + f_geneloc
    print "input out file is: "+ out_file
    print "expression cut off is: "+ str(cut1)
    print "expression mean cut off is: "+ str(cut2)

    #pdb.set_trace()

    # import expression
    f_expr = options.expression_filename
    f_geneloc = options.loc_filename
    M = sp.loadtxt(f_expr, dtype='str', delimiter='\t')
    Mloc = sp.loadtxt(f_geneloc, dtype=object, delimiter='\t')
    assert (Mloc[1:, 0]==M[1:, 0]).all(), 'gene id do not match'

    # filtering genes of the following criteria:
    # 1. on autosomal chromosome
    # 2. expressed in at least half of the samples
    #       (where expressed means more than 5 counts)
    chrom = Mloc[1:, 1]
    Iaut = sp.array([(chrom[i] not in ['M', 'X', 'Y']) for i in range(chrom.shape[0])])
    expr = M[1:][:, 1:].astype(float).T
    Iexpr = ((expr>cut1).mean(0)>cut2)
    Igene = sp.logical_and(Iaut, Iexpr)

    # export data
    RV = {}
    RV['matrix'] = expr[:, Igene]
    RV['col_header'] = {}
    RV['col_header']['gene_ID'] = M[1:, 0][Igene]
    RV['col_header']['gene_chrom'] = chrom[Igene].astype(float)
    RV['col_header']['gene_start'] = Mloc[1:, 2][Igene].astype(float)
    RV['col_header']['gene_end'] = Mloc[1:, 3][Igene].astype(float)
    RV['row_header'] = {}
    RV['row_header']['sample_ID'] = M[0, 1:]

    #pdb.set_trace()


    fout = h5py.File(out_file, 'w')
    smartDumpDictHdf5(RV, fout)
    fout.close()


