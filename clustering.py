### ---------------------------------------------
### Example script that uses the UDPClust class 
### to perform unsupervised density peak clustering
### ---------------------------------------------
### usage: python clustering.py input_file_name dimensionality output_file_name
### 
### NB: the dimensionality is a critical parameter and sometimes not trivial to estimate
###     ask to Elena Facco, for further help on this
### ---------------------------------------------
### Written by Giovanni Pinamonti, SISSA, Trieste, 2016
### ---------------------------------------------

###TODO: fai gli import con try:
#import numpy as np
import sys

### copied from baRNAba
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
    sys.stderr.write('# Python 2.7 is required. Aborting \n')
    sys.exit(1)
else:
    import argparse

import UDPClust as dp

def parse():
    parser=argparse.ArgumentParser(description='Unsupervised Density Peak Clustering\nplease cite dErrico et al., 2016')
    parser.add_argument("-f", dest="fname",help="input file name",default=None,required=True)

    parser.add_argument("--dim", dest="dim",help="intrinsic dimension of data set",required=True,default=None,type=int)

    parser.add_argument("--stride", dest="stride",help="stride to use to perform clustering. The rest of the data will be assigned afterwards",required=False,default=1,type=int)

    parser.add_argument("--delta", dest="delta",help="core set definition parameter",required=False,default=1.0,type=float)

    parser.add_argument("--sens", dest="sens",help="sensibility parameter in the clustering algorithm. Increase it to merge more clusters.",required=False,default=1.0,type=float)

    parser.add_argument("--coring", dest="coring",help="identify core sets",action="store_true",required=False)

    parser.add_argument("--bigdata", dest="bigdata",help="set to true to work with > 20k points (will be slooow)",action="store_true",required=False)

    parser.add_argument("-o",dest="oname",help="output file name", required=False,default=None)

    args=parser.parse_args()

    return args

def main():
    ### copiato da baRNAba
    # check at startup!
    try:
        import numpy as np
    except ImportError:
        sys.stderr.write('# Numpy is not installed \n')
        sys.exit(1)
    ###
    
    args=parse()
    if args.oname==None:
        args.oname=args.fname

    traj=[]
    for line in open(args.fname,'r'):
        traj.append([float(x) for x in line.split()])
    traj=np.array(traj)
    print ('shape of input array =',traj.shape)
    cl=dp.cluster_UDP(args.dim,traj,stride=args.stride,coring=args.coring,delta=args.delta,sens=args.sens,bigdata=args.bigdata,n_jobs=-1)
    print ('Clustering done')
    #fout=open(sys.argv[3],'wb')
    #pickle.dump(cl,fout,-1)
    #fout.close()
    cl.dump_cl(args.oname)
    cl.dump_frames(args.oname)

####################### MAIN ########################

if __name__ == "__main__":
    main()

