### ---------------------------------------------
### Example program that uses the UDPClust class 
### to perform unsupervised density peak clustering
### ---------------------------------------------
### usage: python clustering.py input_file_name dimensionality output_file_name
### 
### NB: the dimensionality is a critical parameter and sometimes not trivial to estimate
###     ask to Elena Facco, for further help on this
### ---------------------------------------------
### Written by Giovanni Pinamonti, SISSA, Trieste, 2016
### ---------------------------------------------

import numpy as np
import sys

#import cPickle as pickle
import UDPClust as dp
fname=sys.argv[1]
dim=int(sys.argv[2])
traj=[]
for line in open(fname,'r'):
    traj.append([float(x) for x in line.split()])
traj=np.array(traj)
print 'shape of input array =',traj.shape
cl=dp.cluster_UDP(dim,traj)
print 'Clustering done'
#fout=open(sys.argv[3],'wb')
#pickle.dump(cl,fout,-1)
#fout.close()
cl.dump_cl(sys.argv[3])
