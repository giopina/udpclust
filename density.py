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

import numpy as np
import sys

import UDPClust as dp
fname=sys.argv[1]
dim=int(sys.argv[2])
traj=[]
for line in open(fname,'r'):
    traj.append([float(x) for x in line.split()])
traj=np.array(traj)
#np.random.shuffle(traj)
print 'shape of input array =',traj.shape

#cl=dp.cluster_UDP(dim,traj)
#print 'Clustering done'
#rho=cl.rho
#filt=cl.filt
from scipy.spatial import distance
from UDP_functions import locknn
dmat=distance.pdist(traj)
print 'dmat computed'
rho,rho_err,filt,Nlist,Nstar=locknn(dmat,len(traj),dim)

fh=open(sys.argv[3]+'_rho.dat','w')
iframe=0
for r in rho:
    fh.write('%d %f %f\n'%(iframe,r,filt[iframe]))
    iframe+=1
fh.close()
