### ---------------------------------------------
### Example script that uses the UDPClust class 
### to perform unsupervised density peak clustering
### given a distance matrix file
### ---------------------------------------------
### usage: python cluster_dmat.py trajfile distancefile dim sensibility outputname
### 
### NB: the dimensionality is a critical parameter and sometimes not trivial to estimate
###     ask to Elena Facco, for further help on this
### ---------------------------------------------
### Written by Giovanni Pinamonti, FU Berlin, 2019
### ---------------------------------------------


import numpy as np
import sys
from udpclust import UDP_modules
from udpclust import UDPClust as dp
fname=sys.argv[1]
traj=np.loadtxt(fname)

fname=sys.argv[2]
dist=np.loadtxt(fname)

dim=int(sys.argv[3])
sens=float(sys.argv[4])

maxknn=500
dmat=np.sort(dist,axis=1)[:,1:maxknn+1]
Nlist=np.argsort(dist,axis=1)[:,1:maxknn+1]


cl_dp=dp.cluster_UDP(dim,traj,dmat=dmat,Nlist=Nlist,maxknn=maxknn,sens=sens)
print(cl_dp.n_clusters)

cl_dp.dump_cl(sys.argv[5])
cl_dp.dump_frames(sys.argv[5])
