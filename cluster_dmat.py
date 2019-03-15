### ---------------------------------------------
### Example script that uses the UDPClust class 
### to perform unsupervised density peak clustering
### given a distance matrix file
### ---------------------------------------------
### usage: python ???
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

dmat=np.sort(dist,axis=1)
Nlist=np.argsort(dist,axis=1)

maxknn=100
dp.cluster_UDP(dim,traj,dmat=dmat[:,:maxknn],Nlist=Nlist[:,:maxknn],maxknn=maxknn,sens=1.0)


