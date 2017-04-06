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
import UDP_modules

import UDPClust as dp
fname=sys.argv[1]
dim=int(sys.argv[2])
traj=[]
for line in open(fname,'r'):
    traj.append([float(x) for x in line.split()])
traj=np.array(traj)
#np.random.shuffle(traj)
print('shape of input array =',traj.shape)


maxknn=496
from scipy.spatial import distance
from scipy.spatial import cKDTree

tree=cKDTree(traj)
dmat,Nlist=tree.query(traj,k=maxknn+1,n_jobs=-1)

Nlist=Nlist[:,1:]
dmat=dmat[:,1:]
Nlist+=1

#print dmat[0]
#print dmat.shape
#order=np.argsort(dmat,axis=1)
#dmat=dmat[order]
#Nlist=Nlist[order]
print('dmat computed')
Npoints=traj.shape[0]
# 1) initialize quantities that will be computed by the fortran subroutine
#   density for each frame in trj_sub

rho=np.zeros(Npoints)
rho_err=np.zeros(Npoints)
#   filter value for each frame in trj_sub (either 0 or 1)
filt=np.zeros(Npoints,dtype=np.int32)
#   error flag
id_err=np.array(0,dtype=np.int32)
#   dimension
dim=np.array(dim,dtype=np.int32)
#   Neighbour list within dc
#Nlist=np.ones((Npoints,maxknn),dtype=np.int32,order='F')*-1
Nlist=np.array(Nlist,dtype=np.int32,order='F')
#   N. of NN taken for comp dens
Nstar=np.zeros(Npoints,dtype=np.int32)
#

# 2) call fortran subroutine
import time
print('fortran locknn')
t0=time.time()
UDP_modules.dp_clustering.get_densities\
    (id_err,dmat,dim,rho,rho_err,filt,Nlist,Nstar)
print('Done!'),
print(time.time()-t0),

fh=open(sys.argv[3]+'_rho.dat','w')
iframe=0
for r in rho:
    fh.write('%d %.6e %.6e\n'%(iframe,r,filt[iframe]))
    iframe+=1
fh.close()
