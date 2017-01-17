### This is old and useless. Probably it doesn't work

### ---------------------------------------------
### Example script that uses the UDPClust class 
### to perform unsupervised density peak clustering
### Use this script if your trajectory is too big to
### be clustered directly.
### ---------------------------------------------
### usage: suppose you have your data set stored in big_input_file
### 1) save a representative subset of your data set in small_input_file
### 2) run 
###        python clustering.py small_input_file big_input_file dimensionality output_name
### ---------------------------------------------
### output: - output_name_cl_idx.dat has the indexes of the frames in small_input_file
###         - output_name_big_cl_idx.dat has the indexes of the frames in big_input_file
### ---------------------------------------------
### NB: the dimensionality is a critical parameter and sometimes not trivial to estimate
###     ask to Elena Facco, for further help on this
### ---------------------------------------------
### Written by Giovanni Pinamonti, SISSA, Trieste, 2016
### ---------------------------------------------

import numpy as np
import sys

#import cPickle as pickle
import UDPClust as dp

f_small=sys.argv[1]
f_big=sys.argv[2]
dim=int(sys.argv[3])
traj_small=[]
for line in open(f_small,'r'):
    traj_small.append([float(x) for x in line.split()])
traj_small=np.array(traj_small)
print 'shape of input array for clustering =',traj_small.shape
cl=dp.cluster_UDP(dim,traj_small,coring=False)
print 'Clustering done'
cl.dump_cl(sys.argv[4])

traj_big=[]
for line in open(f_big,'r'):
    traj_big.append([float(x) for x in line.split()])
traj_big=[np.array(traj_big)]

print 'now assigning'
frames_idx=cl.assign(traj_big)
clusters_big=[[] for i in range(cl.n_clusters)]
iframe=0
for fr in frames_idx[0]:
    clusters_big[fr].append(iframe)
    iframe+=1

fh=open(sys.argv[4]+'_big_cl_idx.dat','w')
icl=0
for clust in clusters_big:
    icl+=1
    fh.write('Cluster %d:\n'%icl)
    for idx in clust:
        fh.write('  %d\n'%idx)
fh.close()



    

