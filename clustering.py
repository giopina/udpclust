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
