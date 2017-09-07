# This script is used to test the MSM
# construction with coring after clustering
# using PyEmma
#
# Please refer to PyEmma documentation for more information

import pyemma
import pyemma.msm as msm
import numpy as np
from udpclust import UDPClust as dp

tica_traj=[]
for i in range(4):
    fname='DATA/test-its-traj'+str(i)+'.dat'
#    tr=tica_traj[i][::25]
    tica_traj.append(np.loadtxt(fname))

cl_dpa=dp.cluster_UDP(dim=6,trj_tot=tica_traj,stride=10)

ctrajs=cl_dpa.get_core_traj()

its=msm.its(ctrajs,lags=range(1,10,1))

np.savetxt('out-msm_its.dat',its.timescales,fmt="%.6e")

