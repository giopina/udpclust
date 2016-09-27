import pyemma
import numpy as np
import pyemma.coordinates as coor
import sys
import pyemma.msm as msm
import my_tools as my
import readwrite as rw
import cPickle as pickle
from scipy.spatial import distance
import gvec_func as gv
import dihedrals_func as dih
import ClustDep2 as dp

directo='/u/sbp/giopina/srnas2/markov/dino/'
trajfiles=[]
for i in ['1','2']:
    for j in ['a','b']:
        trajfiles.append(directo+'run'+i+'/run'+i+j+'/nopbc_tot.xtc')
topfile=directo+'setup.gro'
