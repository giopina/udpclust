
import numpy as np
import math



def larod(dmat,cutoff):
    NFRAMES=np.shape[0]
    if NFRAMES!=np.shape[1]:
        print '#error: C is not a square matrix'

    print "# NFRAMES:", NFRAMES

    rho=numpy.zeros(NFRAMES)
    frames=[ i for i in range(0,NFRAMES)]
    d_max=0.0
    for i in range(NFRAMES):
        for j in range(NFRAMES):
            if dmat[i][j]<cutoff and i!=j:
                rho[i]+=1
            if d>d_max:
                d_max=d
    delta=[]
    for i in range(NFRAMES):
        temp_min=d_max
        for j in range(NFRAMES):
            if rho[j]>rho[i]:
                if dmat[i][j]<temp_min:
                    temp_min=dmat[i][j]
        delta.append(temp_min)

    return rho, delta
        
               
