def clustering(dmat,cutoff):
    import numpy as np
    import math
    NFRAMES=dmat.shape[0]
    if NFRAMES!=dmat.shape[1]:
        print('#error: C is not a square matrix')
        
    print("# NFRAMES:", NFRAMES)

    frames=[ i for i in range(0,NFRAMES)]

    rho=np.sum(dmat<cutoff,axis=0)
    d_max=np.max(dmat)
    delta=[]
    for i in range(NFRAMES):
        temp_min=d_max
        higher_r_idx=np.where(rho>rho[i])[0]
        if higher_r_idx.shape[0]==0:
            delta.append(temp_min)
            continue
        temp_min=np.min(dmat[i,higher_r_idx])
        delta.append(temp_min)

    return rho, delta
        
               
