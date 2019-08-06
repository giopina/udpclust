##########################
### written by giopina ###
##########################

import sys

import numpy as np
def twonn(traj,stride=1,frac=0.95,plot=True):
    """Function that estimates the intrinsic dimension of a dataset,
    as described in 
    Facco, E., d'Errico, M., Rodriguez, A. and Laio, A. 
    "Estimating the intrinsic dimension of datasets by a minimal neighborhood information." 
    Scientific reports 7.1 (2017): 12140.
    -----
    Input arguments:
    traj :: the input data as a numpy array of shape NxK, where N is the number of points, K the number of coordinates for each point
    stride :: set stride>1 if you want to use a subset of the total points. Recommended for N>~10^5
    frac :: fraction of points to consider in the fit.
    ----
    Output:
    returns the estimated intrinsic dimension as an integer, and plots the graph log( 1-F(nu) ) vs log(nu)"""
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        sys.stderr.write('# Scipy is not installed \n')
        sys.exit(1)

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        sys.stderr.write('# matplotlib is not installed \n')
        sys.exit(1)

    trj_sub=traj[::stride] # select strided subset of points
    tree=cKDTree(trj_sub) # K-dimensional tree to speed up NN calculatioj 
    dst,ndx=tree.query(trj_sub,k=3,n_jobs=-1) # compute distances to 1st 2 NN
    nu=np.sort(dst[:,2]/dst[:,1]) # compute  and sort nu (Skip the 1st nn which is the point itself)
    Fnu=np.arange(0,len(nu))/float(len(nu)) # Fnu is the position in the ordered array of nu, divided by the total number of nu

    # plot the plot
    #plt.figure()
    if plot:
        plt.plot(np.log(nu),-np.log(1.-Fnu),ls='',marker='.',ms=1)
        plt.xlabel(r'$\log(\nu)$')
        plt.ylabel(r'$\log(F(\nu))$')
        # fit
        i_min=0
        i_max=int(len(nu)*frac)
        m,b=np.polyfit(np.log(nu[i_min:i_max]),-np.log(1.-Fnu[i_min:i_max]),1)
        plt.title("Fit parameters: $m=%.1f$, $q=%.3f$"%(m,b))
        x=np.array([np.log(nu[0]),np.log(nu[-1])])
        plt.plot(x,m*x+b)
        #    plt.show()
    return int(round(m))

def main():
    try:
        import numpy as np
    except ImportError:
        sys.stderr.write('# Numpy is not installed \n')
        sys.exit(1)

    #input_file_name='data.dat'
    input_file_name=sys.argv[1]
    data=np.loadtxt(input_file_name)
    stride=1
    print(data.shape)
    plt.figure()
    dim=twonn(data,stride)
    print("The intrinsic dimension estimated is %d"%dim)
    plt.show()
    
if __name__ == "__main__":
    main()
