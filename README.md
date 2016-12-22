# INSTALLATION
### to import into python ###
use

	f2py -c -m UDP_modules critfile.f90 UDP_modules.f90

then add the directory to PYTHONPATH or copy the files

     UDP_modules.so
     UDPClust.py

in the working directory

# USAGE
## In your python script / Jupyter notebook
### To call the subroutine:
      import UDPClust as dp
      clustering=dp.cluster_UDP(dim,trj_tot)
##### Input variables:
        dim :: intrinsic dimension of the input data set

        trj_tot :: coordinates of the points in the data set that I want to cluster
                   should be shaped (N.frames)x(N.coords), or be a list of analogue numpy arrays

        dmat :: (Optional) matrix of distances between data points. If not provided by the user will
                           be computed as the euclidean distances between points in trj_tot

        stride :: (default=1) distance matrix calculation and storage is unpractical for N.frames >~ 10^4
                              use a stride>1 in order to perform the clustering only on a subset of the
                              total dataset. The rest of the points will be assigned later

        dump_dmat :: (default=False) set to True to save the distance matrix on disk as udp-dmat.dat

        coring :: (default=True) define the core of each cluster, as the points with 
                  rho/rho_max>exp(-delta)
                  were rho_max is the density in the cluster's center
                  
        sens :: (default=1) sensibility parameter in the clustering algorithm. Increase to merge more clusters

        delta :: (default=1.0) core set definition parameter
        bigdata :: (default=False) set True if you really want to let the program run with >20k points (it's going to be really slow)
        n_jobs :: (default = -1) number of processor to use for cKDTree.query (-1 will use all of them)

##### Optional:
	sens = sensibility parameter for cluster merging (default=1.0)
	delta = parameter for core set definition (default=1.0)
	stride = stride to use for the clustering. The rest of the points will be assigned later to the same cluster of the closest point

##### Note:
The KDTree branch solve the problem of computing the whole distance matrix, so in principle can be used with a very high number of points
Anyway you can use a stride>1

##### OUTPUT VARIABLES
Here the results of clustering are stored

    clustering.frame_cl   # index of the cluster for reach frame
    clustering.cl_idx     # indexes of frames in each cluster
    clustering.n_clusters # number of clusters identified
    clustering.rho        # density for each frame
Other internal variables are

    clustering.filt       # 1 if the computed density is not statistically realiable, 0 if was fine. Points with filt=1 were given the density of the closer filt=0 point
    clustering.cores_idx  # indexes of frames in each cluster's core (i.e. points with high density)


##### Other functions
     clustering.get_core_trajs() # return the discrete trajectories using the coring approach of Hummer and Buchete

     clustering.get_centers()     #Computes the average position of each cluster (Not the best choice for a "center". You should use the argmax(rho) for each cluster

## From terminal




# AUTHORS
This class was written by Giovanni Pinamonti, SISSA, Trieste, 2016

The fortran modules are based on a program written by Alex Rodriguez

Please cite 
d'Errico et al., PNAS, 2016 (soon to be published)