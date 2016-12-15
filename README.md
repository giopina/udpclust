### to import into python ###
use

	f2py -c -m UDP_modules critfile.f90 UDP_modules.f90

then add the directory to PYTHONPATH or copy the files

     UDP_modules.so
     UDPClust.py

in the working directory


### To call the subroutine: ###
      import UDPClust as dp
      clustering=dp.cluster_UDP(dim,trj_tot)

### Input variables: ###
    dim = dimensionality of the dataset
    trj_tot = trajectory to perform the clustering (a subset of the total data set, usually)
          should be a numpy array shaped (N.frames)x(N.coords)

Optional:
	sens = sensibility parameter for cluster merging (default=1.0)
	delta = parameter for core set definition (default=1.0)
	stride = stride to use for the clustering. The rest of the points will be assigned later to the same cluster of the closest point

##### Note: #####
the distance matrix calculation and storage is unpractical for N. points >10^4
if that happen it is recommended to use a stride>1


### OUTPUT VARIABLES ###
Here the results of clustering are stored

    clustering.frame_cl   # index of the cluster for reach frame
    clustering.cl_idx     # indexes of frames in each cluster
    clustering.n_clusters # number of clusters identified
    clustering.rho        # density for each frame
Other internal variables are

    clustering.filt       # 1 if the computed density is not statistically realiable, 0 if was fine. Points with filt=1 were given the density of the closer filt=0 point
    clustering.cores_idx  # indexes of frames in each cluster's core (i.e. points with high density)


### Other functions ###
     clustering.get_core_trajs() # return the discrete trajectories using the coring approach of Hummer and Buchete

     clustering.get_centers()     #Computes the average position of each cluster (Not the best choice for a "center". You should use the argmax(rho) for each cluster


### AUTHORS ###
This class was written by Giovanni Pinamonti, SISSA, Trieste, 2016

The fortran modules are based on a program written by Alex Rodriguez

Please cite 
d'Errico et al., PNAS, 2016 (soon to be published)