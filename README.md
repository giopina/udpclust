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


##### TO DO: #####
1) Add the option to feed directly a distance matrix
2) Add stride and automatic assignment


##### Note: #####
the distance matrix calculation and storage is unpractical for N. points >10^4
if that happen it is recommended to use a subset of the total data set (trj_tot)
to perform the clustering and subsequently use the method assign(trajs) to assign
all the dataset to the clusters


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
clustering.assign(trajs):
	Assigns the frames from a list of trajectories to the clusters identified before
clustering.get_centers():
	Computes the average position of each cluster


### AUTHORS ###
This class was written by Giovanni Pinamonti, SISSA, Trieste, 2016
The fortran modules are based on a program written by Alex Rodriguez
Please cite 
d'Errico et al., PNAS, 2016 (soon to be published)
