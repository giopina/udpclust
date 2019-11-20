# INSTALLATION
### to import into python ###
make
(or use the command: f2py -c -m UDP_modules critfile.f90 UDP_modules.f90)

then add the directory to PYTHONPATH 
(I know, there's surely a better way to do this installation, but I don't know how to do it)

# USAGE
## In your python script / Jupyter notebook
### To call the subroutine:
      from udpclust import UDPClust as dp
      clustering=dp.cluster_UDP(dim,trj_tot)

##### Input variables:
      dim :: intrinsic dimension of the input data set (if unknown can be estimated using the algorithm described in [5])

      trj_tot :: coordinates of the points in the data set that I want to cluster should be shaped (N.frames)x(N.coords), or be a list of such numpy arrays

      dmat :: (Optional) matrix of distances between data points. If not provided by the user will be computed as the euclidean distances between points in trj_tot (I never use it, so this is not 100% tested)

      stride :: (default=1) even with KDTrees distance matrix calculation and storage is unpractical for N.frames >~ 10^5 (use a stride>1 in order to perform the clustering only on a subset of the total dataset. The rest of the points will be assigned later)

      dump_dmat :: (default=False) set to True to save the distance matrix on disk as udp-dmat.dat

      coring :: (default=True) define the core of each cluster, as the points with rho/rho_max>exp(-delta) were rho_max is the density in the cluster's center
                  
      sens :: (default=1) sensibility parameter in the clustering algorithm. Increase to merge more clusters

      delta :: (default=1.0) core set definition parameter
      
      bigdata :: (default=False) set True if you really want to let the program run with >100k points (it's going to be slow and use a lot of memory)
      
      n_jobs :: (default = -1) number of processor to use for cKDTree.query (-1 will use all of them)
      
      i_noise :: (default = 0.0001) random gaussian noise that will be added to your data points prior to the clustering if identical points are found. A warning will be printed. Be careful with it!

      maxknn :: (default = 496) maximum number of nearest neighbors to explore

##### Note:
This branch uses the KDTree implementation which solves the problem of computing the whole distance matrix, so in principle it can be used with a very high number of points (check the input parameter bigdata).
Anyway you can use a stride>1 if you need/want it

##### OUTPUT VARIABLES
Here the results of clustering are stored

    clustering.frame_cl    # index of the cluster for reach frame
    clustering.cl_idx      # indexes of frames in each cluster
    clustering.n_clusters  # number of clusters identified
    clustering.rho         # density for each frame
    clustering.centers_rho # density of the center of each cluster
    clustering.centers_idx # indexes of the center of each cluster  (i.e. points with high density). If you provided multiple input trajectories the idx will refer to a "concatenated trajectory". This can probably be fixed.
    clustering.clustercenters # coordinates of the centers of each cluster

Other internal variables are

    clustering.filt       # 1 if the computed density is not statistically realiable, 0 if it is ok. Points with filt=1 are given the density of the closer filt=0 point

##### Other functions
     clustering.get_core_trajs() # return the discrete trajectories using the coring approach of Hummer and Buchete [1]

     clustering.get_centers()     #Computes the average position of each cluster (Not the best choice for a "center". You should use the argmax(rho) for each cluster


# AUTHORS
This class was written by Giovanni Pinamonti, SISSA, Trieste, 2016

The fortran modules are based on a program written by Alex Rodriguez.

Please cite d'Errico et al. [3] if you use this.

# References

###### Core-MSM approach:
[1] Buchete, Nicolae-Viorel, and Gerhard Hummer. "Coarse master equations for peptide folding dynamics." *The Journal of Physical Chemistry B* 112.19 (2008): 6057-6069.

###### Density peak clustering:
[2] Rodriguez, Alex, and Alessandro Laio. "Clustering by fast search and find of density peaks." *Science* 344.6191 (2014): 1492-1496.

[3] d'Errico, Maria, et al. "Automatic topography of high-dimensional data sets by non-parametric Density Peak clustering." *arXiv preprint arXiv:1802.10549* (2018).

###### Others:
[4] Rodriguez, Alex, et al. "Computing the Free Energy without Collective Variables." *Journal of chemical theory and computation* 14.3 (2018): 1206-1215.

[5] Facco, Elena, et al. "Estimating the intrinsic dimension of datasets by a minimal neighborhood information." *Scientific reports* 7.1 (2017): 12140.
