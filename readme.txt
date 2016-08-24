
###to import into python
use
f2py -c -m DPA_clustering critfile.f90 DPA_modules.f90
then copy the file
DPA_clustering.so
in the directory where you are working

###To call the subroutine:
DPA_clustering.dp_clustering.dp_advance(dmat,frame_cl,rho,filt,dim)

### input variables:
dmat=np.array with distances, triangular form, like the output of scipy.distance.pdist()
dim=dimensionality of the data set

### input/output variables
Here the results of clustering are stored
frame_cl=np.zeros(Npoints,dtype=np.int32) # index of the cluster for reach frame
rho=np.zeros(self.Npoints) # density for each frame
filt=np.zeros(self.Npoints,dtype=np.int32) # 1 if the density is statistically unrealiable, 0 if it's fine
