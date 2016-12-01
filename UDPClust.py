#######################################################################
### This Class performs the unsupervised density peak clustering    ###
### as described in                                                 ###
### D'Errico, Facco, Laio, Rodriguez, pending pubblication, 2016    ###
###                                                                 ###
### The original clustering routines were written by Alex Rodriguez ###
### and adapted for Python compatibility by Giovanni Pinamonti      ###
###                                                                 ###
### This class was written by Giovanni Pinamonti                    ###
### SISSA, Trieste, Italy, 2016                                     ###
#######################################################################

import time
import sys
import numpy as np
from scipy.spatial import distance
import UDP_modules

class cluster_UDP:    
    """Class cluster_UDP
    performs the unsupervised density peak clustering
    as described in   
    D'Errico, Facco, Laio, Rodriguez, pending pubblication, 2016

    internal variables:

      # input parameters
    coring      :: bool    :: True to define core sets
    delta       :: float   :: parameter for core set definition
    sensibility :: float   :: parameter in clusters merging

      # dataset information
    Ntot        :: int     :: total number of data points
    Npoints     :: int     :: number of data points in trj_sub
    ND          :: int     :: number of distances
    dim         :: int     :: intrinsic dimension of data set
    stride      :: int     :: points to skip for clustering
    trj_tot     :: ndarray :: total data set
    trj_sub     :: ndarray :: reduced data set; trj_sub=trj_tot[::stride]
    trj_shape   :: list    :: shape of initial input data (it is a list of tuples or a tuple if the input was a single ndarray)
    

      # clustering in/out variables
    frame_cl_sub    :: ndarray :: index of cluster for each frame in trj_sub
    rho_sub         :: ndarray :: density for each frame in trj_sub
    filt_sub        :: ndarray :: filter value for each frame in trj_sub (either 0 or 1)

      # clustering information
      ### TODO: change this to include information of index of traj+ index of frame.
      ###       Right now will only consider the index of the frame in the concatenated trj_tot,
      ###       which may be unpractical to use if the input data-set is composed by more trajectories.
      ###       This is relevant only for MSM applications
    n_clusters  :: int     :: number of clusters
    cl_idx      :: list    :: frames of trj_tot in each cluster
    frame_cl    :: ndarray :: index of cluster for each frame in trj_tot 
    rho         :: ndarray :: density for each frame in trj_tot
    filt        :: ndarray :: filter value for each frame in trj_tot (either 0 or 1)
    cl_idx      :: list    :: frames of trj_tot in each cluster

    cores_idx   :: list    :: frames (of trj_tot) in the core of every cluster ### TODO: maybe add this for sub/tot


    id_err      :: int     :: error flag
    
    """
    #### this should be the constructor
#    def __init__(self,dmat,dim,trj_tot,stride=1,merge=True):
    def __init__(self,dim,trj_tot,dmat=None,stride=1,dump_dmat=False,coring=True,sens=1.0,delta=1.0):
        """Constructor of the class cluster_UDP:
        input variables

        dim :: intrinsic dimension of the input data set

        trj_tot :: coordinates of the points in the data set that I want to cluster
                   should be shaped (N.frames)x(N.coords), or be a list of such numpy arrays

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
        """

        #### store internal variables
        self.coring=coring
        self.sensibility=sens
        self.dim=dim
        self.delta=delta


        # trj_tot can be a a numpy array shaped (N.frames)x(N.coords)
        #         or a list of numpy arrays
        # usa isinstance(tica_traj,list/np.array) to understand which type it has
        ### TODO: maybe I have to change this because with np.concatenate I'm copying the whole trj_tot
        ###       anyway, it is a small memory consumption compared to the size of dmat and cosidering 
        ###       that I'm storing anyway rho, filt, cl_idx, frame_cl
        if isinstance(trj_tot,list):
            try:
                self.trj_tot=np.concatentate(trj_tot)
            except:
                sys.exit('ERROR: problem in the input data set; shape of the arrays in the trj_tot list is wrong')
            self.trj_shape=[ttt.shape for ttt in trj_tot]
        elif isinstance(trj_tot,np.ndarray):
            self.trj_tot=trj_tot
            self.trj_shape=trj_tot.shape
        else:
            sys.exit('ERROR: problem in the input data set; unrecognized type of variable trj_tot')

#        if stride!=1:
            #print 'stride different from 1 not supported yet'
            #return
        assert stride>0, "ERROR: negative stride not supported!"
        assert isinstance(stride,int), "ERROR: stride must be a positive integer"
        self.stride=stride
        self.trj_sub=trj_tot[::self.stride]
        assert self.trj_sub.shape[0]>1, 'ERROR: stride is too large, the subset contains only one point'

#        assert self.trj_sub.shape[0]<20000, 'WARNING: the size of the distance matrix will be large. Maybe you should decrease the stride'
        if self.trj_sub.shape[0]>20000: print 'WARNING: the size of the distance matrix will be large. Maybe you should decrease the stride'

        ### compute the distance matrix if not provided 
        ###  (in this way it will be deleted at the end of the function. suggested to avoid large memory consumption)
        if dmat==None:
            print 'Computing distances'
            dmat=distance.pdist(self.trj_sub)
        else:
            assert dmat.shape[0]==self.trj_sub.shape[0],"ERROR: trj_tot[::stride] and distance matrix shapes do not match"

        self.ND=len(dmat)
        self.Ntot=self.trj_tot.shape[0]
        self.Npoints=self.trj_sub.shape[0]

        ### dump the distance matrix if needed for dimensionality calculation
        if dump_dmat:
            print 'Writing distance matrix on file','udp-dmat.dat'
            maxmem=1e9 #max memory to use for the string, in bytes
            with open('udp-dmat.dat','w') as fh:
                k=0
                stringa=''
                for i in range(self.Npoints):
                    for j in range(i+1,self.Npoints):
                        d=dmat[k]
                        stringa+='%d %d %f\n' % (i+1,j+1,d)
                        k+=1
                        if len(stringa)*sys.getsizeof('a')>maxmem:
                            fh.write(stringa)
                            stringa=''
                fh.write(stringa)
                stringa=''

        ### perform the clustering
        self.__clustering(dmat)
        ### check for errors
        assert not self.__errorcheck(), 'ERROR: Problem in clustering'

        ### core sets
        if self.coring:
            self.__find_core_sets(R_CORE=np.exp(-self.delta))
        else:
            self.__find_core_sets(R_CORE=1.)


    def __clustering(self,dmat):
        #
        # 1) initialize quantities that will be computed by the fortran subroutine
        #   index of cluster for each frame in trj_sub
        self.frame_cl_sub=np.zeros(self.Npoints,dtype=np.int32)
        #   density for each frame in trj_sub
        self.rho_sub=np.zeros(self.Npoints)
        #   filter value for each frame in trj_sub (either 0 or 1)
        self.filt_sub=np.zeros(self.Npoints,dtype=np.int32)
        #    error flag
        self.id_err=np.array(0,dtype=np.int32)
        #
        # 2) call fortran subroutine
        print 'fortran clustering'
        t0=time.time()
        UDP_modules.dp_clustering.dp_advance\
            (dmat,self.frame_cl_sub,self.rho_sub,self.filt_sub,self.dim,self.id_err,self.sensibility)
        print 'Done!',
        print time.time()-t0,
        print 's; now post-processing'
        # 3) post processing of output
        self.frame_cl_sub-=1 #back to python enumeration in arrays

        ### I do this later with all the points!
        self.cl_idx_sub=[ [] for i in range(np.max(self.frame_cl_sub)+1)] #frames for each cluster
        i=0
        for i_cl in self.frame_cl_sub:
            self.cl_idx_sub[i_cl].append(i)
            i+=1
#        self.cl_idx=[ [] for i in range(np.max(self.frame_cl)+1)] #frames for each cluster
#        i=0
#        for i_cl in self.frame_cl:
#            self.cl_idx[i_cl].append(i)
#            i+=1

        # 4) assign densities of nearest-neighbours to the filtered points
        f1=np.where(self.filt_sub==1)[0]
        f0=np.where(self.filt_sub==0)[0]
        for i in f1:
            dists=distance.cdist(self.trj_sub[f0],np.array([self.trj_sub[i]]))[:,0]
            imin=np.argmin(dists)
            self.rho_sub[i]=self.rho_sub[f0[imin]]


        # 5) assign trj_tot points to the clusters
        #   index of cluster for each frame in trj_tot
        self.frame_cl=np.zeros(self.Ntot,dtype=np.int32)
        #   density for each frame in trj_tot
        self.rho=np.zeros(self.Ntot)
        #   filter value for each frame in trj_tot (either 0 or 1)
        self.filt=np.zeros(self.Ntot,dtype=np.int32)
        #   frames of trj_tot for each cluster
        self.cl_idx=[ [] for i in range(np.max(self.frame_cl_sub)+1)] 

        ### if stride==1 don't recompute distances:
        if self.stride==1:
            print 'assign strided points'
            self.frame_cl=self.frame_cl_sub
        #   density for each frame in trj_tot
            self.rho=self.rho_sub
        #   filter value for each frame in trj_tot (either 0 or 1)
            self.filt=self.filt_sub
            self.cl_idx=[ [] for i in range(np.max(self.frame_cl)+1)] #frames for each cluster
            i=0
            for i_cl in self.frame_cl:
                self.cl_idx[i_cl].append(i)
                i+=1
            self.n_clusters=len(self.cl_idx)
            return
        ###
        
        ### this part may take some time
        t0=time.time()
        for iframe in range(self.Ntot):
            frame=self.trj_tot[iframe] # 
            #                dists=distance.cdist(self.trj_sub,np.array([frame]))[:,0]
            sqdists=distance.cdist(self.trj_sub,np.array([frame]),'sqeuclidean')[:,0] # should be faster
            idx=np.argmin(sqdists)
            self.frame_cl[iframe]=self.frame_cl_sub[idx]
            self.rho[iframe]=self.rho_sub[idx]
            self.filt[iframe]=self.filt_sub[idx]
            icl=self.frame_cl_sub[idx]
            self.cl_idx[icl].append(iframe)
        ###
        self.n_clusters=len(self.cl_idx)
        print time.time()-t0
        print "finished clustering"
        return 
    #END FUNCTION __CLUSTERING

    def __errorcheck(self):
        if self.id_err!=0:
            #TODO: change print into sys.exit( ... )
            if self.id_err==1 :   print "Select case error for distance type" #obsolete?
            elif self.id_err==2 : print "Error opening distance file"#obsolete?
            elif self.id_err==3 : print "Error opening coordinates file"#obsolete?
            elif self.id_err==4 : print "Error on distance file format" #obsolete?
            elif self.id_err==5 : print "Error on coordinates file format" #obsolete?
            elif self.id_err==6 : print "Peridic Boundary Conditions option unrecognized" # obsolete?
            elif self.id_err==7 : print "Dimension calculation option unrecognized" # obsolete?
            elif self.id_err==8 : print "Error opening critical values file" #obsolete?
            elif self.id_err==9 : print "Just one cluster"
            elif self.id_err==10: print "Just one cluster after merging"
            elif self.id_err==11: print "Error in assignation"
            elif self.id_err==12: print "Error in clustering: ig undefined; this may mean that there's a maximum in g_i not identified as a cluster center. Maybe you have identical points."
            elif self.id_err==13: print "Error in distance: There are identical points!"
            else : print "unrecognized error"
            return True
        else: return False

        
    ### CORE SETS IDENTIFICATION
    def __find_core_sets(self,R_CORE=np.exp(-1)):
#        print " Questo e' il cutoff sul rapporto della densita' core con il picco:",R_CORE
        print " Identifying core sets using a cutoff of %s with respect to the peak density" % R_CORE
        self.cores_idx=[  [] for i in range(len(self.cl_idx))]
        k_cl=0
        #for cluster in self.cl_idx_sub:
        #    rhomax=np.max(self.rho_sub[cluster]) #search max in rho_sub instead of rho, to save time
        for cluster in self.cl_idx:
            rhomax=np.max(self.rho[cluster]) #search max in rho_sub instead of rho, to save time
            for ipoint in cluster:
                tmp_rho=self.rho[ipoint]
                if tmp_rho/rhomax>R_CORE:
                    self.cores_idx[k_cl].append(ipoint)
            k_cl+=1


    ### CTRAJS 
    def coring(self,tica_traj=None):
        """Return the core-set discrete trajectory, defined with the "coring" approach (Buchete and Hummer, 2008).
           The idea is to assign each frame which is NOT in a core set, to the last core set visited.

        input:
        tica_traj :: (optional, NOT SUPPORTED YET) trajectory to assign (different from trj_tot)
        """
        #        if self.ctrajs==None: ### Could I really need this sometimes?
        ctrajs=[]
        if isinstance(traj_shape,list):
            idx=0 # this counts the idx in the concatenated list self.trj_tot 
            for itraj in range(len(traj_shape)):
                old_icl=len(self.cores_idx) # fake microstate, where u start all the trj and never enter again
                ct=[]
                for iframe in range(traj_shape[itraj][0]):
                    icl=self.frame_cl[idx]
                    if idx not in self.cores_idx[icl]:
                        icl=old_icl
                    ct.append(icl)
                    old_icl=icl
                    idx+=1
                ctrajs.append(np.array(ct))

        elif isinstance(traj_shape,np.ndarray):
            old_icl=len(self.cores_idx) # fake microstate, where u start all the trj and never enter again
            for iframe in range(traj_shape[itraj][0]):
                icl=self.frame_cl[iframe]
                if iframe not in self.cores_idx[icl]:
                    icl=old_icl
                ctrajs.append(icl)
                old_icl=icl
            ctrajs=np.array(ctrajs) ### TODO: check if it is better to return a single ndarray or a list, for PyEmma compatibility
            
        return ctrajs




    ### CTRAJS 
    ### Here it assigns frames from a list of trajectories to the clusters.
    def assign(self,tica_traj):
        """ Deprecated, you should use stride>1 and then coring()"""
        # this part takes some time
#        if self.ctrajs==None:
        ctrajs=[]
        it=0
        # usa isinstance(tica_traj,list/np.array) to understand which type it has
        for tt in tica_traj:
            old_icl=len(self.cores_idx) # fake microstate, where u start all the trj and never enter again
            ct=[]
            for frame in tt:
                #                dists=distance.cdist(self.trj_tot,np.array([frame]))[:,0]
                sqdists=distance.cdist(self.trj_tot,np.array([frame]),'sqeuclidean')[:,0] # should be faster
                idx=np.argmin(sqdists)
                icl=self.frame_cl_sub[idx]
                
                if self.coring and idx not in self.cores_idx[icl]:
                    icl=old_icl
                ct.append(icl)
                old_icl=icl
            ctrajs.append(np.array(ct))
            it+=1
        return ctrajs

    ### CTRAJS 
    def assign_w_rev(self,tica_traj):
        """ Deprecated, you should use stride>1 and then coring()"""
        # this part takes some time
        ### assign both forward and reverse trajectories
#        if self.ctrajs==None:
        tmp_ctrjs=[]
        it=0
        for tt in tica_traj:
            out_icl=len(self.cores_idx) # fake microstate, when you are not in a core set
            ct=[]
            for frame in tt:
                sqdists=distance.cdist(self.trj_tot,np.array([frame]),'sqeuclidean')[:,0]
                idx=np.argmin(sqdists)
                icl=self.frame_cl_sub[idx]
                if idx not in self.cores_idx[icl]:
                    icl=out_icl
                ct.append(icl)
            tmp_ctrjs.append(np.array(ct))
            it+=1

        ctrajs=[]
        for tmp_ct in tmp_ctrjs:
            ct=[]
            old_icl=out_icl
            for icl in tmp_ct:
                if icl==out_icl:
                    icl=old_icl
                ct.append(icl)
                old_icl=icl
            ctrajs.append(np.array(ct))
            
        ctrajs_rev=[]
        for tmp_ct in tmp_ctrjs:
            ct=[]
            old_icl=out_icl
            for icl in tmp_ct[::-1]:
                if icl==out_icl:
                    icl=old_icl
                ct.append(icl)
                old_icl=icl
            ctrajs_rev.append(np.array(ct))

        return ctrajs,ctrajs_rev


    ### cluster centers in the initial coordinates space
    ### (defined as average of each coord over the points
    ###  in the each core set)
    def get_centers(self):
        """ Deprecated, you should use get_centers_idx"""
#        if self.centers==None:
        Nframes,Ncoords=self.trj_tot.shape
        centers=np.zeros((self.n_clusters,Ncoords))
        for icl in range(self.n_clusters):
            centers[icl]=np.average(self.trj_tot[self.cores_idx[icl]],axis=0)
#        for iframe in range(Nframes):
#            centers[self.frame_cl[iframe]]+=self.trj_tot[iframe]            
        return centers
    ### NEW cluster centers in the initial coordinates space
    ### (defined as points with maximum density)
    def get_centers_idx(self):
        cent_idx=[]
        for cluster in self.cl_idx:
            cent_idx.append(cluster[np.argmax(self.rho[cluster])])
        return np.array(cent_idx)


    def dump_dmat(self,name):
        ### I don't know if this works
        print 'computing distance matrix'
        dmat=distance.pdist(self.trj_tot)
        print 'writing distance matrix on file',name+'_udp-dmat.dat'
        k=0
        stringa=''
        for i in range(self.Npoints):
            for j in range(i+1,self.Npoints):
                d=dmat[k]
                stringa+='%d %d %f\n' % (i+1,j+1,d)
                k+=1
        fh=open(name+'_udp-dmat.dat','w')
        fh.write(stringa)
        fh.close()
        del stringa
    

    def dump_cl(self,name):
        fh=open(name+'_cl_idx.dat','w')
        icl=0
        for clust in self.cl_idx:
            icl+=1
            fh.write('Cluster %d:\n'%icl)
            for idx in clust:
                fh.write('  %d\n'%idx)
        fh.close()
