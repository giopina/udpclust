import sys
import numpy as np
from scipy.spatial import distance
import UDP_modules

class cluster_UDP:    
   
    #### this should be the constructor
#    def __init__(self,dmat,dim,trj_tot,stride=1,merge=True):
    def __init__(self,dim,trj_tot,dmat=None,stride=1,dump_dmat=False,coring=True,sens=1.0,delta=1.0):
        "Constructor of the class cluster_UDP"

        #### store internal variables
        self.coring=coring
        self.sensibility=sens
        ### compute the distance matrix if not provided ( in this way it should be deleted at the end of the function. suggested to avoid large memory consumption)
        if dmat==None:
            dmat=distance.pdist(trj_tot)
        else:
            assert dmat.shape[0]==trj_tot.shape[0],"trj_tot and distance matrix shapes do not match"
        self.trj_tot=trj_tot #trajectory on thich I made the clustering (a subset of the total data set, usually)
                             #should be shaped (N.frames)x(N.coords)
        if stride!=1:
            print 'stride different from 1 not supported yet'
            return

        self.ND=len(dmat)
        self.Npoints=len(trj_tot)
        self.dim=dim
        self.delta=delta

### dump the distance matrix if needed for dimensionality calculation
        if dump_dmat:
            print 'writing distance matrix on file','udp-dmat.dat'
            maxmem=1e9 #max memory to use for the string, in bytes
            fh=open('udp-dmat.dat','w')
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
            fh.close()

        ### perform the clustering
        self.__clustering(dmat)
        ### check for errors
        assert not self.__errorcheck(), 'Problem in clustering'
        ### assign densities of nearest-neighbours to the filtered points
        f1=np.where(self.filt==1)[0]
        f0=np.where(self.filt==0)[0]
        for i in f1:
            dists=distance.cdist(trj_tot[f0],np.array([trj_tot[i]]))[:,0]
            imin=np.argmin(dists)
            self.rho[i]=self.rho[f0[imin]]
        if self.coring:
            ### find core sets
            self.__find_core_sets(R_CORE=np.exp(-self.delta))
        else:
            self.__find_core_sets(R_CORE=1.)

    def __clustering(self,dmat):
        # 1) initialize quantities that will be computed by the fortran subroutine
        #index of cluster for each frame in trj_tot
        self.frame_cl=np.zeros(self.Npoints,dtype=np.int32)
        #density for each frame in trj_tot
        self.rho=np.zeros(self.Npoints)
        #filter value for each frame in trj_tot (either 0 or 1)
        self.filt=np.zeros(self.Npoints,dtype=np.int32)
        # 2) Fortran subroutine for clustering
        self.id_err=np.array(0,dtype=np.int32)
        #fortran subroutine
        UDP_modules.dp_clustering.dp_advance\
            (dmat,self.frame_cl,self.rho,self.filt,self.dim,self.id_err,self.sensibility)

        self.frame_cl-=1 #back to python enumeration in arrays
        self.cl_idx=[ [] for i in range(np.max(self.frame_cl)+1)] #frames for each cluster
        self.n_clusters=len(self.cl_idx)
        i=0
        for i_cl in self.frame_cl:
            self.cl_idx[i_cl].append(i)
            i+=1

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
        for cluster in self.cl_idx:
            rhomax=np.max(self.rho[cluster])
            for ipoint in cluster:
                tmp_rho=self.rho[ipoint]
                if tmp_rho/rhomax>R_CORE:
                    self.cores_idx[k_cl].append(ipoint)
            k_cl+=1

    ### CTRAJS 
    ### Here it assigns frames from a list of trajectories to the clusters.
    ### TODO: if coring =True (default) it will assign using the "coring" approach (Buchete and Hummer, 2008)
    def assign(self,tica_traj):
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
                icl=self.frame_cl[idx]
                
                if self.coring and idx not in self.cores_idx[icl]:
                    icl=old_icl
                ct.append(icl)
                old_icl=icl
            ctrajs.append(np.array(ct))
            it+=1
        return ctrajs

    ### CTRAJS 
    def assign_w_rev(self,tica_traj):
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
                icl=self.frame_cl[idx]
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
