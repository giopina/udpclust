import sys
import numpy as np
from scipy.spatial import distance
from scipy.special import gamma
#import UDP_modules
from critfile import critV

#from UDP_modules import HPSORT

def locknn(dmat,Nele,dimint,maxknn=496):
    # dmat: distance matrix in triangular ofrm
    # Nele: Number of points (maybe useless)
    # dimint: dimensionality of the data set
    # maxknn: maximum number of neighbours to explore (default=496)
    #
    # outputs will be:
    # Rho(Nele)         : Density
    # Rho_err(Nele)     : Density error
    # filt(Nele)      : pnts with anomalous dens
    # Nlist(Nele,maxknn): Neighbour list within dc. 2nd dim maybe can be reduced but maybe not (before was =limit)
    # Nstar(Nele)       : N. of NN taken for comp dens
    #
    #use critfile !!!

    fh=open('nstar.tmp','w')
    
    Rho=np.zeros(Nele)
    Rho_err=np.zeros(Nele)
    Nstar=np.zeros(Nele,dtype=np.int32)
    Nlist=np.zeros((Nele,maxknn),dtype=np.int32)
    minknn=8 # minimum number of neighbours to explore

    id_err=0

    c_idx = lambda i,j: i*Nele + j - i*(i+1)/2 - i - 1

    limit=np.min([maxknn,Nele/2])
    if limit%4 != 0:
        limit=limit+4-limit%4

    # get prefactor for Volume calculation
    prefactor=np.pi**(dimint*0.5)/gamma(dimint*0.5+1)
    # looks ok up to the 5th digit
 
    import time
    t1=0
    t2=0
    t3=0
   
    for i in range(Nele):
#        if i%(Nele/50)==0: print i
        Vols=np.ones(Nele)*9.9e99 #float
        t_start=time.time()
        ### THE SLOW PART IS ACTUALLY HERE ###
        for j in range(Nele):
            if i!=j:
                #Vols[j]=prefactor*dmat[c_idx(i,j)]**dimint
                Vols[j]=dmat[c_idx(i,j)]
        #Vols=distance.squareform(dmat)[:,i] ### this looks even slower
        ######################################
        t1+=time.time()-t_start
        iVols=np.argsort(Vols)
        #print 'primo loop'
        Nlist[i,:]=iVols[:maxknn]
        #t1+=time.time()-t_start

        t_start=time.time()
        # get nstar
        viol=False
        k=minknn
        n=1
        ordvols=prefactor*Vols[iVols[:limit]]**dimint

        Nstar[i]=limit
        for k in range(minknn,limit+1,4):
            rhg=float(k)/ordvols[k-1]
            dL=np.abs(4.*rhg*(ordvols[k-1]-ordvols[3*k/4-1]\
                              -ordvols[k/4-1])/float(k))
            if dL > critV[(k-minknn)/4]:
                Nstar[i]=k-4 ### ha senso?
                break
                
        if Nstar[i]<minknn : Nstar[i]=minknn ### puo' succedere..?
        if Nstar[i]<1: print '1st',Nstar[i]
        ### aggiungo vicini che sono 
        ### a una distanza computazionalmente indistinguibile
        ### (ma ha senso sta roba?)
        kadd=(not (Nstar[i] == limit))
        while kadd:
            if np.abs(dmat[c_idx(i,Nlist[i,Nstar[i]-1])]-\
                      dmat[c_idx(i,Nlist[i,Nstar[i]])])\
                      <9.99e-99:
                Nstar[i]=Nstar[i]+1
                if Nstar[i] == limit: kadd=False
            else:
                kadd=False
        if Nstar[i]==0: print '2nd',Nstar[i]
        partit=np.ones(4,dtype=np.int32)*int(Nstar[i]/4)
        partit[:Nstar[i]%4]+=1
        t2+=time.time()-t_start
        fh.write('%d %d\n'%(i, Nstar[i]))

        x=(np.cumsum(partit)-partit[0]+0.5*partit)*4./Nstar[i]
        rhg=float(Nstar[i])/ordvols[Nstar[i]-1] # Rho without fit
        
        t_start=time.time()
        # get inv rho of the four quarters
        rh=np.zeros(4)
        for i4 in range(4):
            rh[i4]=( ordvols[np.sum(partit[:i4+1])-1] - ordvols[np.sum(partit[:i4])-1] ) / float(partit[i4])
        # make the quadratic fit rhj=1/rho+C*j^2
        xmean=np.mean(x)
        ymean=np.mean(rh)
        b=np.sum((x-xmean)**2)
        a=np.sum((x-xmean)*(rh-ymean))
        a=a/b
        rjfit=ymean-a*xmean
        # Perform jacknife resampling for estimate the error (it includes statistical
        # error and curvature error) 
        rjk=np.zeros(4)
        for i4 in range(4):
            x_temp=np.concatenate((x[:i4],x[i4+1:]))
            rh_temp=np.concatenate((rh[:i4],rh[i4+1:]))
            xmean=np.mean(x_temp)
            ymean=np.mean(rh_temp)
            b=np.sum((x_temp-xmean)**2)
            a=np.sum((x_temp-xmean)*(rh_temp-ymean))
            a=a/b
            rjk[i4]=ymean-a*xmean


        rjaver=np.mean(rjk)
        Rho[i]=4.*rjfit-3.*rjaver
        Rho_err[i]=0.75*np.sum((rjk-rjaver)**2)
        Rho[i]=1./Rho[i]
        if Nstar[i]==0: print '3rd',Nstar[i]
        Rho_err[i]=max(Rho[i]/np.sqrt(float(Nstar[i])),Rho[i]**2*np.sqrt(Rho_err[i]))

        del (Vols,iVols)
    t3+=time.time()-t_start
    print 
    print t1,t2,t3, 's'
    print
    print 'filtriamo'
    t1=0
    t2=0
    t3=0
    # Filter with neighbours density (Iterative version)
    viol=True
    niter=0
    filt=np.zeros(Nele,dtype=np.int32)
    while (viol):
#        if niter%10==0: print niter
        niter=niter+1
        viol=False
        nfilt=0
   
        for i in range(Nele):
            # compute avg rhox in NN and std.
            ### come puo' essere Rho minore di zero...!?
            if filt[i]==1: continue
            t_start=time.time()
            if Rho[i]<0 or Rho_err[i] > Rho[i]:
                filt[i]=1
                viol=True
                t1+=time.time()-t_start
                continue
            t1+=time.time()-t_start
            t_start=time.time()
            #if len(np.intersect1d(Nlist[i,:],no_filt))==0:
            unfilt_nn=np.where(filt[Nlist[i,:]]<1)[0]
                #np.delete(no_filt,i)
            if len(unfilt_nn)==0:
                filt[i]=1
                viol=True
                t2+=time.time()-t_start
                continue
            t2+=time.time()-t_start
            t_start=time.time()
            a=np.mean(Rho[unfilt_nn])
            b=np.std(Rho[unfilt_nn])
            if (Rho[i]-a)**2 > b**2+Rho_err[i]**2:
                filt[i]=1
                viol=True
            t3+=time.time()-t_start
    
    print 'TOT T',t1,t2,t3,'s'

    fh.close()

    return  Rho, Rho_err, filt, Nlist, Nstar
