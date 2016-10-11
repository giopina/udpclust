import sys
import numpy as np
from scipy.spatial import distance
from scipy.special import gamma
#import UDP_modules

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
    
    Rho=np.zeros(Nele)
    Rho_err=np.zeros(Nele)
    filt=np.zeros(Nele)
    minknn=8 # minimum number of neighbours to explore

    id_err=0

    c_idx = lambda i,j: i*Nele + j - i*(i+1)/2 - i - 1

    limit=np.min(maxknn,Nele/2)
    if limit%4 != 0:
        limit=limit+4-limit%4

    # get prefactor for Volume calculation
    prefactor=np.pi**(dimint*0.5)/gamma(dimin*0.5+1)
    # check if it's correct
    #    if dimint%2==0:
    #        k=dimint/2
    #        m=1
    #        for i in range(1,k+1):
    #          m=m*i
    #        prefactor=np.pi**k/(float(m))
    #    else:
    #       k=(dimint-1)/2
    #       m=1
    #       for i in range(1,k+1):
    #          m=m*i
    #       n=m
    #       for i in range(k+1,dimint+1):
    #          n=n*i
    #       prefactor=2.*float(m)*(4.*np.pi)**k/(float(n))

    #    iVols=-np.ones(Nele) #int

    for i in range(Nele):
        Vols=np.ones(Nele)*9.9e99 #float
        for j in range(Nele):
            if i!=j:
                Vols[j]=prefactor*dmat[c_idx(i,j)]**dimint
                iVols=np.argsort(Vols)

        Nlist[i,:]=iVols[:maxknn]

        # get nstar
        viol=False
        k=minknn
        n=1
        for k in range(minknn,limit+1,4):
            #do while (.not.viol)
            rhg=float(k)/Vols[k]
            dL=np.abs(4.*rhg*(Vols[k]-Vols[3*k/4]-Vols[k/4])/float(k))
            if dL > critV[k/4]:
                Nstar[i]=k-4 ! ### ha senso?
                break
        if Nstar[i]<minknn : Nstar(i)=minknn ### puo' succedere..?
        
        ### aggiungo vicini che sono 
        ### a una distanza computazionalmente indistinguibile
        ### (ma ha senso sta roba?)
        kadd=(not (Nstar[i] == limit))
        while kadd:
            if np.abs(dmat[c_idx(i,Nlist[i,Nstar[i]])]-\
                      dmat[c_idx(i,Nlist[i,Nstar[i]+1])])\   
                      <9.99e-99:
                Nstar[i]=Nstar[i]+1
                if Nstar[i] == limit: kadd=False
            else:
                kadd=False

        partit=np.ones(4)*Nstar[i]/4
        partit[:Nstar[i]%4]+=1

        x=(np.cumsum(partit)-partit[0]+0.5*partit)*4./Nstar[i]
#        x[0]=0.5*partit[0]
#        x[1]=partit[0]+0.5*partit[1]
#        x[2]=partit[0]+partit[1]+0.5*partit[2]
#        x[3]=partit[0]+partit[1]+partit[2]+0.5*partit[3]
#        x*=4./Nstar[i]
        rhg=float(Nstar[i])/Vols[iVols[Nstar[i]]] # Rho without fit
        
        # get inv rho of the four quarters
        #            j=Nstar(i)/4
        #            a=dfloat(j)
        for i in range(4):
            rh[i]=( Vols[iVols[ np.sum(partit[:i+1]) ]] - Vols[iVols[ np.sum(partit[:i]) ]] ) / float(partit[i])
        #rh[0]=Vols[iVols[partit[0]]]/float(partit[0])
        #rh[1]=(Vols[iVols[partit[0]+partit[1]]]-Vols[iVols[partit[0]]])/float(partit[1])
        #rh[2]=(Vols[iVols[partit[0]+partit[1]+partit[2]]]-Vols[iVols[partit[0]+partit[1]]])/float(partit[2])
        #rh[3]=(Vols[iVols[np.sum(partit)]]-Vols[iVols[np.sum[partit[:3]])/float(partit[3])

        # make the quadratic fit rhj=1/rho+C*j^2
        xmean=np.mean(x)
        ymean=np.mean(rh)
        b=np.sum((x-xmean)**2)
        a=np.sum((x-xmean)*(rh-ymean))
        a=a/b
        rjfit=ymean-a*xmean
        # Perform jacknife resampling for estimate the error (it includes statistical
        # error and curvature error) 
        for i in range(4):
            x_temp=np.concatenate(x[:i],x[i+1:])
            rh_temp=np.concatenate(rh[:i],rh[i+1:])
            xmean=np.mean(x_temp)
            ymean=np.mean(rh_temp)
            b=np.sum((x_temp-xmean)**2)
            a=np.sum((x_temp-xmean)*(rh_temp-ymean))
            a=a/b
            rjk[i]=ymean-a*xmean


        rjaver=np.mean(rjk)
        Rho[i]=4.*rjfit-3.*rjaver
        Rho_err[i]=0.75*np.sum((rjk-rjaver)**2)
        Rho[i]=1./Rho[i]
        Rho_err[i]=max(Rho[i]/sqrt(float(Nstar[i])),Rho[i]*Rho[i]*sqrt(Rho_err[i]))

        del (Vols,iVols)

        # Filter with neighbours density (Iterative version)
        #    filt(:)=.false.
        viol=True
        niter=0
        while (viol):
            niter=niter+1
            viol=False
            nfilt=0
            for i in range(Nele):
                # compute average density in the neighborhood and standard dev.

### VOY POR AQUI
          if (.not.filt(i)) then
             ! ### come puo' essere Rho minore di zero...!?
             if ((Rho(i).lt.0).or.(Rho_err(i).gt.Rho(i))) then
                filt(i)=.true.
                viol=.true.
             else
                a=0.
                b=0.
                n=0
                ! ### questo prob si puo' fare senza loop ma bisogna pensarci (c'e' da considerare i filt...)
                do j=1,Nstar(i) 
                   if (.not.filt(Nlist(i,j))) then
                      a=a+Rho(Nlist(i,j))
                      b=b+Rho(Nlist(i,j))**2
                      n=n+1
                   endif
                enddo
                if (n.gt.0) then
                   a=a/dfloat(n)              ! average
                   if (a*a.le.b/float(n)) then
                      b=dsqrt(b/dfloat(n)-a*a)    ! std. dev.
                      if (((Rho(i)-a)).gt.sqrt(b*b+Rho_err(i)*Rho_err(i))) then
                         filt(i)=.true.
                         viol=.true.
                      endif
                   endif
                else
                   filt(i)=.true.
                   viol=.true.
                endif
             endif
          else
             nfilt=nfilt+1
          endif
       enddo
    enddo
    return
