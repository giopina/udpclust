//
// Created by marscher on 1/11/17.
//

#include <limits>

#include "udpclust.h"
#include "heap_sort.h"


// mark survivors, based on filter vector
int UDPClust::get_survivors() {
    int Nsurv = 0;
    for (size_t i = 0; i < survivors.size(); ++i) {
        if (!filter[i]) {
            Nsurv += 1;
            survivors[Nsurv] = i;
        }
    }
    return Nsurv;
}

// This function allows to get the distances in a matrix like style
double UDPClust::gDist(size_t i, size_t j) {
    size_t k, l, m;
    l = std::max(i, j);
    m = std::min(i, j);
    k = (m - 1) * Nele - (m * m + m) / 2 + l;
    return dist_mat[k];
}

/*
 *     subroutine clustering(id_err)
      implicit none
      integer :: id_err
      !! Local variables
      integer :: i,j,k
      integer :: ii,jj,kk
      integer :: ig
      integer :: l
      logical :: idmax
      real*8,allocatable :: Rho_prob(:)   ! Probability of having maximum Rho
      real*8,allocatable :: Rho_copy(:)
      integer,allocatable :: iRho(:),ordRho(:)
      integer,allocatable :: eb(:,:)    ! Border elements
      !integer,allocatable :: survivors(:)  ! ###
      real*8 :: d,dmin
      logical :: extend
      id_err=0

      !! Identify centers: delta>dc eqv. Max(rho) within dc
      Cluster(:)=0              !
      Nclus=0
      ! Here I compute the probability of having density rho, g_i
      allocate (Rho_prob(Nele))
      Rho_prob(:)=0.
      call get_survivors() ! ###
      ! ### change this with survivors
      do ii=1,Nsurv
         i=survivors(ii)
         do jj=1,Nsurv
            j=survivors(jj)
            Rho_prob(i)=Rho_prob(i)-log(1.d0+exp(2.*(Rho(j)-Rho(i))/sqrt(Rho_err(i)**2+Rho_err(j)**2)))
         enddo
      enddo

      ! copy of rho (not efficient, but clarifies the code) ### !!!
      allocate (Rho_copy(Nele))
      allocate (iRho(Nele))
      Rho_copy(:)=-Rho_prob(:)
      call HPSORT(Nele,Rho_copy,iRho) ! iRho contains the order in density (iRho(1) is the element with highest Rho...)
      deallocate (Rho_copy)
      allocate (ordRho(Nele))
      do i=1,Nele
         ordRho(iRho(i))=i                 ! ordRho is the complementary of iRho. Given an element, ordRho returns its order in density
      enddo

      ! ###
      ! Now I'm getting the clusters
      ! ### change this with survivors
      do ii=1,Nsurv
         i=survivors(ii)
         idmax=.true.
         j=1
         do while (idmax .and. (j.le.Nstar(i)))
            if ((ordRho(i).gt.ordRho(Nlist(i,j))).and.(.not.filter(Nlist(i,j))))  idmax=.false. ! ### I could probably also change this
            j=j+1
         enddo
         if (idmax) then
            Nclus=Nclus+1
            Cluster(i)=Nclus
         endif
      enddo
      ! ###


      allocate (Centers(Nclus))
      do i=1,Nele
         if (Cluster(i).ne.0) then
            Centers(Cluster(i))=i
         endif
      enddo
      if (Nclus.gt.1) then
         ! Assign not filtered
         ! ### change it with survivors
         do i=1,Nele
            ig=-1
            j=iRho(i)
            if (.not.filter(j).and.Cluster(j).eq.0) then
               dmin=9.9d99
               do k=1,i-1
                  l=iRho(k) ! ### questa cosa mi da davvero un vantaggio?
                  if (.not.filter(l)) then
                     if (gDist(j,l).le.dmin) then
                        ig=l
                        dmin=gDist(j,l) !
                     endif
                  endif
               enddo
               if (ig.eq.-1) then
                  id_err=12
                  RETURN
               else
                  Cluster(j)=Cluster(ig)
               endif
            endif
         enddo


         ! Assign filtered to the same Cluster as its nearest unfiltered neighbour
         ! what happens if all neighbors are filtered
         do i=1,Nele
            ig=-1
            if (Cluster(i).eq.0) then
               dmin=9.9d99
               do j=1,Nstar(i) ! find the min d in not filt elements
                  l=Nlist(i,j)
                  if ((Cluster(l).ne.0).and.(.not.filter(l))) then
                     d=gDist(i,l)
                     if (d.le.dmin) then
                        dmin=d
                        ig=l
                     endif
                  endif
               enddo
               if (dmin.gt.9.8d99) then
                  do j=1,Nele ! find the min d in not filter elements
                     if ((Cluster(j).ne.0).and.(.not.filter(j))) then
                        d=gDist(i,j)
                        if (d.le.dmin) then
                           dmin=d
                           ig=j
                        endif
                     endif
                  enddo
               endif
               if (ig.eq.-1) then
                  id_err=12
                  RETURN
               endif
               Cluster(i)=Cluster(ig)
            endif
         enddo
         ! find border densities
         allocate (Bord(Nclus,Nclus),Bord_err(Nclus,Nclus),eb(Nclus,Nclus))
         Bord(:,:)=-9.9D99
         Bord_err(:,:)=0.
         eb(:,:)=0


         do i=1,Nele ! si puo' fare il loop solo su i filter?
            ig=-1
            if (filter(i)) CYCLE
            ! ### non so farlo bene...
            !do j=1,Nstar(i)
            !   l=Nlist(i,j)
            !   if (cluster(l).ne.cluster(i)) then
            !
            !   endif
            !enddo
            dmin=9.9d99
            do j=1,Nstar(i)
               l=Nlist(i,j)
               if (filter(l)) CYCLE
               if (cluster(l).eq.cluster(i)) CYCLE
               d=gDist(i,l)
               if (d.lt.dmin) then
                  dmin=d
                  ig=l
               endif
            enddo
            if (dmin.gt.9.8d99) CYCLE
            extend=.true.
            if (ig.eq.-1) then
               id_err=12
               RETURN
            endif
            do k=1,Nstar(i)
               if(filter(Nlist(i,k))) CYCLE
               if (cluster(Nlist(i,k)).eq.cluster(i)) then
                  if (gDist(Nlist(i,k),ig).lt.dmin) then
                     extend=.false.
                     EXIT
                  endif
               endif
            enddo
            if (extend) then
               if (Rho_prob(i).gt. Bord(cluster(i),cluster(ig))) then ! this if is useless? no it's not
                  Bord(cluster(i),cluster(ig))=Rho_prob(i)
                  Bord(cluster(ig),cluster(i))=Rho_prob(i)
                  Bord_err(cluster(i),cluster(ig))=Rho_err(i)
                  Bord_err(cluster(ig),cluster(i))=Rho_err(i)
                  eb(cluster(i),cluster(ig))=i
                  eb(cluster(ig),cluster(i))=i
               endif
            endif
         enddo ! i=1,Nele

         ! ### altro cluster (? non capisco sto commento che ho fatto...)
         do i=1,Nclus-1
            do j=i+1,Nclus
               if (eb(i,j).ne.0) then
                  Bord(i,j)=Rho(eb(i,j))
                  Bord(j,i)=Rho(eb(i,j))
               else
                  Bord(i,j)=0.
                  Bord(j,i)=0.
               endif
            enddo
         enddo

         deallocate (Rho_prob)
         deallocate (iRho)
         deallocate(ordRho)
         ! Info per graph pre automatic merging
         allocate (cent(Nclus))
         allocate (cent_err(Nclus))

         do i=1,Nclus
            cent(i)=Rho(Centers(i))
            cent_err(i)=Rho_err(Centers(i))
            ! Modify centers in such a way that get more survival
            do j=1,Nele
               if (.not.filter(j)) then
                  if ((cluster(j).eq.i).and.((Rho(j)-Rho_err(j)).gt.(cent(i)-cent_err(i)))) then
                     cent(i)=Rho(j)
                     cent_err(i)=Rho_err(j)
                  endif
               endif
            enddo
         enddo
      else
         Cluster(:)=1
         id_err=9
      endif
      return
    end subroutine clustering
 */

void UDPClustering::clustering() {
    size_t i, j, k;
    size_t ii, jj, kk;
    size_t ig;
    size_t l;
    bool idmax;
    VecDouble Rho_prob(Nele);  // Probability of having maximum Rho
    VecDouble Rho_copy;
    VecInt iRho, ordRho;
    VecInt Centers;
    //integer, allocatable::eb(:,:)    !Border   elements
    double d, dmin;
    bool extend;

    int Nclus = 0;
    int Nsurv = get_survivors();

    /*
     do ii=1,Nsurv
         i=survivors(ii)
         do jj=1,Nsurv
            j=survivors(jj)
            Rho_prob(i)=Rho_prob(i)-log(1.d0+exp(2.*(Rho(j)-Rho(i))/sqrt(Rho_err(i)**2+Rho_err(j)**2)))
         enddo
      enddo
     */
    // Here I compute the probability of having density rho, g_i
    for (ii = 0; ii < Nsurv; ++ii) {
        i = survivors[ii];
        for (jj = 0; jj < Nsurv; ++jj) {
            j = survivors[jj];
            Rho_prob[i] = Rho_prob[i] - std::log(1.0 + std::exp(2. * (Rho[j] - Rho[i]))
                                                       / std::sqrt(Rho_err[i] * Rho_err[i] + Rho_err[j] * Rho_err[j]));
        }
    }

    /**
      ! copy of rho (not efficient, but clarifies the code) ### !!!
      allocate (Rho_copy(Nele))
      allocate (iRho(Nele))
      Rho_copy(:)=-Rho_prob(:)
      call HPSORT(Nele,Rho_copy,iRho) ! iRho contains the order in density (iRho(1) is the element with highest Rho...)
      deallocate (Rho_copy)
      allocate (ordRho(Nele))
      do i=1,Nele
         ordRho(iRho(i))=i                 ! ordRho is the complementary of iRho. Given an element, ordRho returns its order in density
      enddo
     */
     //TODO: init iRho!
    Rho_copy = Rho; // copy
    heap_sort::sort(Rho_copy.pointer, Rho_copy.size());

    for (size_t i; i < Rho_copy.size();++i) {
        Rho_copy[i] -= Rho_prob[i];
    }
    //  ordRho is the complementary of iRho. Given an element, ordRho returns its order in density
    for (size_t i ; i < Nele; ++i) {
        ordRho[iRho[i]] = i;
    }

    /**
     *   ! ###
      ! Now I'm getting the clusters
      ! ### change this with survivors
      do ii=1,Nsurv
         i=survivors(ii)
         idmax=.true.
         j=1
         do while (idmax .and. (j.le.Nstar(i)))
            if ((ordRho(i).gt.ordRho(Nlist(i,j))).and.(.not.filter(Nlist(i,j))))  idmax=.false. ! ### I could probably also change this
            j=j+1
         enddo
         if (idmax) then
            Nclus=Nclus+1
            Cluster(i)=Nclus
         endif
      enddo
      ! ###
     */
    for (ii =0; ii < Nsurv; ++ii) {
        i = survivors[ii];
        idmax = true;
        j = 1;
        while (idmax && j <= Nstar[i]) {
            if ((ordRho[i] > ordRho[Nlist[i, j]]) && (!filter[Nlist[i, j]])) {
                idmax = False;
                j += 1;
            }
        }
        if (idmax) {
            Nclus += 1;
            cluster[i] = Nclus;
        }
    }

    /*
      allocate (Centers(Nclus))
      do i=1,Nele
         if (Cluster(i).ne.0) then
            Centers(Cluster(i))=i
         endif
      enddo
     */
    for (i=0; i < Nele; ++i) {
        if (Cluster[i] != 0) {
            Centers[Cluster[i]] = i;
        }
    }

    /*
    if (Nclus.gt.1) then
     ! Assign not filtered
     ! ### change it with survivors
     do i=1,Nele
        ig=-1
        j=iRho(i)
        if (.not.filter(j).and.Cluster(j).eq.0) then
           dmin=9.9d99
           do k=1,i-1
              l=iRho(k) ! ### questa cosa mi da davvero un vantaggio?
              if (.not.filter(l)) then
                 if (gDist(j,l).le.dmin) then
                    ig=l
                    dmin=gDist(j,l) !
                 endif
              endif
           enddo
           if (ig.eq.-1) then
              id_err=12
              RETURN
           else
              Cluster(j)=Cluster(ig)
           endif
        endif
     enddo
     */
    /// assign not filtered
    /// TODO: change it with survivors
    if (Nclus <= 1) {
        // TODO: assign to -1;
        Cluster[:] = -1;
        throw ONLY_ONE_CLUSTER;
    }
    for (i =0; i < Nele; ++i) {
        ig = -1;
        j = iRho[i];
        if (!filter[j] && Cluster[j] == 0) {
            dmin = std::numeric_limits<double>::max();
            for (k = 1; k < i - 1; ++k) {
                l = iRho[k];
                if (!filter[l] && gDist(j, l) <= dmin) {
                    ig = l;
                    dmin = gDist(j, l);
                }
            }
            if (ig == -1) {
                throw IG_UNDEFINED;
            } else {
                Cluster[j] = Cluster[ig];
            }
        }
    }

    /*
     ! find border densities
     allocate (Bord(Nclus,Nclus),Bord_err(Nclus,Nclus),eb(Nclus,Nclus))
     Bord(:,:)=-9.9D99
     Bord_err(:,:)=0.
     eb(:,:)=0
     */
    VecDouble2d Bord(Nclus, Nclus));
    VecDouble2d Bord_err(Nclus, Nclus);
    VecInt2d eb(Nclus, Nclus);

    dmin = std::numeric_limits<double>::max();
    /*do i=1,Nele ! si puo' fare il loop solo su i filter?*/
    for (i = 0; i < Nele; ++i) {
        /*
        ig=-1
        if (filter(i)) CYCLE
        dmin=9.9d99
        do j=1,Nstar(i)
           l=Nlist(i,j)
           if (filter(l)) CYCLE
           if (cluster(l).eq.cluster(i)) CYCLE
           d=gDist(i,l)
           if (d.lt.dmin) then
              dmin=d
              ig=l
           endif
        enddo
         */
        ig = -1;
        if (filter[i]) continue;
        for (j = 0; j < Nstar[i]; ++j) {
            l = Nlist[i, j];
            if (filter[l]) continue;
            if (cluster[l] == cluster[i]) continue;

            d = gDist(i, l);
            if (d < dmin) {
                dmin = d;
                ig = l;
            }
        }

        /*
        if (dmin.gt.9.8d99) CYCLE
        extend=.true.
        if (ig.eq.-1) then
           id_err=12
           RETURN
        endif
         */
        if (dmin > std::numeric_limits<double>::max() - 1) continue;
        extend = true;
        if (ig == -1) {
            throw IG_UNDEFINED;
        }

        /*
        do k=1,Nstar(i)
           if(filter(Nlist(i,k))) CYCLE
           if (cluster(Nlist(i,k)).eq.cluster(i)) then
              if (gDist(Nlist(i,k),ig).lt.dmin) then
                 extend=.false.
                 EXIT
              endif
           endif
        enddo
         */
        for (k=0; k < Nstar[i]; ++k) {
            if (filter[Nlist[i, k]]) continue;
            if (cluster[Nlist[i, k]] == cluster[i]) {
                extend = false;
                break;
            }
        }

        /*
         if (extend) then
               if (Rho_prob(i).gt. Bord(cluster(i),cluster(ig))) then ! this if is useless? no it's not
                  Bord(cluster(i),cluster(ig))=Rho_prob(i)
                  Bord(cluster(ig),cluster(i))=Rho_prob(i)
                  Bord_err(cluster(i),cluster(ig))=Rho_err(i)
                  Bord_err(cluster(ig),cluster(i))=Rho_err(i)
                  eb(cluster(i),cluster(ig))=i
                  eb(cluster(ig),cluster(i))=i
               endif
            endif
         */
        if (extend && Rho_prob[i] > Bord[cluster[i], cluster[ig]]) {
            Bord[cluster[i], cluster[ig]] = Rho_prob[i];
            Bord[cluster[ig], cluster[i]] = Rho_prob[i];

            Bord_err[cluster[i], cluster[ig]] = Rho_err[i];
            Bord_err[cluster[ig], cluster[i]] = Rho_err[i];
        }
    } //  enddo ! i=1,Nele

/*
         ! ### altro cluster (? non capisco sto commento che ho fatto...)
         do i=1,Nclus-1
            do j=i+1,Nclus
               if (eb(i,j).ne.0) then
                  Bord(i,j)=Rho(eb(i,j))
                  Bord(j,i)=Rho(eb(i,j))
               else
                  Bord(i,j)=0.
                  Bord(j,i)=0.
               endif
            enddo
         enddo
     */
    for (int i =0; i < Nclus -1; ++i) {
        for (int j = i+1; j < Nclus; ++j) {
            if (eb[i, j] != 0) {
                Bord[i, j] = Bord[j, i] = Rho[eb[i,j]];
            } else {
                Bord[i,j] = Bord[j, i] = 0;
            }
        }
    }

    /*
         deallocate (Rho_prob)
         deallocate (iRho)
         deallocate(ordRho)
         ! Info per graph pre automatic merging
         allocate (cent(Nclus))
         allocate (cent_err(Nclus))
         do i=1,Nclus
            cent(i)=Rho(Centers(i))
            cent_err(i)=Rho_err(Centers(i))
            ! Modify centers in such a way that get more survival
            do j=1,Nele
               if (.not.filter(j)) then
                  if ((cluster(j).eq.i).and.((Rho(j)-Rho_err(j)).gt.(cent(i)-cent_err(i)))) then
                     cent(i)=Rho(j)
                     cent_err(i)=Rho_err(j)
                  endif
               endif
            enddo
         enddo
      else
         Cluster(:)=1
         id_err=9
      endif
      return
    end subroutine clustering
     */

    /// Info per graph pre automatic merging
    cent.reserve(Nclus);
    cent_err.reserve(Nclus);

    for(i=0; i < Nclus; ++i) {
        cent[i] = Rho[Centers[i]];
        cent_err[i] = Rho_err[Centers[i]];
        /// Modify centers in such a way that get more survival
        for (j=0; j < Nele; ++j) {
            if(! filter[j] && (cluster[j] == i) && (Rho[j] - Rho_err[j] > (cent[i] - cent_err[i]) )) {
                cent[i] = Rho[j];
                cent_err[i] = Rho_err[j];
            }
        }
    }

}

void UDPClustering::get_densities() {
    /*
     * *
 *   subroutine get_densities(id_err,dist_mat,Nele,dimint,Rho,Rho_err,filter,Nlist,Nstar)
    use critfile
    implicit none

    integer,parameter :: maxknn=496   ! maximum number of neighbours to explore
    integer,parameter :: minknn=8     ! minimum number of neighbours to explore
    !!Global variables
    real*8,intent(in) :: dist_mat(Nele*(Nele-1)/2)       ! Distance matrix !###
    integer,intent(in) :: Nele                   ! Number of elements
    !    integer,intent(in) :: ND                     ! Number of distances
    integer,intent(in) :: dimint                 ! integer of dimset (avoid real*8 calc)

    ! These variables are used in densities calculation and then passed to clustering
    real*8,intent(inout) :: Rho(Nele)        ! Density
    real*8,intent(inout) :: Rho_err(Nele)    ! Density error
    logical,intent(inout) :: filter(Nele)  ! Pnt with anomalous dens
    integer,intent(inout) :: Nlist(Nele,maxknn) ! Neighbour list within dc ### dimension 2 maybe can be reduced but maybe not (before was =limit)
    integer,intent(inout) :: Nstar(Nele)   ! N. of NN taken for comp dens


    integer :: id_err
    !! Local variables
    integer :: limit
    integer :: i,j,k,m,n
    real*8 :: is,js,ks,ms,ns
    !    integer :: i1 ! ### my integer for loops
    integer :: kadd
    integer :: niter,nfilter
    real*8,allocatable :: Vols(:)
    integer,allocatable :: iVols(:)
    real*8, parameter :: pi=3.14159265359
    real*8 :: prefactor
    real*8 :: rhg,dL,rjaver,rjfit
    !    real*8,dimension(4) :: rh
    !    real*8 :: rjaver,rjfit
    real*8,allocatable :: x(:),rh(:),rjk(:)
    logical :: viol
    real*8 :: xmean,ymean,a,b,c            !  FIT
    !    real*8, dimension(4) :: x
    !    integer :: partit(4)
    integer :: Npart, partGood,savNstar,fin
    real*8 :: slope,yintercept
    real*8 :: temp_err,temp_rho



    id_err=0

    limit=min(maxknn,nint(0.5*Nele))
    if (mod(limit,4).ne.0) then
       limit=limit+4-mod(limit,4)
    endif

    ! get prefactor for Volume calculation

    if (mod(dimint,2).eq.0) then
       k=dimint/2
       m=1
       do i=1,k
          m=m*i
       enddo
       prefactor=pi**k/(dfloat(m))
    else
       k=(dimint-1)/2
       ms=0.
       do i=1,k
          ms=ms+dlog(dfloat(i))
       enddo
       ns=ms
       do i=k+1,dimint
          ns=ns+dlog(dfloat(i))
       enddo
       prefactor=2.*dexp(ms-ns+k*dlog(4*pi))
    endif

    allocate (Vols(Nele))
    allocate (iVols(Nele))

!    write(12345,*) Nele
!    open(22,file='cacca.dat')
    do i=1,Nele
       Vols(:)=9.9E9
       do j=1,Nele
          if (i.ne.j) then
             if(i.eq.j) then
                id_err=13
                RETURN
             endif
             Vols(j)=prefactor*(gDist(i,j))**dimint
          endif
       enddo
       call HPSORT(Nele,Vols,iVols) !sort Vols, iVols is the permutation
       do j=1,limit
          Nlist(i,j)=iVols(j)
       enddo
       ! get nstar
       viol=.false.
       k=minknn
       n=1
       do while (.not.viol)
          rhg=dfloat(k)/Vols(k)
          dL=dabs(4.*(rhg*(Vols(k)-Vols(3*k/4)-Vols(k/4)))/dfloat(k))
          if (dL.gt.critV(n)) then
             viol=.true.
             cycle
          endif
          n=n+1
          k=k+4
          if (k.gt.limit) viol=.true.
       enddo
       Nstar(i)=k-4 ! ### ha senso?
       if (Nstar(i).lt.minknn) Nstar(i)=minknn ! ### puo' succedere..?

       ! ### ##########################3
       Rho_err(i)=-9.9d99
       rhg=dfloat(Nstar(i))/Vols(Nstar(i)) ! Rho without fit
       Rho(i)=rhg
       savNstar=Nstar(i)
       Npart=4
       fin=Nstar(i)/2
       do while (Npart.le.fin)
          if (mod(Nstar(i),Npart).lt.(Nstar(i)/4)) then
             if (mod(Nstar(i),Npart).ne.0) Nstar(i)=Nstar(i)-mod(Nstar(i),Npart)
             ! get inv rho of the partition
             allocate (x(Npart),rh(Npart),rjk(Npart))
             j=Nstar(i)/Npart
             a=dfloat(j)
             n=0
             do k=1,Npart
                n=n+j
                x(k)=dfloat(k)
                if (k.eq.1) then
                   rh(k)=Vols(n)/a
                else
                   rh(k)=(Vols(n)-Vols(n-j))/a
                endif
             enddo
             xmean=sum(x(:))/dfloat(Npart)
             ymean=sum(rh(:))/dfloat(Npart)
             b=0.
             c=0.
             do k=1,Npart
                a=x(k)-xmean
                b=b+a*(rh(k)-ymean)
                c=c+a*a
             enddo
             a=b/c
             slope=a
             rjfit=ymean-a*xmean
             yintercept=rjfit
             ! Perform jacknife resampling for estimate the error (it includes statistical
             ! error and curvature error)
             do n=1,Npart
                xmean=(sum(x(:))-x(n))/dfloat(Npart-1)
                ymean=(sum(rh(:))-rh(n))/dfloat(Npart-1)
                b=0.
                c=0.
                do k=1,Npart
                   if (k.ne.n) then
                      a=x(k)-xmean
                      b=b+a*(rh(k)-ymean)
                      c=c+a*a
                   endif
                enddo
                a=b/c
                rjk(n)=ymean-a*xmean
             enddo
             rjaver=sum (rjk(:))/dfloat(Npart)
             temp_rho=dfloat(Npart)*rjfit-dfloat(Npart-1)*rjaver
             temp_err=0.
             do k=1,Npart
                temp_err=temp_err+(rjk(k)-rjaver)**2
             enddo
             temp_err=dfloat(Npart-1)*temp_err/dfloat(Npart)
             temp_rho=1./temp_rho
             temp_err=max(temp_rho/sqrt(float(Nstar(i))),temp_rho*temp_rho*sqrt(temp_err))
             if (temp_err.gt.Rho_err(i)) then
                Rho(i)=temp_rho
                Rho_err(i)=temp_err
                partGood=Npart
             endif
             deallocate (x,rh,rjk)
          endif
          Npart=Npart+1
          Nstar(i)=savNstar
       enddo
       ! ### print
!       write (22,'(i6,1x,i3,1x,3(es28.18,1x),i3)') i,Nstar(i),Rho(i),Rho_err(i),rhg,partGood
    enddo
!    close(22)
    deallocate (Vols,iVols)

    ! Filter with neighbours density (Iterative version)
    filter(:)=.false.
    viol=.true.
    niter=0
    do i=1,Nele
       if ((Rho(i).lt.0).or.(Rho_err(i).gt.Rho(i)).or.(Rho(i).gt.1D308)) then
          filter(i)=.true.
       endif
    enddo

    do while (viol)
       niter=niter+1
       viol=.false.
       nfilter=0
       do i=1,Nele
          ! compute average density in the neighborhood and standard dev.
          if (.not.filter(i)) then
             ! ### come puo' essere Rho minore di zero...!?
             if ((Rho(i).lt.0).or.(Rho_err(i).gt.Rho(i))) then
                filter(i)=.true.
                viol=.true.
             else
                a=0.
                b=0.
                n=0
                ! ### questo prob si puo' fare senza loop ma bisogna pensarci (c'e' da considerare i filter...)
                do j=1,Nstar(i)
                   if (.not.filter(Nlist(i,j))) then
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
                         filter(i)=.true.
                         viol=.true.
                      endif
                   endif
                else
                   filter(i)=.true.
                   viol=.true.
                endif
             endif
          else
             nfilter=nfilter+1
          endif
       enddo
    enddo
    return

  end subroutine get_densities
 */

}