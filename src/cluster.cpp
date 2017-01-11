//
// Created by marscher on 1/11/17.
//

#include <udpclust.h>

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
      ! if (Nclus.lt.1) then
      !   Cluster(:)=1
      !   id_err=9
      !   RETURN
      !endif
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
    std::vector<double> Rho_prob(Nele);  // Probability of having maximum Rho
    std::vector<double> Rho_copy(Nele);
    integer, allocatable::iRho(:),ordRho(:)
    integer, allocatable::eb(:,:)    !Border
    elements
    //!integer,allocatable :: survivors(:)  ! ### // member var
    double d, dmin;
    bool extend;


    int Nsurv = self->get_survivors();

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
    std::vector<double> Rho_prob(Nele);

    for (ii = 0; ii < Nsurv; ++ii) {
        i = survivors[ii];
        for (jj = 0; jj < Nsurv; ++jj) {
            j = survivors[jj];
            Rho_prob[i] = Rho_prob[i] - std::log(1.0 + std::exp(2. * (Rho[j] - Rho[i]))
                                                       / std::sqrt(Rho_err[i] * Rho_err[i] + Rho_err[j] * Rho_err[j]));
        }
    }

}