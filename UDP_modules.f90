! ##########################################################
! ###             Modules for UDP clustering             ###
! ##########################################################

! ### I edited ths to make it modular and python compatible
! ### starting on April 5th 2016
! ### I'm removing subroutines that I don't need (graph, noise, etc)
! ### on April 6th 2016.
! ###
! ### 
! ### by Giovanni Pinamonti
! ### PhD student @ SISSA, Trieste
! ### Sector of Molecular and Statistical Biophysics

! 
!  Program that effectuates Cluster analysis in fortran
!  Fully automated clustering by accurate non-parametric density estimation
!  (D'Errico et al, publication pendent, 2015) 
! 
!
! TODO:: Define integer and real types in an accurate way (use KIND in a module)
!           a) normal integer gives errors when il number of elements is big
!           b) investigate the difference between real and real*8


module dp_clustering
  implicit none

contains  

  !############################################
  !### MAIN CLUSTERING ALGORITHM SUBROUTINE ###
  !############################################
  subroutine dp_advance(dist_mat,Cluster,Rho,filter,dimint,Nele,ND,id_err,sensibility)
    implicit none
    !#####################
    integer,intent(inout) :: id_err
    !!Global variables
    real*8,intent(in) :: dist_mat(ND)       ! Distance matrix !###
    integer,intent(in) :: Nele                   ! Number of elements
    integer,intent(in) :: ND                     ! Number of distances
    integer,intent(in) :: dimint                 ! integer of dimset (avoid real*8 calc)
    real*8, intent(in) :: sensibility

    integer,parameter :: maxknn=496   ! maximum number of neighbours to explore
    ! These variables are used in densities calculation and then passed to clustering
    real*8,intent(inout) :: Rho(Nele)        ! Density
    real*8 :: Rho_err(Nele)    ! Density error
!    real*8,allocatable :: dc(:)         ! Dist for density calculation
    logical,intent(inout) :: filter(Nele)  ! Pnt with anomalous dens
    integer :: Nlist(Nele,maxknn) ! Neighbour list within dc
    integer :: Nstar(Nele)   ! N. of NN taken for comp dens
       
    integer,allocatable :: Centers(:) ! Centers of the peaks
    integer,intent(inout) :: Cluster(Nele) ! Cluster ID for the element
    integer :: Nclus                  ! Number of Cluster
    ! These seems to be used for merging 
    real*8,allocatable :: Bord(:,:)     ! Border Densities
    real*8,allocatable :: Bord_err(:,:) ! Border Densities Error

    real*8,allocatable :: cent(:)       ! Center Density
    real*8,allocatable :: cent_err(:)   ! Center Error
    ! Underscore m implies data after automatic mergin
    integer,allocatable :: Cluster_m(:) ! Cluster ID for the element
    integer :: Nclus_m                  ! Number of Cluster merged

    id_err=0
    
    call get_densities(id_err,dist_mat,Nele,dimint,Rho,Rho_err,filter,Nlist,Nstar) ! ### my version
    call clustering(id_err)                      ! get Clusters
    call merging(id_err) ! Generate a new matrix without peaks within border error  
    Cluster=Cluster_m
    !    stop
    return

  contains
 
    subroutine clustering(id_err)
      implicit none
      integer :: id_err
      !! Local variables
      integer :: i,j,k
      integer :: ig
      integer :: l
      logical :: idmax
      real*8,allocatable :: Rho_prob(:)   ! Probability of having maximum Rho
      real*8,allocatable :: Rho_copy(:)
      integer,allocatable :: iRho(:)
      integer,allocatable :: eb(:,:)    ! Border elements
      real*8 :: d,dmin
      logical :: extend
      !!
      id_err=0
      !      allocate (Cluster(Nele))
      !! Identify centers: delta>dc eqv. Max(rho) within dc
      Cluster(:)=0              !
      Nclus=0
      ! Here I compute the probability of having density rho, g_i
      allocate (Rho_prob(Nele))
      Rho_prob(:)=0.
      do i=1,Nele
         if (.not.filter(i)) then
            do j=1,Nele
               if (.not.filter(j)) then
                  Rho_prob(i)=Rho_prob(i)-log(1.d0+exp(2.*(Rho(j)-Rho(i))/sqrt(Rho_err(i)**2+Rho_err(j)**2)))
               endif
            enddo
         endif
      enddo
      ! Now I'm getting the clusters
      do i=1,Nele
         if (.not.filter(i)) then
            idmax=.true.
            j=1
            do while (idmax .and. (j.le.Nstar(i)))
               if ((Rho_prob(i).le.Rho_prob(Nlist(i,j))).and.(.not.filter(Nlist(i,j))))  idmax=.false.
               j=j+1
            enddo
            if (idmax) then
               Nclus=Nclus+1
               Cluster(i)=Nclus
            endif
         endif
      enddo
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

         ! copy of rho (not efficient, but clarifies the code) ### !!!
         allocate (Rho_copy(Nele))
         allocate (iRho(Nele))
         Rho_copy(:)=-Rho_prob(:)

         call HPSORT(Nele,Rho_copy,iRho) ! iRho contains the order in density (iRho(1) is the element with highest Rho...)
         deallocate (Rho_copy)

         ! Assign not filtered
         do i=1,Nele
            j=iRho(i)
            if (.not.filter(j).and.Cluster(j).eq.0) then
               dmin=9.9E29
               do k=1,i-1
                  l=iRho(k) ! ### questa cosa mi da davvero un vantaggio?
                  if (.not.filter(l)) then
                     if (gDist(j,l).le.dmin) then
                        ig=l
                        dmin=gDist(j,l) !
                     endif
                  endif
               enddo
               Cluster(j)=Cluster(ig)
            endif
         enddo


         ! Assign filtered to the same Cluster as its nearest unfiltered neighbour
         ! what happens if all neighbors are filtered
         do i=1,Nele
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

               Cluster(i)=Cluster(ig)
            endif
         enddo
         ! find border densities

         allocate (Bord(Nclus,Nclus),Bord_err(Nclus,Nclus),eb(Nclus,Nclus))
         Bord(:,:)=-9.9D99
         Bord_err(:,:)=0.
         eb(:,:)=0
         
         ! ### questo si puo' fare in maniera migliore? Magari senza dc che viene usato solo qua e non sembra utile.
         !do i=1,Nele ! si puo' fare il loop solo su i filter?
         !   if (.not.filter(i)) then ! what if it is not filtered?
         !      dmin=9.9d99
         !      do j=1,Nele
         !         if (.not.filter(j)) then
         !            if (cluster(j).ne.cluster(i)) then
         !               d=gDist(i,j)
         !               if (d.lt.dmin) then
         !                  dmin=d
         !                  ig=j
         !               endif
         !            endif
         !         endif
         !      enddo
         !      if (dmin.le.dc(i)) then
         !         iref=i
         !         k=0
         !         extend=.true.
         !         do while ( (k.lt.Nstar(i)).and.extend)
         !            k=k+1
         !            if (cluster(Nlist(i,k)).eq.cluster(i)) then
         !               if (gDist(Nlist(i,k),ig).lt.dmin) extend=.false.
         !            endif
         !         enddo
         !         if (extend) then
         !            if (Rho_prob(iref).gt. Bord(cluster(i),cluster(ig))) then
         !               Bord(cluster(i),cluster(ig))=Rho_prob(iref)
         !               Bord(cluster(ig),cluster(i))=Rho_prob(iref)
         !               Bord_err(cluster(i),cluster(ig))=Rho_err(iref)
         !               Bord_err(cluster(ig),cluster(i))=Rho_err(iref)
         !               eb(cluster(i),cluster(ig))=iref
         !               eb(cluster(ig),cluster(i))=iref
         !            endif
         !         endif
         !      endif !dmin.le.dc(i)
         !   endif !filter
         !enddo ! i=1,Nele
         ! ######################3
         ! ### my version
         do i=1,Nele ! si puo' fare il loop solo su i filter?
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
            do k=1,Nstar(i)
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
    !
    !
    subroutine merging(id_err)
      implicit none
      integer :: id_err
      !! Local variables
      integer :: i,j,k,n,l,alive,dead,niter ! niter is useless
      logical :: Survive (Nclus)
      logical :: change
      integer :: Nbarr
      integer,allocatable :: Bcorr (:,:)
      real*8,allocatable :: Barrier (:)
      real*8,allocatable :: Barrier_err (:)
      integer,allocatable :: iBarrier(:)
      real*8 :: c1,c2,b12
      integer,allocatable :: M2O(:) !conversion from merged to original cluster number
      integer :: O2M(Nclus) ! Conversion from original cluster number to its equivalent in merged

      id_err=0
      Nbarr=(Nclus*Nclus-Nclus)/2 ! n. of contacts between clusters
      allocate (Barrier(Nbarr))
      allocate (Barrier_err(Nbarr))
      allocate (iBarrier(Nbarr))
      allocate (Bcorr(Nbarr,2))
      n=0
      do i=1,Nclus-1
         do j=i+1,Nclus
            n=n+1
            Barrier(n)= Bord(i,j)
            Barrier_err(n)= Bord_err(i,j)
            Bcorr(n,1)=i
            Bcorr(n,2)=j
         enddo
      enddo
      if (Nbarr.gt.1) then
         call HPSORT(Nbarr,Barrier,iBarrier)
      else
         iBarrier(1)=1 !e gli altri?
      endif
      Survive(:)=.true.
      change=.true.
      allocate (Cluster_m(Nele))
      Cluster_m(:)=Cluster(:)
      niter=0
      do while (change)
         niter=niter+1
         change=.false.
         mdo:  do n=Nbarr,1,-1
            k=iBarrier(n)
            i=Bcorr(k,1)
            j=Bcorr(k,2)
            if ((Bord(i,j).gt.0.).and.(i.ne.j)) then
               if (Survive(i).and.Survive(j)) then
                  c1=cent(i)-sensibility*cent_err(i)
                  c2=cent(j)-sensibility*cent_err(j)
                  b12=Bord(i,j)+sensibility*Bord_err(i,j)
                  if ((c1.lt.b12).or.(c2.lt.b12)) then
                     change=.true.
                     Bord(i,j)=0.
                     Bord_err(i,j)=0.
                     if (c1.gt.c2) then
                        alive=i
                        dead=j
                     else
                        alive=j
                        dead=i
                     endif
            !forse queste due linee sono la differenza con alex
                     Bord(alive,alive)=0.
                     Bord_err(alive,alive)=0.
                     Survive(dead)=.false.
                     do k=1,Nclus
                        if (Survive(k)) then
                           if (Bord(i,k).gt.Bord(j,k)) then
                              Bord(alive,k)=Bord(i,k)
                              Bord(k,alive)=Bord(k,i)
                              Bord_err(alive,k)=Bord_err(i,k)
                              Bord_err(k,alive)=Bord_err(k,i)
                           else
                              Bord(alive,k)=Bord(j,k)
                              Bord(k,alive)=Bord(k,j)
                              Bord_err(alive,k)=Bord_err(j,k)
                              Bord_err(k,alive)=Bord_err(k,j)
                           endif
                        endif
                     enddo
                     do l=1,Nele
                        if (Cluster_m(l).eq.dead) Cluster_m(l)=alive
                     enddo
                     if (n.gt.1) then
                        do l=n-1,1,-1
                           k=iBarrier(l)
                           if (Bcorr(k,1).eq.dead) Bcorr(k,1)=alive
                           if (Bcorr(k,2).eq.dead) Bcorr(k,2)=alive
                        enddo
                     endif
                     exit mdo
                  endif
               endif
            endif
         enddo mdo
      enddo

      ! get dictionary

      Nclus_m=0
      do i=1,Nclus
         if (Survive(i)) Nclus_m=Nclus_m+1
      enddo
      allocate (M2O(Nclus_m))
      n=0
      O2M(:)=-1
      do i=1,Nclus
         if (Survive(i)) then
            n=n+1
            M2O(n)=i
            O2M(i)=n
         endif
      enddo

      ! get survival characteristics
      if (Nclus_m.gt.1) then
!         allocate (Bord_m(Nclus_m,Nclus_m),Bord_err_m(Nclus_m,Nclus_m))
!         allocate (cent_m(Nclus_m),cent_err_m(Nclus_m))
!         allocate (Centers_m(Nclus_m))
         do i=1,Nele
            Cluster_m(i)=O2M(Cluster_m(i))
         enddo
!         do i=1,Nclus_m
!            do j=i+1,Nclus_m
!               Bord_m(i,j)=Bord(M2O(i),M2O(j)) 
!               Bord_err_m(i,j)=Bord_err(M2O(i),M2O(j)) 
!               Bord_m(j,i)=Bord(M2O(i),M2O(j)) 
!               Bord_err_m(j,i)=Bord_err(M2O(i),M2O(j)) 
!            enddo
!            cent_m(i)=cent(M2O(i))
!            cent_err_m(i)=cent_err(M2O(i))
!            Centers_m(i)=Centers(M2O(i))
!         enddo
      else
         id_err=10
         Cluster_m(:)=1
      endif

      return
    end subroutine merging
    !
    ! This function allows to get the distances in a matrix like style
    real*8 function gDist(i,j)
      integer :: i,j,k,l,m
      l=max(i,j)
      m=min(i,j)
      k=(m-1)*Nele-(m*m+m)/2+l
      gDist=dist_mat(k)
      return
    end function gDist
    !
    !
  end subroutine dp_advance

  subroutine get_densities(id_err,dist_mat,Nele,dimint,Rho,Rho_err,filter,Nlist,Nstar)
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
    !### cambiato da Alex
    !    if (mod(dimint,2).eq.0) then
    !       k=dimint/2
    !       m=1
    !       do i=1,k
    !          m=m*i
    !       enddo
    !       prefactor=pi**k/(dfloat(m))
    !    else
    !       k=(dimint-1)/2
    !       m=1
    !       do i=1,k
    !          m=m*i
    !       enddo
    !       n=m
    !       do i=k+1,dimint
    !          n=n*i
    !       enddo
    !       prefactor=2.*dfloat(m)*(4.*pi)**k/(dfloat(n))
    !    endif
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

    write(12345,*) Nele
!    open(22,file='cacca.dat')
    do i=1,Nele
       Vols(:)=9.9E9
       do j=1,Nele
          if (i.ne.j) then
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

       ! ### kadd funge da boolean. Qui sto aggiungendo vicini che sono 
       ! ### a una distanza computazionalmente indistinguibile al mio calcolo
       ! ### (ma ha senso sta roba?) Comunque nel nuovo prog di alex non c'e'
       !kadd=1
       !if (Nstar(i).eq.limit) kadd=0
       !do while (kadd.eq.1)
       !   if ((abs(gDist(i,Nlist(i,Nstar(i)))-gDist(i,Nlist(i,Nstar(i)+1)))).lt.9.99D-99) then
       !      Nstar(i)=Nstar(i)+1
       !      if (Nstar(i).eq.limit) kadd=0
       !   else
       !      kadd=0
       !   endif
       !enddo
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


  contains
    !
    real*8 function gDist(i,j)
      integer :: i,j,k,l,m
      l=max(i,j)
      m=min(i,j)
      k=(m-1)*Nele-(m*m+m)/2+l
      gDist=dist_mat(k)
      return
    end function gDist
    !
  end subroutine get_densities
  
  
  SUBROUTINE HPSORT(N,RA,s_order)
    implicit none
    integer N,s_order(N)
    real*8 RA(N)
    integer L,IR,I,J,sord
    real*8 RRA
    do i=1,n
       s_order(i)=i
    enddo
    L=N/2+1
    IR=N
    !The index L will be decremented from its initial value during the
    !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
    !will be decremented from its initial value down to 1 during the
    !"retirement-and-promotion" (heap selection) phase.
10  continue
    if(L > 1)then
       L=L-1
       RRA=RA(L)
       sord=s_order(L)
    else
       RRA=RA(IR)
       sord=s_order(IR)
       RA(IR)=RA(1)
       s_order(ir)=s_order(1)
       IR=IR-1
       if(IR.eq.1)then
          RA(1)=RRA
          s_order(1)=sord
          return
       end if
    end if
    I=L
    J=L+L
20  if(J.le.IR)then
       if(J < IR)then
          if(RA(J) < RA(J+1))  J=J+1
       end if
       if(RRA < RA(J))then
          RA(I)=RA(J)
          s_order(i)=s_order(j)
          I=J; J=J+J
       else
          J=IR+1
       end if
       goto 20
    end if
    RA(I)=RRA
    s_order(I)=sord
    goto 10
  END subroutine HPSORT

end module dp_clustering
