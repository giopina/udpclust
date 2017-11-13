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
  subroutine dp_advance(dist_mat,Cluster,Rho,Rho_err,filter,dimint,Nlist,Nstar,Nele,id_err,sensibility,maxknn)
    !$ use omp_lib
    implicit none
    integer,intent(in) :: maxknn   ! maximum number of neighbours to explore
    integer,intent(inout) :: id_err
    !!Global variables
    real*8,intent(in) :: dist_mat(Nele,maxknn)       ! Distance matrix 
    integer,intent(in) :: Nele                   ! Number of elements
    integer,intent(in) :: dimint                 ! integer of dimset (avoid real*8 calc)
    real*8, intent(in) :: sensibility

    ! These variables are used in densities calculation and then passed to clustering
    real*8,intent(in) :: Rho(Nele)        ! Density
    real*8,intent(in) :: Rho_err(Nele)    ! Density error
    logical,intent(in) :: filter(Nele)  ! Pnt with anomalous dens
    integer,intent(in) :: Nlist(Nele,maxknn) ! Neighbour list within dc
    integer,intent(in) :: Nstar(Nele)   ! N. of NN taken for comp dens
    integer :: survivors(Nele) 
    integer :: Nsurv           
       
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
    !call get_densities(id_err,dist_mat,Nele,dimint,Rho,Rho_err,filter,Nlist,Nstar,maxknn) 
    call clustering(id_err)                      ! get Clusters
    if(Nclus.gt.1) then
       if(sensibility.gt.0.0) then
          call merging(id_err) ! Generate new clusters without peaks within border error  
          Cluster=Cluster_m
       endif
    endif
    return

  contains
    subroutine get_survivors()
      integer :: i
      Nsurv=0
      do i=1,Nele
         if (.not.filter(i)) then
            Nsurv=Nsurv+1
            survivors(Nsurv)=i
         endif
      enddo
    end subroutine get_survivors
    

    subroutine clustering(id_err)
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
      real*8 :: d,dmin
      logical :: extend      
      logical :: viol,newass
      integer :: Nnoass

      id_err=0
      !! Identify centers: delta>dc eqv. Max(rho) within dc
      Cluster(:)=0              !
      Nclus=0
      ! Here I compute the probability of having density rho, g_i
      allocate (Rho_prob(Nele))
      Rho_prob(:)=0.
      call get_survivors() 
      !$OMP PARALLEL DO private(i,j) shared(survivors,Rho_prob,Rho,Nsurv,Rho_err)
      do ii=1,Nsurv
         i=survivors(ii)
         do jj=1,Nsurv
            j=survivors(jj)
            Rho_prob(i)=Rho_prob(i)-log(1.d0+exp(2.*(Rho(j)-Rho(i))/sqrt(Rho_err(i)**2+Rho_err(j)**2)))
         enddo
      enddo
      !$OMP END PARALLEL DO

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

      ! Now I'm getting the clusters
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


      allocate (Centers(Nclus))
      do i=1,Nele
         if (Cluster(i).ne.0) then
            Centers(Cluster(i))=i
         endif
      enddo
      if (Nclus.lt.1) then
         ! TODO: I'm actually not sure this can happen
         Cluster(:)=-1
         id_err=9
         return
      endif
      if (Nclus.eq.1) then      
         return
      endif
   
      ! Assign not filtered
      ! TODO: change it with survivors
      do j=1,Nele
         i=iRho(j)
         if ((.not.filter(i)).and.Cluster(i).eq.0) then
            ig=-1
            dmin=9.9d99
            do k=1,Nstar(i)
               l=Nlist(i,k)
               if (.not.filter(l)) then
                  if(Rho_prob(i).lt.Rho_prob(l)) then ! ### TODO check this
                     if (dist_mat(i,k).le.dmin) then
                        ig=l
                        dmin=dist_mat(i,k) 
                     endif
                  endif
               endif
            enddo
            if (ig.eq.-1) then
               id_err=12
               RETURN
            else
               Cluster(i)=Cluster(ig)
            endif
         endif
      enddo

      ! Assign filtered to the same Cluster as its nearest unfiltered neighbour
      ! what happens if all neighbors are filtered?
      ! this can happen, at least in principle. What should I do then?
      ! option 1: assign the point to the same cluster of the closest filtered point (you should do this in the end, but it will depend on the order you consider the points...
      viol=.false.
      do i=1,Nele
         if (Cluster(i).eq.0) then
            ig=-1
            if(.not.filter(i)) then
               id_err=12
               RETURN
            endif
            dmin=9.9d99
            do k=1,maxknn ! find the min d in not filt elements
               l=Nlist(i,k)
               if ((Cluster(l).ne.0).and.(.not.filter(l))) then
                  d=dist_mat(i,k)
                  if (d.le.dmin) then
                     dmin=d
                     ig=l
                  endif
               endif
            enddo
            if (ig.ne.-1) then
               Cluster(i)=Cluster(ig)
            else
               viol=.true.
            endif
         endif
      enddo
      do while (viol) 
         ! TODO: change this so does not depend on the order
         viol=.false.
         NEWASS=.FALSE.
         Nnoass=0
         do i=1,Nele
            if (Cluster(i).eq.0) then
               dmin=9.9d99
               ig=-1
               ! if (dmin.gt.9.8d99) then !### what if all the neighbours are filtered?
               ! do k=1,Nstar(i) ! find the min d in filter elements
               do k=1,maxknn ! find the min d in filter elements
                  j=Nlist(i,k)
                  if ((Cluster(j).ne.0)) then
                     d=dist_mat(i,k)
                     if (d.le.dmin) then
                        dmin=d
                        ig=j
                     endif
                  endif
               enddo
               if (ig.ne.-1) then
                  Cluster(i)=Cluster(ig)
                  NEWASS=.TRUE.
               else
                  viol=.true.
                  Nnoass=Nnoass+1
               endif
            endif ! endif cluster(i)==0
         enddo
         ! check if this stopped or got stuck. It should not happen
         IF(.NOT.NEWASS) THEN
            WRITE(*,*) 'ERROR: cannot assign some filtered points to clusters'
            id_err=12
            RETURN

         endif
      enddo

      ! find border densities
      allocate (Bord(Nclus,Nclus),Bord_err(Nclus,Nclus),eb(Nclus,Nclus))
      Bord(:,:)=-9.9D99
      Bord_err(:,:)=0.
      eb(:,:)=0

      do i=1,Nele ! si puo' fare il loop solo su i filter?
         if (filter(i)) CYCLE
         ig=-1
         dmin=9.9d99
         do k=1,Nstar(i)
            l=Nlist(i,k)
            if (filter(l)) CYCLE
            if (cluster(l).eq.cluster(i)) CYCLE
            d=dist_mat(i,k)
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
            l=Nlist(i,k)
            if(filter(l)) CYCLE
            if(cluster(l).ne.cluster(i)) CYCLE
            d=9.9d99
            do jj=1,maxknn
               j=Nlist(l,jj)
               if (j.eq.ig) THEN
                  d=dist_mat(l,jj)
                  EXIT
               ENDIF
            enddo
            if (d.lt.dmin) then 
               extend=.false.
               EXIT
            endif
         enddo
         if (extend) then
            if (Rho_prob(i).gt. Bord(cluster(i),cluster(ig))) then
               Bord(cluster(i),cluster(ig))=Rho_prob(i)
               Bord(cluster(ig),cluster(i))=Rho_prob(i)
               Bord_err(cluster(i),cluster(ig))=Rho_err(i)
               Bord_err(cluster(ig),cluster(i))=Rho_err(i)
               eb(cluster(i),cluster(ig))=i
               eb(cluster(ig),cluster(i))=i
            endif
         endif
      enddo ! i=1,Nele

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
            if (.not.filter(j)) then !
               if ((cluster(j).eq.i).and.((Rho(j)-Rho_err(j)).gt.(cent(i)-cent_err(i)))) then !
                  cent(i)=Rho(j)
                  cent_err(i)=Rho_err(j) !
               endif
            endif
         enddo
      enddo
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
         iBarrier(1)=1 
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

      if (Nclus_m.eq.1) then
         Cluster_m(:)=1
         return
      endif
      if (Nclus_m.lt.1) then
         ! TODO: not sure this can happen
         id_err=10
         Cluster_m(:)=-1
         return
      endif
      ! get survival characteristics      
      do i=1,Nele
         Cluster_m(i)=O2M(Cluster_m(i))
      enddo
      return
    end subroutine merging
    !
  end subroutine dp_advance

  subroutine get_densities(id_err,dist_mat,Nele,dimint,Rho,Rho_err,filter,Nlist,Nstar,maxknn)
    use critfile
    !$ use omp_lib
    implicit none

    integer,intent(in) :: maxknn   ! maximum number of neighbours to explore
    integer,parameter :: minknn=8     ! minimum number of neighbours to explore
    !!Global variables
    real*8,intent(in) :: dist_mat(Nele,maxknn)        !
    integer,intent(in) :: Nele                   ! Number of elements
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
    integer :: kadd
    integer :: niter,nfilter
    !real*8,allocatable :: Vols(:)
    real*8 :: Vols(Nele,maxknn)
    real*8, parameter :: pi=3.14159265359
    real*8 :: prefactor
    !real*8 :: rhg,dL,rjaver,rjfit
    real*8 :: rhg,Dk,rjaver,rjfit
    real*8,allocatable :: x(:),rh(:),rjk(:)
    logical :: viol
    real*8 :: xmean,ymean,a,b,c            !  FIT
    real*8 :: xsum,ysum,x2sum,xysum
    integer :: Npart, partGood,savNstar,fin
    real*8 :: slope,yintercept
    real*8 :: temp_err,temp_rho
    real*8 :: dimreal

    id_err=0


    ! ### TODO: check this. Now that I'm not dividing by 4 limit makes less sense...
    limit=min(maxknn,nint(0.5*Nele))
    if (mod(limit,4).ne.0) then
       limit=limit+4-mod(limit,4)
    endif
    !write(*,*) limit,maxknn ! ### what is the meaning of this limit ?

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
    !write(*,*) "prefactor", prefactor
    dimreal=FLOAT(dimint)

    Vols=prefactor*dist_mat**dimint
    
    do i=1,Nele
       ! ### get nstar     
       do k=minknn,maxknn         
          j=Nlist(i,k)
          Dk= -2*( log(Vols(i,k)) + log(Vols(j,k)) - 2*log(Vols(i,k)+Vols(j,k)) + log(4.) )
          if (Dk.gt.23.928) then
             exit
          endif
          if (k.gt.limit) exit
       enddo
       Nstar(i)=k-1 ! ### ha senso?
       if (Nstar(i).lt.minknn) Nstar(i)=minknn ! ### puo' succedere..?
    write(12345,*) i,Nstar(i)
    enddo

    do i=1,Nele
       ! ### get effective rho
       Rho_err(i)=-9.9d99
       rhg=dfloat(Nstar(i))/Vols(i,Nstar(i)) ! Rho without fit
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
                   rh(k)=Vols(i,n)/a
                else
                   rh(k)=(Vols(i,n)-Vols(i,n-j))/a
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
             xsum=sum(x(:))
             ysum=sum(rh(:))
             x2sum=sum(x**2)
             xysum=sum(x*rh)
             do n=1,Npart
                xmean=(xsum-x(n))/dfloat(Npart-1)
                ymean=(ysum-rh(n))/dfloat(Npart-1)
                c=x2sum-x(n)**2
                c=c-xmean*xmean*(Npart-1)
                b=xysum-x(n)*rh(n)
                b=b-xmean*ymean*(Npart-1)
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
    enddo

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
