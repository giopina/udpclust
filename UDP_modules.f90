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
  subroutine dp_advance(dist_mat,Cluster,Rho,Rho_err,Nlist,Nstar,Nele,id_err,sensibility,maxknn)
    !$ use omp_lib
    implicit none
    integer,intent(in) :: maxknn   ! maximum number of neighbours to explore
    integer,intent(inout) :: id_err
    !!Global variables
    real*8,intent(in) :: dist_mat(Nele,maxknn)       ! Distance matrix 
    integer,intent(in) :: Nele                   ! Number of elements
    real*8, intent(in) :: sensibility

    ! These variables are used in densities calculation and then passed to clustering

    real*8,intent(in) :: Rho(Nele)        ! LOGARITHM of Density
    real*8,intent(in) :: Rho_err(Nele)    ! error on LOGARITHM of density
    integer,intent(in) :: Nlist(Nele,maxknn) ! Neighbour list within dc
    integer,intent(in) :: Nstar(Nele)   ! N. of NN taken for comp dens
       
    integer,allocatable :: Centers(:) ! Centers of the peaks
    integer,intent(inout) :: Cluster(Nele) ! Cluster ID for the element
    integer :: Nclus                  ! Number of Cluster
    ! These seems to be used for merging 
    real*8,allocatable :: Bord(:,:)     ! Border Densities
    real*8,allocatable :: Bord_err(:,:) ! Border Densities Error

    real*8,allocatable :: cent(:)       ! Center Density
    real*8,allocatable :: cent_err(:)   ! Center Error
    ! Underscore m implies data after automatic mergin
    !integer,allocatable :: Cluster_m(:) ! Cluster ID for the element
    integer :: Cluster_m(Nele) ! Cluster ID for the element
    integer :: Nclus_m                  ! Number of Cluster merged
    logical :: verbose=.FALSE.
    
    id_err=0
    if(verbose) write(*,*) 'clustering'
    call clustering(id_err)                      ! get Clusters
    if(Nclus.gt.1) then
       if(sensibility.gt.0.0) then
          if(verbose) write(*,*) 'finding borders'
          call find_borders(id_err)
          if(verbose) write(*,*) 'merging clusters'
          call merging(id_err) ! Generate new clusters without peaks within border error
          !Cluster=Cluster_m
          !deallocate(Cluster_m)
       endif
    endif
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
      integer,allocatable :: iRho(:),ordRho(:)
      real*8 :: d,dmin

      id_err=0
      !! Identify centers: delta>dc eqv. Max(rho) within dc
      Cluster(:)=0              !
      Nclus=0
      ! Here I compute the probability of having density rho, g_i
      allocate (Rho_prob(Nele))
      Rho_prob(:)=0.
      !
      do i=1,Nele
         Rho_prob(i)=Rho(i)-Rho_err(i)
      enddo
      !write(*,*) 'Rho_prob computed'
      !
      ! copy of rho (not efficient, but clarifies the code) ### !!!
      allocate (Rho_copy(Nele))
      !write(*,*) Nele
      allocate (iRho(Nele))
      !write(*,*) iRho
      Rho_copy(:)=-Rho_prob(:)
      call HPSORT(Nele,Rho_copy,iRho) ! iRho contains the order in density (iRho(1) is the element with highest Rho...)
      deallocate (Rho_copy)      
      allocate (ordRho(Nele))
      do i=1,Nele
         ordRho(iRho(i))=i                 ! ordRho is the complementary of iRho. Given an element, ordRho returns its order in density
      enddo

      ! Now I'm getting the clusters
      open(123,file="my_centers_premerge.dat")
      do i=1,Nele
         idmax=.true.
         j=1
         do while (idmax .and. (j.le.Nstar(i)))
            if ((ordRho(i).gt.ordRho(Nlist(i,j))))  idmax=.false.
            j=j+1
         enddo
         if (idmax) then
            j=1
            do while  (idmax .and. (j.lt.ordRho(i)))
               k=1
               l=iRho(j)
               do while (idmax.and.(k.le.Nstar(l)))
                  if (Nlist(l,k).eq.i) idmax=.false.
                  k=k+1
               enddo
               j=j+1
            enddo
            if (idmax) then
               Nclus=Nclus+1
               Cluster(i)=Nclus
            endif
         endif
         write(123,*) Cluster(i)
      enddo
      close(123)

     

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

      ! Assign points to clusters
      do j=1,Nele
         i=iRho(j)
         if (Cluster(i).eq.0) then
            ig=-1
            dmin=9.9d9
            !do k=1,Nstar(i)
            do k=1,maxknn
               l=Nlist(i,k)
               if(Rho_prob(i).lt.Rho_prob(l)) then
                  if (dist_mat(i,k).le.dmin) then
                     ig=l
                     dmin=dist_mat(i,k) 
                  endif
               endif
            enddo
            if (ig.eq.-1) then
               ! ### TODO: check this
               write(*,*) '*** ig=-1 *** unassigned points!'
               id_err=12
               RETURN
            else
               Cluster(i)=Cluster(ig)
            endif
         endif
      enddo

      open(1234,file='my_cluster_premerge.dat')
      do i=1,Nele
         write(1234,*) Cluster(i)
      enddo
      close(1234)
      
      deallocate (Rho_prob)
      deallocate (iRho)
      deallocate(ordRho)


      return
    end subroutine clustering
    !
    subroutine find_borders(id_err)
      integer :: id_err
      !! Local variables
      integer :: i,j,k
      integer :: jj
      integer :: ig
      integer :: l,icl
      integer :: eb(Nclus,Nclus)    ! Border elements
      real*8 :: d
      real*8 :: dmin(Nclus)
      integer :: imin(Nclus)
      
      logical :: extend      
      logical :: viol,newass
      integer :: Nnoass

      ! find border densities
      !allocate (Bord(Nclus,Nclus),Bord_err(Nclus,Nclus),eb(Nclus,Nclus))
      allocate (Bord(Nclus,Nclus),Bord_err(Nclus,Nclus))
      Bord(:,:)=-9.9D9
      Bord_err(:,:)=0.
      eb(:,:)=0
      !allocate (dmin(Nclus),imin(Nclus))
      do i=1,Nele
         ig=-1
         dmin(:)=9.9d9 !### controlla questo
         imin(:)=-1
         do k=1,Nstar(i)
            l=Nlist(i,k)
            if (cluster(l).eq.cluster(i)) CYCLE
            d=dist_mat(i,k)
            if (d.lt.dmin(cluster(l))) then
               dmin(cluster(l))=d
               imin(cluster(l))=l
               ig=1
            endif
         enddo
         !if (dmin.gt.9.8d99) CYCLE
         if (ig.eq.-1) CYCLE
         extend=.true.
         !if (ig.eq.-1) CYCLE
         !   id_err=12
         !   RETURN
         !endif
         do k=1,Nstar(i)
            l=Nlist(i,k)
            if(cluster(l).ne.cluster(i)) CYCLE
            !d=9.9d99
            d=9.9d9
            do jj=1,maxknn
               j=Nlist(l,jj)
               ig=imin(cluster(j))
               if (j.eq.ig) THEN
                  d=dist_mat(l,jj)
                  if(d.lt.dmin(cluster(j))) imin(cluster(j))=-1
               ENDIF
            enddo
            !if (d.lt.dmin) then 
            !   extend=.false.
            !   EXIT
            !endif
            ! ### TODO: here one might want to add a check to see if the cycle over the NN of i can be stopped
         enddo ! k
         do icl=1,Nclus
            ig=imin(icl)
            if (ig.gt.-1) then
               if ((Rho(i)-Rho_err(i)).gt. Bord(cluster(i),icl)) then
                  Bord(cluster(i),icl)=Rho(i)-Rho_err(i)
                  Bord(icl,cluster(i))=Rho(i)-Rho_err(i)
                  Bord_err(cluster(i),icl)=Rho_err(i)
                  Bord_err(icl,cluster(i))=Rho_err(i)
                  eb(cluster(i),icl)=i
                  eb(icl,cluster(i))=i
               endif
            endif
         enddo !icl
      enddo ! i=1,Nele
      open (22,file="my_borders_premerge.dat")
      do i=1,Nclus-1
         do j=i+1,Nclus
            if (eb(i,j).ne.0) then
               Bord(i,j)=Rho(eb(i,j))
               Bord(j,i)=Rho(eb(i,j))
            else
               Bord(i,j)=0.
               Bord(j,i)=0.
            endif
            write (22,'(i6,1x,i6,1x,i6,1x,es18.10,1x,es18.10)') i,j,eb(i,j),Bord(i,j),Bord_err(i,j)
         enddo
      enddo
      close(22)
      ! Info per graph pre automatic merging
      allocate (cent(Nclus))
      allocate (cent_err(Nclus))
      do i=1,Nclus
         cent(i)=Rho(Centers(i))
         cent_err(i)=Rho_err(Centers(i))
      enddo
      !deallocate(dmin)
      !deallocate(imin)
      !deallocate(eb)
      return
    end subroutine find_borders
    !
    subroutine merging(id_err)
      implicit none
      integer :: id_err
      !! Local variables
      integer :: i,j,k,n,l,alive,dead,niter ! niter is useless
      logical :: Survive (Nclus)
      logical :: change
      integer :: Nbarr
      !integer :: Nbarr=(Nclus*Nclus-Nclus)/2
      integer :: Bcorr ((Nclus*Nclus-Nclus)/2,2)
      real*8 :: Barrier ((Nclus*Nclus-Nclus)/2)
      real*8 :: Barrier_err ((Nclus*Nclus-Nclus)/2)
      integer :: iBarrier((Nclus*Nclus-Nclus)/2)
      real*8 :: c1,c2,b12,b1,b2,f1,f2,f12
      integer,allocatable :: M2O(:) !conversion from merged to original cluster number
      integer :: O2M(Nclus) ! Conversion from original cluster number to its equivalent in merged

      id_err=0
      Nbarr=(Nclus*Nclus-Nclus)/2 ! n. of contacts between clusters
      !allocate (Barrier(Nbarr))
      !allocate (Barrier_err(Nbarr))
      !allocate (iBarrier(Nbarr))
      !allocate (Bcorr(Nbarr,2))
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
      !allocate (Cluster_m(Nele))
      Cluster_m(:)=Cluster(:) !### this is probably useless!
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
                  c1=(cent(i)-Bord(i,j))/(Bord_err(i,j)+cent_err(i))
                  c2=(cent(j)-Bord(i,j))/(Bord_err(i,j)+cent_err(j))

                  !c1=cent(i)-sensibility*cent_err(i)
                  !c2=cent(j)-sensibility*cent_err(j)
                  !b12=Bord(i,j)+sensibility*Bord_err(i,j)
                  b12=min(c1,c2)
                  f1=cent(i)
                  f2=cent(j)
                  f12=Bord(i,j)
                  b1=abs(f12-f1)
                  b2=abs(f12-f2)
                  if ((b12.lt.sensibility).or.(b1.lt.0.0).or.(b2.lt.0.0)) then
                  !if ((c1.lt.b12).or.(c2.lt.b12)) then
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
      enddo !while change
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

      Cluster=Cluster_m
      deallocate(Bord)
      deallocate(Bord_err)
      deallocate(cent)
      deallocate(cent_err)
      !deallocate(Barrier)
      !deallocate(Barrier_err)
      !deallocate(iBarrier)
      !deallocate(Bcorr)

      deallocate(M2O)
      
      return
    end subroutine merging
    !
  end subroutine dp_advance


  subroutine get_densities(id_err,dist_mat,Nele,dim,Rho,Rho_err,Nlist,Nstar,maxknn)
    !$ use omp_lib
    implicit none

    integer,intent(in) :: maxknn   ! maximum number of neighbours to explore
    real*8,intent(in) :: dist_mat(Nele,maxknn)        !
    integer,intent(in) :: Nele                   ! Number of elements

    integer,intent(in) :: dim                 ! integer of dimset (avoid real*8 calc)
    real*8,intent(inout) :: Rho(Nele)        ! Density
    real*8,intent(inout) :: Rho_err(Nele)    ! Density error
    integer,intent(inout) :: Nlist(Nele,maxknn) ! Neighbour list within dc ### dimension 2 maybe can be reduced but maybe not (before was=limit)
    integer,intent(inout) :: Nstar(Nele)   ! N. of NN taken for comp dens
    
    !internal variables
    integer :: i
    real*8 :: F
    integer :: id_err
    real*8 :: Vols(Nele,maxknn)
    real*8 :: min_F
    id_err=0

    Vols=prefactor(dim)*dist_mat**dim
    call get_k(id_err,Vols,Nele,dim,Nlist,Nstar,maxknn)

    min_F=1.0
    do i=1,Nele
       F=free_energy(Nstar(i),Vols(i,:))
       !Rho(i)=exp(-F)
       Rho(i)=F ! I'm computing the log of density...
       if (F.lt.min_F) min_F=F
       Rho_err(i)=dsqrt(dfloat(4*Nstar(i)+2)/dfloat(Nstar(i)*(Nstar(i)-1))) ! I can do this outside!
    enddo
    Rho(:)=Rho(:)-min_F+1.0 ! This way the minimum Rho is equal to 1. (for hierarchies it does not matter so much)

    return

  end subroutine get_densities

  subroutine get_k(id_err,Vols,Nele,dim,Nlist,Nstar,maxknn)
    !$ use omp_lib
    implicit none

    integer,intent(in) :: maxknn   ! maximum number of neighbours to explore
    integer,parameter :: minknn=4     ! minimum number of neighbours to explore
    !!Global variables
    integer,intent(in) :: Nele                   ! Number of elements
    integer,intent(in) :: dim                 ! integer of dimset (avoid real*8 calc)
    integer,intent(inout) :: Nlist(Nele,maxknn) ! Neighbour list within dc ### dimension 2 maybe can be reduced but maybe not 
    integer,intent(inout) :: Nstar(Nele)   ! N. of NN taken for comp dens
    real*8, intent(in) :: Vols(Nele,maxknn) ! the ordered spherical volumes corresponding to nn radii
    integer :: id_err
    !! Local variables
    integer :: limit
    integer :: i,j,k
    real*8 :: Dk

    id_err=0
    limit=min(maxknn,int(Nele/2))
    do i=1,Nele
       ! ### get nstar     
       do k=minknn,maxknn-1
          j=Nlist(i,k+1)
          Dk= -2*k*( log(Vols(i,k)*Vols(j,k)/(Vols(i,k)+Vols(j,k))**2) + log(4.) )
          if (Dk.gt.23.928) then
             exit
          endif
          if (k.gt.limit-1) exit
       enddo
       Nstar(i)=k ! ### ha senso?
       if (Nstar(i).lt.minknn) Nstar(i)=minknn ! ### puo' succedere..?
    enddo
    return

  end subroutine get_k

  real*8 function free_energy(Nstar,VV)
    ! This function computes the optimal "Free energy" aka -log(Rho)
    ! using the procedure described in Rodriguez et al. JCTC 2018
    !
    implicit none
    integer, intent(in) :: Nstar
    integer :: id_err
    !! Local variables
    integer :: i,j,k,m
    integer :: niter
    real*8,intent(in) :: VV(:)
    real*8 :: dL
    real*8 :: vi(Nstar)
    logical :: viol
    real*8 :: a,F
    real*8 :: H(2,2),Hinv(2,2)
    real*8 :: func,sigma,t,jf
    real*8 :: tt,gb,ga,sa,sb,a_err
    real*8 :: stepmax !! This variable controls the maximum variation on the log(rho) accepted during the optimization process
    real*8 :: r
    ! Variables Specific of variable kernel
    real*8 :: g
    ! Variable that accounts for vi=0 --> then use kNN instead for density (error
    ! still PAk)
    logical :: kNN
    !real*8 :: Rho_err

    ! ### Here I'm computing the volumes of the spherical shells between neighbors
    Vi(1)=VV(1)
    kNN=.false.
    do j=2,Nstar
       Vi(j)=VV(j)-VV(j-1)
       if ((Vi(j).eq.0).and.(.not.kNN)) then !Can this really happen? What does it mean?
          write (*,*) "Point density estimated with k-NN" ! Am I really not doing anything if there are 2 neighbors at the same distance??
          kNN=.true.
       endif
    enddo
    
    F=log(float(Nstar)/VV(Nstar)) ! Starting guess for F
    a=0.  !
    if (.not.kNN) then
       stepmax=0.1*abs(F)
       !H=Hessian(a,b,Nstar)
       gb=float(Nstar)
       ga=float(Nstar+1)*float(Nstar)/2.
       H(1,1)=0.
       H(1,2)=0.
       H(2,2)=0.
       do j=1,Nstar
          jf=float(j)
          tt=vi(j)*exp(F+a*jf)
          gb=gb-tt
          ga=ga-jf*tt
          H(1,1)=H(1,1)-tt
          H(1,2)=H(1,2)-jf*tt
          H(2,2)=H(2,2)-jf*jf*tt
       enddo
       H(2,1)=H(1,2)
       Hinv=matinv2(H)
       func=100.
       niter=0
       do while ((func>1D-3).and.(niter.lt.1000))
          sb=(Hinv(1,1)*gb+Hinv(1,2)*ga)
          sa=(Hinv(2,1)*gb+Hinv(2,2)*ga)
          niter=niter+1
          sigma=0.1
          if (abs(sigma*sb).gt.stepmax) then
             sigma=abs(stepmax/sb)
          endif
          F=F-sigma*sb
          a=a-sigma*sa
          gb=float(Nstar)
          ga=float(Nstar+1)*float(Nstar)/2.
          H(1,1)=0. 
          H(1,2)=0. 
          H(2,2)=0. 
          do j=1,Nstar
             jf=float(j)
             tt=vi(j)*exp(F+a*jf)
             gb=gb-tt
             ga=ga-jf*tt
             H(1,1)=H(1,1)-tt
             H(1,2)=H(1,2)-jf*tt
             H(2,2)=H(2,2)-jf*jf*tt
          enddo
          H(2,1)=H(1,2)
          Hinv=matinv2(H)
          if ((abs(a).le.tiny(a)).or.(abs(F).le.tiny(F))) then
             func=max(gb,ga)
          else
             func=max(abs(gb/F),abs(ga/a))
          endif
       enddo
       H(:,:)=-H(:,:)
       Hinv=matinv2(H)
       a_err=sqrt(Hinv(2,2))
    endif !.not.kNN
    ! ### I should probably move this outside of the function...
    free_energy=F
    if (ISNAN(free_energy)) then
       write (*,*) "Density NaN at point",i 
       return
    endif
    if ((free_energy.gt.huge(free_energy)).or.(free_energy.lt.-huge(free_energy))) then
       write (*,*) "Density at point", free_energy
       return
    endif
    !Rho_err=dsqrt(dfloat(4*Nstar+2)/dfloat(Nstar*(Nstar-1))) ! I can do this outside!
    return

  end function free_energy

  pure function Hessian_L(a,F,Nstar,vi) result(H)
    ! Computes the Hessian of the log-likelihood
    ! L(F,A|{v_i,l}_l<=Nstar)
    ! as in Eq. 4 of Rodriguez et al. JCTC 2018
    real*8 :: H(2,2) !the Hessian
    real*8, intent(in) :: a,F
    integer :: j
    integer, intent(in) ::Nstar
    real*8, intent(in):: vi(Nstar)
    real*8 :: jf,tt
       
    H(1,1)=0. !gbb
    H(1,2)=0. !gab
    H(2,2)=0. !gaa
    do j=1,Nstar
       jf=float(j)
       tt=vi(j)*exp(F+a*jf)
       H(1,1)=H(1,1)-tt
       H(1,2)=H(1,2)-jf*tt
       H(2,2)=H(2,2)-jf*jf*tt
    enddo
    H(2,1)=H(1,2)
  end function Hessian_L

  pure function matinv2(A) result(B)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    real*8, intent(in) :: A(2,2)   !! Matrix
    real*8             :: B(2,2)   !! Inverse matrix
    real*8             :: detinv
    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
  end function matinv2

  real*8 function prefactor(dim)
    ! get prefactor for Volume calculation
    integer, intent(in) :: dim
    integer :: k,m,i,ms,ns
    real*8, parameter :: pi=3.14159265359
    if (mod(dim,2).eq.0) then
       k=dim/2
       m=1
       do i=1,k
          m=m*i
       enddo
       prefactor=pi**k/(dfloat(m))
    else
       k=(dim-1)/2
       ms=0.
       do i=1,k
          ms=ms+dlog(dfloat(i))
       enddo
       ns=ms
       do i=k+1,dim
          ns=ns+dlog(dfloat(i))
       enddo
       prefactor=2.*dexp(ms-ns+k*dlog(4*pi))
    endif

    return
  end function prefactor
  
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
