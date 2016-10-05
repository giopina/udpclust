
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
!


module dp_clustering
  implicit none

contains  

  !############################################
  !### MAIN CLUSTERING ALGORITHM SUBROUTINE ###
  !############################################
  subroutine dp_advance(dist_mat,Cluster,Rho,filter,dimint,Nele,ND,id_err,merge)
!   implicit none
    !#####################
    integer,intent(in) :: merge
    integer,intent(inout) :: id_err
    !!Global variables
    real*8,intent(in) :: dist_mat(ND)       ! Distance matrix !###
    integer,intent(in) :: Nele                   ! Number of elements
    integer,intent(in) :: ND                     ! Number of distances

    integer,intent(in) :: dimint                 ! integer of dimset (avoid real*8 calc)
    integer,parameter :: maxknn=496   ! maximum number of neighbours to explore
    integer,parameter :: minknn=8     ! minimum number of neighbours to explore
    !integer,parameter :: critdim=123  ! dimension of the array with the critical values
    integer :: limit
    real*8, parameter  :: Sharpness=1.0 ! Filter criteria, the lower, the smoother
    real*8,intent(inout) :: Rho(Nele)        ! Density
    real*8,allocatable :: Rho_err(:)    ! Density error
    real*8,allocatable :: Rho_prob(:)   ! Probability of having maximum Rho
    real*8,allocatable :: dc(:)         ! Distance for density calculation
    logical,intent(inout) :: filter(Nele)  ! Point with anomalous density
    logical,allocatable :: backgr(:)  ! Point assigned to background noise
    logical,allocatable :: scarti(:)  ! Points filtered or assigned to background noise
    integer,allocatable :: Nlist(:,:) ! Neighbour list within dc
    integer,allocatable :: Nstar(:)   ! Number of neighbours taken for computing density 
   integer,allocatable :: Centers(:) ! Centers of the peaks
    integer,intent(inout) :: Cluster(Nele) ! Cluster ID for the element
    integer :: Nclus                  ! Number of Cluster
    integer,allocatable :: candidates_B(:)     ! Cluster of the Border candidates
    real*8,allocatable :: Bord(:,:)     ! Border Densities
    real*8,allocatable :: Bord_err(:,:) ! Border Densities Error
    integer,allocatable :: eb(:,:)    ! Border elements
    integer,allocatable :: Pop(:)     ! Cluster Population
    real*8,allocatable :: cent(:)       ! Center Density
    real*8,allocatable :: cent_err(:)   ! Center Error
    ! Underscore m implies data after automatic merging
    integer :: Nclus_m                  ! Number of Cluster merged
    real*8,allocatable :: Bord_m(:,:)     ! Border Densities
    real*8,allocatable :: Bord_err_m(:,:) ! Border Densities Error
    integer,allocatable :: Pop_m(:)     ! Cluster Population
    real*8,allocatable :: cent_m(:)       ! Center Density
    real*8,allocatable :: cent_err_m(:)   ! Center Error
    integer,allocatable :: Centers_m(:) ! Centers of the peaks
    integer,allocatable :: Cluster_m(:) ! Cluster ID for the element
    ! 
    integer:: iglob,n,ai,bi
    !####################################################

    id_err=0

    call get_densities_and_dc(id_err)            ! from the distances, compute densities and dc
!    open (22,file="my_densities.dat")
!    write(22,*) id_err
    call clustering(id_err)                      ! get Clusters
!    write(22,*) id_err
!    close(22)
    if (merge.eq.1) then
       call merging(id_err) ! Generate a new matrix without peaks within border error  
       Cluster=Cluster_m
    endif
    !    stop
    return

  contains

    subroutine get_densities_and_dc(id_err)
      use critfile
      implicit none
      integer :: id_err
      !! Local variables
      integer :: i,j,k,m,n
      integer :: i1,i2 ! ### my integer for loops
      integer :: kadd
      integer :: niter,nfilter
      real*8,allocatable :: Vols(:)
      integer,allocatable :: iVols(:)
      !real*8 :: critV(critdim)
      real*8, parameter :: pi=3.14159265359
      real*8 :: prefactor
!      real*8 :: rhg,rh1,rh2,rh3,rh4,dL
      real*8 :: rhg,dl
      real*8,dimension(4) :: rh
!      real*8 :: rjk1,rjk2,rjk3,rjk4,rjaver,rjfit
      real*8 :: rjaver,rjfit
      real*8,dimension(4) :: rjk
      logical :: viol
      real*8 :: xmean,ymean,a,b,c                                     !  FIT
!      real*8 :: x1,x2,x3,x4                                           ! STUFF LINEAR NOT EQUALLY POPULATED PARTITIONS
      real*8, dimension(4) :: x
      integer :: partit(4)
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
         m=1
         do i=1,k
            m=m*i
         enddo
         n=m
         do i=k+1,dimint
            n=n*i
         enddo
         prefactor=2.*dfloat(m)*(4.*pi)**k/(dfloat(n))
      endif

      allocate (Nstar(Nele))
      allocate (Rho_err(Nele))
      allocate (Nlist(Nele,limit))
      allocate (Vols(Nele))
      allocate (iVols(Nele))

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
         ! ### (ma ha senso sta roba?)
         kadd=1
         if (Nstar(i).eq.limit) kadd=0
         do while (kadd.eq.1)
            if ((abs(gDist(i,Nlist(i,Nstar(i)))-gDist(i,Nlist(i,Nstar(i)+1)))).lt.9.99D-99) then
               Nstar(i)=Nstar(i)+1
               if (Nstar(i).eq.limit) kadd=0
            else
               kadd=0
            endif
         enddo
!         partit(:)=0 ! ### che e'?!
!         do j=1,Nstar(i)
!            k=mod(j,4)
!            if (k.eq.0) k=4 ! ### !
!            partit(k)=partit(k)+1
!         enddo
         ! ### direi che il ciclo di qui sopra si puo' fare in due righe ovvero
         partit(:)=Nstar(i)/4 !-MOD(Nstar(i),4)
         partit(:MOD(Nstar(i),4))=partit(:MOD(Nstar(i),4))+1
         !
         ! ### ho capito il senso. Dovrebbe essere per ovviare al fatto che Nstar non e' per forza un multiplo di 4.
         ! ### non sono sicuro se questo sia il modo migliore per farlo...
         x(1)=dfloat(partit(1))/dfloat(Nstar(i))*0.5
         x(2)=dfloat(partit(1))/dfloat(Nstar(i))+dfloat(partit(2))/dfloat(Nstar(i))*0.5
         x(3)=dfloat(partit(1)+partit(2))/dfloat(Nstar(i))+dfloat(partit(3))/dfloat(Nstar(i))*0.5
         x(4)=dfloat(partit(1)+partit(2)+partit(3))/dfloat(Nstar(i))+dfloat(partit(4))/dfloat(Nstar(i))*0.5
         x=x*4.
         rhg=dfloat(Nstar(i))/Vols(Nstar(i)) ! Rho without fit
         if (Nstar(i).le.0) then
            Rho(i)=rhg
            Rho_err(i)=Rho(i)/dsqrt(dfloat(Nstar(i)))
         else
            ! get inv rho of the four quarters
!            j=Nstar(i)/4
!            a=dfloat(j)
            rh(1)=Vols(partit(1))/dfloat(partit(1))
            rh(2)=(Vols(partit(1)+partit(2))-Vols(partit(1)))/dfloat(partit(2))
            rh(3)=(Vols(partit(1)+partit(2)+partit(3))-Vols(partit(1)+partit(2)))/dfloat(partit(3))
            rh(4)=(Vols(partit(1)+partit(2)+partit(3)+partit(4))-Vols(partit(1)+partit(2)+partit(3)))/dfloat(partit(4))
            ! make the quadratic fit rhj=1/rho+C*j^2
            xmean=0.25*sum(x)
            ymean=0.25*sum(rh)
            b=SUM((x-xmean)**2)
            a=SUM((x-xmean)*(rh-ymean))
            a=a/b
            rjfit=ymean-a*xmean
            ! Perform jacknife resampling for estimate the error (it includes statistical
            ! error and curvature error) 
            ! 4 out
            do i1=1,4
               xmean=(SUM(x)-x(i1))/dfloat(3)
               ymean=(SUM(rh)-rh(i1))/dfloat(3)
               b=SUM((x-xmean)**2)-(x(i1)-xmean)**2
               a=SUM((x-xmean)*(rh-ymean))-(x(i1)-xmean)*(rh(i1)-ymean)
               a=a/b
               !###
               rjk(i1)=ymean-a*xmean
            enddo
            rjaver=0.25*SUM(rjk)
            Rho(i)=4.*rjfit-3.*rjaver
!            Rho_err(i)=0.75*((rjk1-rjaver)**2+(rjk2-rjaver)**2+(rjk3-rjaver)**2+(rjk4-rjaver)**2)
            Rho_err(i)=0.75*SUM((rjk-rjaver)**2)
            Rho(i)=1./Rho(i)
            Rho_err(i)=max(Rho(i)/sqrt(float(Nstar(i))),Rho(i)*Rho(i)*sqrt(Rho_err(i)))
         endif
         ! density output: N, n*,density,Error in density,density from knn formula
!         write (22,'(i6,1x,i3,1x,3(f28.9,1x))') i,Nstar(i),Rho(i),Rho_err(i),rhg 
         ! ### ### COS'E' rhg????????? ### ### tipo il raggio...
      enddo
      deallocate (Vols,iVols)
!      close (22)
      !      allocate (filter(Nele))
      allocate (backgr(Nele))
      allocate (scarti(Nele))

      backgr(:)=.false.
      scarti(:)=.false.
      ! Filter with neighbours density (Iterative version)
      filter(:)=.false.
      viol=.true.
      niter=0
      do while (viol)
         niter=niter+1
         viol=.false.
         nfilter=0
         do i=1,Nele
            ! compute average density in the neighborhood and standard dev.
            if (.not.filter(i)) then
               if ((Rho(i).lt.0).or.(Rho_err(i).gt.Rho(i))) then
                  filter(i)=.true.
                  viol=.true.
               else
                  a=0.
                  b=0.
                  n=0
                  do j=1,Nstar(i) 
!                     if (i.eq.5060) write (6,'(i5,i3,i5,f20.5,l,f20.5)') i,j,&
!                          Nlist(i,j),Rho(Nlist(i,j)),filter(Nlist(i,j)),gDist(i,Nlist(i,j))
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
!      write (51,*) "Filter iterations=",niter
!      write (51,*) "Filtered elements=",nfilter
      scarti(:)=filter(:)
!      open (22,file='filter.dat')
!      do i=1,Nele
!         if (filter(i)) then
!            write (22,*) 1
!         else
!            write (22,*) 0
!         endif
!      enddo
!      close (22)
      !Get dc
      allocate (dc(Nele))
!      open (21,file="dc.dat")
      do i=1,Nele
         j=Nlist(i,Nstar(i))
         dc(i)=gDist(i,j)
!         write (21,*) dc(i)
      enddo
!      close (21)
      return
211   id_err=8
      return
    end subroutine get_densities_and_dc

    subroutine clustering(id_err)
      implicit none
      integer :: id_err
      !! Local variables
      integer :: i,j,k
      integer :: ig,jg
      integer :: nassign
      integer :: nonfilter
      integer :: iref
      integer :: l
      logical :: idmax
      real*8,allocatable :: Rho_copy(:)
      integer,allocatable :: iRho(:)
      integer,allocatable :: ordRho(:) 
      real*8 :: d,dmin
      logical :: extend
      !!
      id_err=0
      !      allocate (Cluster(Nele))
      allocate (candidates_B(Nele))
      !! Identify centers: delta>dc eqv. Max(rho) within dc
      Cluster(:)=0
      Nclus=0
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
      do i=1,Nele
         if (.not.scarti(i)) then
            idmax=.true.
            j=1
            do while (idmax .and. (j.le.Nstar(i)))
               if ((Rho_prob(i).le.Rho_prob(Nlist(i,j))).and.(.not.scarti(Nlist(i,j))))  idmax=.false.
               j=j+1
            enddo
            if (idmax) then
               Nclus=Nclus+1
               Cluster (i)=Nclus
            endif
         endif
      enddo
!      write (51,*) "Nclus:",Nclus
      allocate (Centers(Nclus))
      nonfilter=0
      do i=1,Nele
         if (Cluster(i).ne.0) then
            Centers(Cluster(i))=i
 !           write (51,*) "Cluster:",Cluster(i),"Center:",i
         endif
         if (.not.filter(i)) nonfilter=nonfilter+1
      enddo
      nassign=Nclus
!      write(22,*) id_err,Nclus !###
      if (Nclus.gt.1) then

         ! copy of rho (not efficient, but clarify the code)

         allocate (Rho_copy(Nele))
         allocate (iRho(Nele))
         Rho_copy(:)=-Rho_prob(:)

         call HPSORT(Nele,Rho_copy,iRho) ! iRho contains the order in density (iRho(1) is the element with highest Rho...)
         deallocate (Rho_copy)

         allocate (ordRho(Nele))
         do i=1,Nele
            ordRho(iRho(i))=i                 ! ordRho is the complementary of iRho. Given an element, ordRho returns its order in density
         enddo

         ! Assign not filtered

!         write (51,*) "Assignation of not filtered points"
         do i=1,Nele
            j=iRho(i)
            if (.not.filter(j).and.Cluster(j).eq.0) then
               dmin=9.9E29
               do k=1,i-1
                  l=iRho(k)
                  if (.not.filter(l)) then
                     if (gDist(j,l).le.dmin) then
                        ig=l
                        dmin=gDist(j,l)
                     endif
                  endif
               enddo
               Cluster(j)=Cluster(ig)
            endif
         enddo


         ! Assign filtered to the same Cluster as its nearest unfiltered neighbour (and
         ! count cluster population)
         allocate (Pop(Nclus))
         Pop(:)=0
         do i=1,Nele
            if (Cluster(i).eq.0) then
               dmin=9.9d99
               do j=1,Nele
                  if ((Cluster(j).ne.0).and.(.not.filter(j))) then
                     d=gDist(i,j)
                     if (d.le.dmin) then
                        dmin=d
                        jg=j
                     endif
                  endif
               enddo
               Cluster(i)=Cluster(jg)
            endif
            Pop(Cluster(i))=Pop(Cluster(i))+1
         enddo
         ! find border densities

         allocate (Bord(Nclus,Nclus),Bord_err(Nclus,Nclus),eb(Nclus,Nclus))
         Bord(:,:)=-9.9D99
         Bord_err(:,:)=0.
         eb(:,:)=0
         candidates_B(:)=0

         do i=1,Nele
            if (.not.scarti(i)) then
               dmin=9.9d99
               do j=1,Nele
                  if (.not.scarti(j)) then
                     if (cluster(j).ne.cluster(i)) then
                        d=gDist(i,j)
                        if (d.lt.dmin) then
                           dmin=d
                           jg=j
                        endif
                     endif
                  endif
               enddo
               if (dmin.le.dc(i)) then
                  iref=i
                  k=0
                  extend=.true.
                  do while ( (k.lt.Nstar(i)).and.extend)
                     k=k+1
                     if (cluster(Nlist(i,k)).eq.cluster(i)) then
                        if (gDist(Nlist(i,k),jg).lt.dmin) extend=.false.
                     endif
                  enddo
                  if (extend) then
                     candidates_B(i)=cluster(jg)
                     if (Rho_prob(iref).gt. Bord(cluster(i),cluster(jg))) then
                        Bord(cluster(i),cluster(jg))=Rho_prob(iref)
                        Bord(cluster(jg),cluster(i))=Rho_prob(iref)
                        Bord_err(cluster(i),cluster(jg))=Rho_err(iref)
                        Bord_err(cluster(jg),cluster(i))=Rho_err(iref)
                        eb(cluster(i),cluster(jg))=iref
                        eb(cluster(jg),cluster(i))=iref
                     endif
                  endif
               endif
               !   endif
            endif
         enddo
!         open (22,file="Borders_pre_merge")
         do i=1,Nclus-1
            do j=i+1,Nclus
               if (eb(i,j).ne.0) then
                  Bord(i,j)=Rho(eb(i,j))
                  Bord(j,i)=Rho(eb(i,j))
               else
                  Bord(i,j)=0.
                  Bord(j,i)=0.
               endif
!               write (22,*) i,j,Bord(i,j),Bord_err(i,j)
            enddo
         enddo
!         close (22)

         deallocate (Rho_prob)
         deallocate (iRho)
         deallocate (ordRho)
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
      integer :: i,j,k,n,l,alive,dead,niter
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
      Nbarr=(Nclus*Nclus-Nclus)/2
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
         do n=1,Nbarr
            k=iBarrier(n)
            i=Bcorr(k,1)
            j=Bcorr(k,2)
!            write (51,'(a,2i3,f16.4)') "BARR",i,j,Bord(i,j)
         enddo
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
                  c1=cent(i)-cent_err(i)
                  c2=cent(j)-cent_err(j)
                  b12=Bord(i,j)+Bord_err(i,j)
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
!                     write (51,'(i4,1x,a,i4,a,i4)') n,"MERGED:",dead,"->",alive
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
!         write (51,*) "MERGING ITER:",niter
      enddo

      ! get dictionary

      Nclus_m=0
      do i=1,Nclus
         if (Survive(i)) Nclus_m=Nclus_m+1
      enddo
!      write (51,*) "TOTAL CLUSTERS AFTER MERGING:",Nclus_m
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
         allocate (Pop_m(Nclus_m))
         allocate (Bord_m(Nclus_m,Nclus_m),Bord_err_m(Nclus_m,Nclus_m))
         allocate (cent_m(Nclus_m),cent_err_m(Nclus_m))
         allocate (Centers_m(Nclus_m))
         Pop_m(:)=0
!         open (21,file="cluster_merged")
         do i=1,Nele
            Cluster_m(i)=O2M(Cluster_m(i))
!            write (21,*) i,Cluster_m(i)
            Pop_m(Cluster_m(i))=Pop_m(Cluster_m(i))+1
         enddo
!         close (21)
         do i=1,Nclus_m
            do j=i+1,Nclus_m
               Bord_m(i,j)=Bord(M2O(i),M2O(j)) 
               Bord_err_m(i,j)=Bord_err(M2O(i),M2O(j)) 
               Bord_m(j,i)=Bord(M2O(i),M2O(j)) 
               Bord_err_m(j,i)=Bord_err(M2O(i),M2O(j)) 
            enddo
            cent_m(i)=cent(M2O(i))
            cent_err_m(i)=cent_err(M2O(i))
            Centers_m(i)=Centers(M2O(i))
         enddo

!         open (21,file="Properties_after_merge")
!         do i=1,Nclus_m
!            write (21,*) i,cent_m(i),cent_err_m(i),Centers_m(i)
!         enddo
!         do i=1,Nclus_m-1
!            do j=i+1,Nclus_m
!               write (21,*) i,j,Bord_m(i,j),Bord_err_m(i,j)
!            enddo
!         enddo
!         close (21)
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
10    continue
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
20    if(J.le.IR)then
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
  end subroutine dp_advance

end module dp_clustering
