//
// Created by marscher on 1/11/17.
//

#ifndef UDPCLUST_UDPCLUST_H
#define UDPCLUST_UDPCLUST_H

#include "types.h"

enum errors {
   /* if self.id_err==1 :   print "Select case error for distance type" #obsolete?
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
    else : print "unrecognized error"*/
    ONLY_ONE_CLUSTER=9,
    ONLY_ONE_CLUSTER_AFTER_MERGING=9,
    ERROR_IN_ASSIGNATION=11,
    IG_UNDEFINED=12
};

class UDPClustering {
public:
    UDPClustering(VecDouble &dist_mat, int max_knn) :
            Nele(dist_mat.size()),
            max_knn(max_knn),
            dist_mat(dist_mat),
            Nlist(Nele),
            Nstar(Nele)
    {

    }

protected:
    int max_knn;
    VecBool filter;
    VecBool survivors;
    VecDouble dist_mat;
    size_t Nele;
    VecInt Nlist;   // Neighbour list within dc
    VecInt Nstar;   // N. of NN taken for comp dens

    /**
     *     real*8,allocatable :: cent(:)       ! Center Density
    real*8,allocatable :: cent_err(:)   ! Center Error
    ! Underscore m implies data after automatic mergin
    integer,allocatable :: Cluster_m(:) ! Cluster ID for the element
    integer :: Nclus_m                  ! Number of Cluster merged
     */
    VecDouble cent; //  Center Density
    VecDouble cent_err; // Center error

    VecInt Cluster_m; //  Cluster ID for the element
    int Nclus_m; //  Number of Cluster merged

    // main routine
    void clustering();

    // mark survivors, based on filter vector
    int get_survivors();

    double gDist(size_t i, size_t j);

    void find_border_densities();
};

/*############################################
!### MAIN CLUSTERING ALGORITHM SUBROUTINE ###
!############################################*/
/** Arguments:
 * dist_mat,Cluster,Rho,filter,dimint,Nele,ND,id_err,sensibility

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

    integer :: survivors(Nele) ! ###
    integer :: Nsurv           ! ###

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
 */
// dp_advance has some subroutines, so is dp_advance actually a class in c++?
void dp_advance(double* dist_mat, double* Rho, VecBool filter,
                int Nele, int ND, int dimint, float sensibility, int maxknn) {
    // create class
    instance = UDPClustering();
    /*

    call get_densities(id_err,dist_mat,Nele,dimint,Rho,Rho_err,filter,Nlist,Nstar) ! ### my version
    call clustering(id_err)                      ! get Clusters
    if(sensibility.gt.0.0) then
       call merging(id_err) ! Generate a new matrix without peaks within border error
       Cluster=Cluster_m
    endif
    !    stop
    return
     */
    dens = get_densities();

    instance.clustering(dens);
}

//subroutine get_densities(id_err,dist_mat,Nele,dimint,Rho,Rho_err,filter,Nlist,Nstar)
void get_densities();

//TODO: what is s_order?
//   SUBROUTINE HPSORT(N,RA,s_order)
void heap_sort();

#endif //UDPCLUST_UDPCLUST_H
