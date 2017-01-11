//
// Created by marscher on 1/11/17.
//

#ifndef UDPCLUST_UDPCLUST_H
#define UDPCLUST_UDPCLUST_H

#include <vector>

class UDPClustering {
public:
    UDPClustering(std::vector<double>& dist_mat, int max_knn) :
            Nele(dist_mat.size()) {

    }

protected:
    int max_knn;
    std::vector<bool> filter;
    std::vector<bool> survivors;
    std::vector<double> dist_mat;
    size_t Nele;

    void clustering();

    // mark survivors, based on filter vector
    int get_survivors() {
        int Nsurv = 0;
        for (size_t i =0; i < survivors.size(); ++i) {
            if (! filter[i]) {
                Nsurv += 1;
                survivors[Nsurv] = i;
            }
        }
        return Nsurv;
    }

    // This function allows to get the distances in a matrix like style
    double gDist(size_t i, size_t j) {
        int i, j, k, l, m;
        l = max(i, j);
        m = min(i, j);
        k = (m - 1) * Nele - (m * m + m) / 2 + l;
        return dist_mat[k];
    }


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
void dp_advance(double* dist_mat, double* Rho, std::vector<bool>& filter,
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
    get_densities();

    instance.clustering();
}

//subroutine get_densities(id_err,dist_mat,Nele,dimint,Rho,Rho_err,filter,Nlist,Nstar)
void get_densities();


//   SUBROUTINE HPSORT(N,RA,s_order)
void heap_sort();

#endif //UDPCLUST_UDPCLUST_H
