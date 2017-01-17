//
// Created by marscher on 1/11/17.
//

#ifndef UDPCLUST_UDPCLUST_H
#define UDPCLUST_UDPCLUST_H

#include "types.h"

namespace udpclust {


enum errors {
    ONLY_ONE_CLUSTER = 9,
    ONLY_ONE_CLUSTER_AFTER_MERGING = 10,
    ERROR_IN_ASSIGNATION = 11,
    IG_UNDEFINED = 12
};


class UDPClustering {
public:

    UDPClustering(double *dist_mat,
                  double *Rho,
                  int *Cluster,
                  bool *filter,
                  int dimint,
                  int *Nlist,
                  size_t Nele,
                  double sensibility,
                  int max_knn,
                  double *Rho_error,
                  int *Nstar
    );

    /// calculate densities of of points
    void get_densities();

    /// clustering routine
    void clustering();

    void merging();

protected:
    size_t dimint;
    size_t min_knn;
    size_t max_knn;

    VecDoubleExt Rho; // density for all points
    VecDoubleExt Rho_err;

    VecIntExt Cluster; // cluster ids

    VecBoolExt filter;
    VecInt survivors; // internal

    VecDoubleExt2d dist_mat;
    double sensibility;

    size_t Nele; // number of input points
    size_t Nclus; // number of clusters found; set by the clustering method.

    VecInt2d Nlist;   // Neighbour list within dc
    VecInt Nstar;   // N. of NN taken for comp dens

    VecDouble cent; //  Center Density
    VecDouble cent_err; // Center error

    VecDouble2d Bord; // border
    VecDouble2d Bord_err;

    VecInt Cluster_m; //  Cluster ID for the element
    int Nclus_m; //  Number of Cluster merged

    // mark survivors, based on filter vector
    int get_survivors();
};

/*############################################
!### MAIN CLUSTERING ALGORITHM SUBROUTINE ###
!############################################*/
void dp_advance(double* dist_mat, int* Cluster, double* Rho, bool* filter, int dimint,
                int* Nlist, int Nele, int* id_err, double sensibility, int maxknn);

/*
//TODO: needed?
// overload with no args, returns the instance.
UDPClustering dp_advance();

int clustering(int *);

int merging(int *);
 */

void get_densities(int* id_err, double* dist_mat, int Nele, int dimint, double* Rho, double* Rho_err,
                   bool* filter, int* Nlist, int* Nstar, int maxknn);

} // end of namespace udpclust

#endif //UDPCLUST_UDPCLUST_H
