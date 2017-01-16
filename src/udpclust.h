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

//TODO: document
double critV[] = {2.91993, 2.54147, 2.12565, 1.96186, 1.77684, 1.64803, 1.57533, 1.5067, 1.41495, 1.34643, 1.27684,
                  1.22541, 1.18433, 1.1471, 1.09867, 1.06576, 1.03915, 1.01645, 0.952258, 0.905418, 0.898981, 0.875249,
                  0.87206, 0.836936, 0.819858, 0.800768, 0.797562, 0.781447, 0.766688, 0.738993, 0.732752, 0.728099,
                  0.70811, 0.684839, 0.684345, 0.679015, 0.67268, 0.653041, 0.647386, 0.627234, 0.624929, 0.611507,
                  0.61933, 0.611234, 0.608107, 0.587233, 0.577294, 0.577192, 0.560923, 0.560888, 0.546671, 0.548981,
                  0.547454, 0.536353, 0.525234, 0.518936, 0.522186, 0.512919, 0.509488, 0.502661, 0.495502, 0.494555,
                  0.489118, 0.485334, 0.476344, 0.479741, 0.465942, 0.459014, 0.464955, 0.461053, 0.453414, 0.454967,
                  0.446396, 0.444694, 0.437131, 0.443723, 0.435582, 0.433105, 0.428171, 0.423817, 0.421205, 0.427278,
                  0.423665, 0.421346, 0.412239, 0.410545, 0.410412, 0.411643, 0.410552, 0.401932, 0.400034, 0.394407,
                  0.389816, 0.391372, 0.386547, 0.378603, 0.380782, 0.373165, 0.378837, 0.374406, 0.372651, 0.371155,
                  0.366195, 0.363802, 0.363094, 0.360809, 0.356549, 0.351748, 0.352785, 0.351686, 0.353992, 0.351931,
                  0.348852, 0.348677, 0.344735, 0.34412, 0.338458, 0.337939, 0.337716, 0.333917, 0.331751, 0.326056,
                  0.326999};

class UDPClustering {
public:
    /**
     *
     * @param dist_mat distance matrix
     * @param max_knn  maximum nearest neighbours
     */
    UDPClustering(double *dist_mat,
                  double *Rho,
                  int *Cluster,
                  bool *filter,
                  int dimint,
                  int *Nlist,
                  size_t Nele,
            // int* id_err,
                  double sensibility,
                  int max_knn,
            // optional, because we can call get_densities to alloc this itself.
                  double *Rho_error = nullptr,
                  int *Nstar = nullptr
    ) :
            dist_mat(dist_mat, Nele, max_knn),
            Rho(Rho, Nele),
            Cluster(Cluster, Nele),
            Nele(Nele),
            dimint(dimint),
            min_knn(8), // TODO: this could also be a parameter
            max_knn(max_knn),
            sensibility(sensibility),
            Nlist(Nlist, Nele, max_knn),
            Nstar(Nstar, Nele),
            Rho_err(Rho_error, Nele),
            survivors(VecInt(Nele)) {
    }

    /// calculate densities of of points
    void get_densities();

    /// clustering routine
    void clustering();

    void merging();

protected:
    size_t dimint;
    size_t min_knn;
    size_t max_knn;

    VecDouble Rho; // density for all points
    VecDouble Rho_err;

    VecInt Cluster; // cluster ids

    VecBool filter;
    VecInt survivors;

    VecDouble2d dist_mat;
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
