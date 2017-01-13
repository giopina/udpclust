//
// Created by marscher on 1/11/17.
//

#ifndef UDPCLUST_UDPCLUST_H
#define UDPCLUST_UDPCLUST_H

#include "types.h"
namespace udpclust {


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
    ONLY_ONE_CLUSTER_AFTER_MERGING=10,
    ERROR_IN_ASSIGNATION=11,
    IG_UNDEFINED=12
};

class UDPClustering {
public:
    /**
     *
     * @param dist_mat distance matrix
     * @param max_knn  maximum nearest neighbours
     */
    UDPClustering(double* dist_mat, size_t Nele, size_t dimint, double sensibility, size_t max_knn) :
            dist_mat(dist_mat, Nele, max_knn),
            Nele(Nele),
            dimint(dimint),
            min_knn(8),
            max_knn(max_knn),
            sensibility(sensibility),
            Nlist(Nele, max_knn),
            Nstar(Nele)
    {
    }

protected:
    size_t dimint;
    size_t max_knn;
    size_t min_knn;

    VecDouble Rho; // density for all points
    VecDouble Rho_err;

    VecInt Cluster; // cluster ids

    VecBool filter;
    VecInt survivors;

    VecDouble2d dist_mat;
    double sensibility;

    size_t Nele;
    VecInt2d Nlist;   // Neighbour list within dc
    VecInt Nstar;   // N. of NN taken for comp dens

    VecDouble cent; //  Center Density
    VecDouble cent_err; // Center error

    VecDouble2d Bord; // border
    VecDouble2d Bord_err;

    VecInt Cluster_m; //  Cluster ID for the element
    int Nclus_m; //  Number of Cluster merged


    /// calculate densities of of points
    void get_densities();
    /// clustering routine
    void clustering();

    void merging(int Nclus);

    // mark survivors, based on filter vector
    int get_survivors();

    void find_border_densities();
};

/*############################################
!### MAIN CLUSTERING ALGORITHM SUBROUTINE ###
!############################################*/
int dp_advance(int*, double* dist_mat, int dimint, int Nele, double* rho, double* rho_err, char* filter, int* Nlist, float sensibility, int maxknn) {
    try {
        // create class
        auto instance = UDPClustering(dist_mat, );
        instance.get_densities();
        instance.clustering(); // perform clustering
        if (sensibility > 0.0) {
            return instance.merging(); //  Generate a new matrix without peaks within border error

        }
    } catch (int val) {
        return val;
    }
}

// overload with no args, returns the instance.
UDPClustering dp_advance() {

}
}

#endif //UDPCLUST_UDPCLUST_H
