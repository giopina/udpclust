//
// Created by marscher on 1/11/17.
//

#include <limits>
#include <cmath>

#include "udpclust.h"
#include "heap_sort.h"

// mark survivors, based on filter vector
int UDPClustering::get_survivors() {
    int Nsurv = 0;
    for (size_t i = 0; i < survivors.size(); ++i) {
        if (!filter[i]) {
            Nsurv += 1;
            survivors[Nsurv] = i;
        }
    }
    return Nsurv;
}

void UDPClustering::clustering() {
    int i, j, k;
    size_t ii, jj, kk;
    size_t ig;
    size_t l;
    bool idmax;
    VecDouble Rho_prob(Nele);  // Probability of having maximum Rho
    VecDouble Rho_copy;
    VecInt iRho, ordRho;
    VecInt Centers;
    double d, dmin;
    bool extend;

    int Nclus = 0;
    // Here I compute the probability of having density rho, g_i
    int Nsurv = get_survivors();
    //TODO: omp
    for (ii = 0; ii < Nsurv; ++ii) {
        i = survivors[ii];
        for (jj = 0; jj < Nsurv; ++jj) {
            j = survivors[jj];
            Rho_prob[i] = Rho_prob[i] - std::log(1.0 + std::exp(2. * (Rho[j] - Rho[i]))
                                                       / std::sqrt(Rho_err[i] * Rho_err[i] + Rho_err[j] * Rho_err[j]));
        }
    }

    //TODO: init iRho!
    Rho_copy = Rho; // copy
    heap_sort::sort(Rho_copy.data(), Rho_copy.size());

    for (i; i < Rho_copy.size(); ++i) {
        Rho_copy[i] -= Rho_prob[i];
    }
    //  ordRho is the complementary of iRho. Given an element, ordRho returns its order in density
    for (i; i < Nele; ++i) {
        ordRho[iRho[i]] = i;
    }

    /// Now I'm getting the clusters
    for (ii = 0; ii < Nsurv; ++ii) {
        i = survivors[ii];
        idmax = true;
        j = 1;
        while (idmax && j <= Nstar[i]) {
            if ((ordRho[i] > ordRho[Nlist(i, j)]) && (!filter[Nlist(i, j)])) {
                idmax = false;
                j += 1;
            }
        }
        if (idmax) {
            Nclus += 1;
            cluster[i] = Nclus;
        }
    }

    for (i = 0; i < Nele; ++i) {
        if (Cluster[i] != 0) {
            Centers[Cluster[i]] = i;
        }
    }

    /// assign not filtered
    /// TODO: change it with survivors
    if (Nclus <= 1) {
        // TODO: assign to -1;
        Cluster[:] = -1;
        throw ONLY_ONE_CLUSTER;
    }
    for (i = 0; i < Nele; ++i) {
        ig = -1;
        j = iRho[i];
        if (!filter[j] && Cluster[j] == 0) {
            dmin = std::numeric_limits<double>::max();
            for (k = 1; k < i - 1; ++k) {
                l = iRho[k];
                if (!filter[l] && gDist(j, l) <= dmin) {
                    ig = l;
                    dmin = gDist(j, l);
                }
            }
            if (ig == -1) {
                throw IG_UNDEFINED;
            } else {
                Cluster[j] = Cluster[ig];
            }
        }
    }


    /// find border densities
    VecDouble2d Bord(Nclus, Nclus));
    VecDouble2d Bord_err(Nclus, Nclus);
    VecInt2d eb(Nclus, Nclus);

    dmin = std::numeric_limits<double>::max();
    for (i = 0; i < Nele; ++i) {
        ig = -1;
        if (filter[i]) continue;
        for (j = 0; j < Nstar[i]; ++j) {
            l = Nlist[i, j];
            if (filter[l]) continue;
            if (cluster[l] == cluster[i]) continue;

            d = gDist(i, l);
            if (d < dmin) {
                dmin = d;
                ig = l;
            }
        }

        if (dmin > std::numeric_limits<double>::max() - 1) continue;
        extend = true;
        if (ig == -1) {
            throw IG_UNDEFINED;
        }

        for (k = 0; k < Nstar[i]; ++k) {
            if (filter[Nlist[i, k]]) continue;
            if (cluster[Nlist[i, k]] == cluster[i]) {
                extend = false;
                break;
            }
        }

        if (extend && Rho_prob[i] > Bord[cluster[i], cluster[ig]]) {
            Bord[cluster[i], cluster[ig]] = Rho_prob[i];
            Bord[cluster[ig], cluster[i]] = Rho_prob[i];

            Bord_err[cluster[i], cluster[ig]] = Rho_err[i];
            Bord_err[cluster[ig], cluster[i]] = Rho_err[i];
        }
    } //  enddo ! i=1,Nele

    for (int i = 0; i < Nclus - 1; ++i) {
        for (int j = i + 1; j < Nclus; ++j) {
            if (eb[i, j] != 0) {
                Bord[i, j] = Bord[j, i] = Rho[eb[i, j]];
            } else {
                Bord[i, j] = Bord[j, i] = 0;
            }
        }
    }

    /// Info per graph pre automatic merging
    cent.reserve(Nclus);
    cent_err.reserve(Nclus);

    for (i = 0; i < Nclus; ++i) {
        cent[i] = Rho[Centers[i]];
        cent_err[i] = Rho_err[Centers[i]];
        /// Modify centers in such a way that get more survival
        for (j = 0; j < Nele; ++j) {
            if (!filter[j] && (cluster[j] == i) && (Rho[j] - Rho_err[j] > (cent[i] - cent_err[i]))) {
                cent[i] = Rho[j];
                cent_err[i] = Rho_err[j];
            }
        }
    }

}

void UDPClustering::get_densities(int dimint) {
    /*
     * *
 *   subroutine get_densities(id_err,dist_mat,Nele,dimint,Rho,Rho_err,filter,Nlist,Nstar)
*/
    int i, j, k, m, n;
    double is, js, ks, ms, ns;

    const double pi = std::acos(-1);
    double prefactor;
    VecDouble Vols(max_knn);
    VecInt iVols(max_knn);
    bool viol = false;

    size_t limit = std::min(max_knn, (int) 0.5 * Nele);

    if (limit % 4 != 0) {
        limit += 4 - (limit % 4);
    }
    /// get prefactor for Volume calculation
    if (dimint % 2) {
        k = dimint / 2;
        m = 1;
        for (i = 0; i < k; ++i) {
            m *= i;
        }
        prefactor = std::pow(pi, k) / (double) m;
    } else {
        k = (dimint - 1) / 2;
        ms = 0;
        for (i = 0; i < k; ++i) {
            ms += std::log(i);
        }
        ns = ms;
        for (i = k + 1; i < dimint; ++i) {
            ns += std::log(i);
        }
        prefactor = 2*std::exp(ms -ns + k*std::log(4*pi));
    }

    // TODO: what for?
    double dimreal = (float) dimint;

    for (i=0; i < Nele; ++i) {
        Vols.clear();
        iVols.clear();
        for (j=0; j< max_knn; ++j) { // TODO:skip the first neighbour!
            Vols[j] = prefactor * std::pow(dist_mat(i, j), dimint);
        }
        viol=false;
        k = min_knn;
        n = 1;
        while(! viol) {
            rhg = k / Vols[k];
            dL = std::abs(4* (rhg*(Vols[k]-Vols[3*k/4])/float(k)));
            if (dL > critV[n]) {
                viol = true;
                break;
            }
            n+=1;
            k+=4;
            if (k > limit){
                viol=true;
            }
        }
        Nstar[i] = k-4;
        if (Nstar[i] < min_knn) {
            Nstar[i] = min_knn;
        }
    }
}