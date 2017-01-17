//
// Created by marscher on 1/11/17.
//

#include <limits>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include "udpclust.h"

namespace {
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
}

namespace udpclust {

/*############################################
!### MAIN CLUSTERING ALGORITHM SUBROUTINE ###
!############################################*/
//subroutine dp_advance(dist_mat,Cluster,Rho,filter,dimint,Nlist,Nele,id_err,sensibility,maxknn)
void dp_advance(double *dist_mat, int *Cluster, double *Rho, bool *filter, int dimint,
                int *Nlist, int Nele, int *id_err, double sensibility, int maxknn) {
    int result = 0;
    try {
        // create class
        //TODO: Rho_err not handed in, so we create it here temporarily
        VecDouble Rho_err(Nele);
        VecInt Nstar(Nele);

        UDPClustering instance(dist_mat, Rho, Cluster, filter,
                               dimint, Nlist, Nele, sensibility,
                               maxknn, Rho_err.data(), Nstar.data());
        instance.get_densities();
        instance.clustering(); // perform clustering
        if (sensibility > 0.0) {
            instance.merging(); // Generate a new matrix without peaks within border error
        }
    } catch (int error_id) {
        result = error_id;
    }
    *(id_err) = result;
}

void get_densities(int *id_err, double *dist_mat, int Nele, int dimint, double *Rho, double *Rho_err,
                   bool *filter, int *Nlist, int *Nstar, int maxknn) {
    int result = 0;
    try {
        float sensibility = 0.0;
        printf("1***********\n");
        UDPClustering instance(dist_mat, Rho, nullptr, filter, dimint, Nlist, Nele, sensibility, maxknn, Rho_err,
                                      Nstar);
        printf("2***********\n");

        instance.get_densities(); // writes the pointers Rho, rho_err and filter
        printf("3***********\n");
    } catch (int error_id) {
        result = error_id;
    }
    *(id_err) = result;
}

UDPClustering::UDPClustering(
        double *dist_mat,
        double *Rho,
        int *Cluster,
        bool *filter,
        int dimint,
        int *Nlist,
        size_t Nele,
        double sensibility,
        int max_knn,
        // optional, because we can call get_densities to alloc this itself.
        double *Rho_error,
        int *Nstar
) :
        dist_mat(VecDoubleExt2d(dist_mat, Nele, max_knn),
        Rho(Rho, Nele),
        Cluster(Cluster, Nele),
        filter(filter, Nele),
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

// mark survivors, based on filter vector
int UDPClustering::get_survivors() {
    int Nsurv = 0;
    for (int i = 0; i < survivors.size(); ++i) {
        if (!filter[i]) {
            Nsurv += 1;
            survivors[Nsurv] = i;
        }
    }
    return Nsurv;
}

void UDPClustering::clustering() {
    size_t i, j, k, l;
    int ii, jj;
    int ig;
    bool idmax;
    VecDouble Rho_prob(Nele);  // Probability of having maximum Rho
    VecInt iRho, ordRho;
    VecInt Centers;
    double d, dmin;
    bool extend;

    Nclus = 0;
    // Here I compute the probability of having density rho, g_i
    int Nsurv = get_survivors();
    Rho_prob = 0;
    //TODO: omp
    for (ii = 0; ii < Nsurv; ++ii) {
        i = survivors[ii];
        for (jj = 0; jj < Nsurv; ++jj) {
            j = survivors[jj];
            Rho_prob[i] -= std::log(1.0 + std::exp(2. * (Rho[j] - Rho[i]))
                                          / std::sqrt(Rho_err[i] * Rho_err[i] + Rho_err[j] * Rho_err[j]));
        }
    }

    /// iRho contains the order in density (iRho(1) is the element with highest Rho...)
    iRho = VecInt(Rho.size());
    for (i = 0; i < iRho.size(); ++i) {
        iRho[i] = i;
    }
    // iRho is the array which sorts Rho descending
    std::sort(iRho.data(), iRho.data() + iRho.size(), [this](int a, int b) {
        return this->Rho[a] > this->Rho[b];
    });
    std::cout << "Rho[iRho[0]]=" << Rho[iRho[0]] << "Rho[iRho[iRho.size()-1]] = " << Rho[iRho[iRho.size() - 1]]
              << std::endl;
    assert(Rho[iRho[0]] > Rho[iRho[1]]);

    //  ordRho is the complementary of iRho. Given an element, ordRho returns its order in density
    for (i = 0; i < Nele; ++i) {
        ordRho[iRho[i]] = i;
    }

    /// Now I'm getting the clusters
    for (ii = 0; ii < Nsurv; ++ii) {
        i = survivors[ii];
        idmax = true;
        for (j = 0; (idmax && (j <= Nstar[i])); ++j) {
            if ((ordRho[i] > ordRho[Nlist(i, j)]) && (!filter[Nlist(i, j)])) {
                idmax = false;
            }
        }
        if (idmax) {
            Nclus += 1;
            Cluster[i] = Nclus;
        }
    }

    Centers = VecInt(Nele);
    for (i = 0; i < Nele; ++i) {
        if (Cluster[i] != 0) {
            Centers[Cluster[i]] = i;
        }
    }

    /// TODO: change it with survivors
    if (Nclus <= 1) {
        Cluster = 1;
        throw ONLY_ONE_CLUSTER;
    }

    /// assign not filtered
    for (i = 0; i < Nele; ++i) {
        i = iRho[j];
        if (!filter[j] && Cluster[j] == 0) {
            ig = -1;
            //dmin = std::numeric_limits<double>::max();
            for (k = 0; k < Nstar[i]; ++k) {
                l = Nlist(i, k);
                if (!filter[l] && Rho_prob[i] < Rho_prob[l]) {
                    ig = l;
                    //dmin = dist_mat(i, k);
                }
            }
            if (ig == -1) {
                throw IG_UNDEFINED;
            } else {
                Cluster[j] = Cluster[ig];
            }
        }
    }

    /*
     * Assign filtered to the same Cluster as its nearest unfiltered neighbour
     * what happens if all neighbors are filtered?
     * this can happen, at least in principle. What should I do then?
     * option 1: assign the point to the same cluster of the closest filtered point (you should do this in the end, but it will depend on the order you consider the points...
     */

    bool viol = false;
    for (i = 0; i < Nele; ++i) {
        ig = -1;
        if (!filter[i]) {
            throw IG_UNDEFINED;
        }
        dmin = std::numeric_limits<double>::max();
        for (k = 0; k < max_knn; ++k) {
            l = Nlist(i, k);
            if (Cluster[l] != 0 && !filter[l]) {
                d = dist_mat(i, k);
                if (d < dmin) {
                    dmin = d;
                    ig = l;
                }
            }
        }
        if (ig != -1) {
            Cluster[i] = Cluster[ig];
        } else {
            viol = true;
        }
    }

    int Nnoass = 0;
    while (viol) {
        viol = false;
        bool NEWASS = false;
        for (i = 0; i < Nele; ++i) {
            if (Cluster[i] == 0) {
                dmin = std::numeric_limits<double>::max();
                ig = -1;

                for (k = 0; k < max_knn; ++k) {
                    j = Nlist(i, k);
                    if (Cluster[j] != 0) {
                        d = dist_mat(i, k);
                        if (d < dmin) {
                            dmin = d;
                            ig = j;
                        }
                    }
                }
                if (ig != -1) {
                    Cluster[i] = Cluster[ig];
                    NEWASS = true;
                } else {
                    viol = true;
                    Nnoass += 1;
                }
            }
        }
        if (!NEWASS) {
            std::cerr << "ERROR: cannot assign " << Nnoass << " filtered points to clusters." << std::endl;
            throw IG_UNDEFINED;
        }
    }

    //// find border densities
    //VecDouble2d Bord(Nclus, Nclus);
    //VecDouble2d Bord_err(Nclus, Nclus);
    VecInt2d eb(Nclus, Nclus);
    eb = 0;

    for (i = 0; i < Nele; ++i) {
        ig = -1;
        dmin = std::numeric_limits<double>::max();
        //if (filter[i]) continue;
        for (j = 0; j < Nstar[i]; ++j) {
            l = Nlist(i, j);
            if (filter[l]) continue;
            if (Cluster[l] == Cluster[i]) continue;

            d = dist_mat(i, l);
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
            l = Nlist(i, k);
            if (filter[l]) continue;
            if (Cluster[l] == Cluster[i]) {
                extend = false;
                break;
            }
            d = std::numeric_limits<double>::max();
            for (jj = 0; jj < max_knn; ++jj) {
                j = Nlist(l, jj);
                if (j == ig) {
                    d = dist_mat(l, jj);
                    break;
                }
            }
            if (d < dmin) {
                extend = false;
                break;
            }
        }
        //TODO: why the heck the cluster variable is lowercase in fortran?
        if (extend && Rho_prob[i] > Bord(Cluster[i], Cluster[ig])) {
            Bord(Cluster[i], Cluster[ig]) = Rho_prob[i];
            Bord(Cluster[ig], Cluster[i]) = Rho_prob[i];

            Bord_err(Cluster[i], Cluster[ig]) = Rho_err[i];
            Bord_err(Cluster[ig], Cluster[i]) = Rho_err[i];

            eb(Cluster[i], Cluster[ig]) = i;
            eb(Cluster[ig], Cluster[i]) = i;
        }
    } //  end for loop over all elements.

    for (i = 0; i < Nclus - 1; ++i) {
        for (j = i + 1; j < Nclus; ++j) {
            if (eb(i, j) != 0) {
                Bord(i, j) = Bord(j, i) = Rho[eb(i, j)];
            } else {
                Bord(i, j) = Bord(j, i) = 0;
            }
        }
    }

    /// Info per graph pre automatic merging
    cent = VecDouble(Nclus);
    cent_err = VecDouble(Nclus);

    for (i = 0; i < Nclus; ++i) {
        cent[i] = Rho[Centers[i]];
        cent_err[i] = Rho_err[Centers[i]];
        /// Modify centers in such a way that get more survival
        for (j = 0; j < Nele; ++j) {
            if (!filter[j] && (Cluster[j] == i) && (Rho[j] - Rho_err[j] > (cent[i] - cent_err[i]))) {
                cent[i] = Rho[j];
                cent_err[i] = Rho_err[j];
            }
        }
    }

}


/**
 * calculates the densities for all points (Nele) and the errors on that (Rho and Rho_err).
 * Sets the filtered points (eg. outliers)
 * @param dimint
 */
void UDPClustering::get_densities() {
    int i, j, k, m, n;
    double /*is, js, ks,*/ ms, ns;

    const double pi = std::acos(-1);
    double prefactor;

    //double slope;
    double temp_rho, temp_err;

    bool viol = false;
    double xmean, ymean, a, b, c;          // FIT
    double xsum, ysum, x2sum, xysum;

    int Npart, savNstar, fin;
    double rhg, dL, rjaver, rjfit;

    VecDouble Vols(max_knn);
    VecDouble x, rh, rjk;
    VecInt iVols(max_knn);

    size_t limit = std::min(max_knn, (size_t) 0.5 * Nele);

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
        prefactor = 2 * std::exp(ms - ns + k * std::log(4 * pi));
    }

    for (i = 0; i < Nele; ++i) {
        iVols = 0; // reset all elements = -1
        for (j = 0; j < max_knn; ++j) { // TODO:skip the first neighbour!
            Vols[j] = prefactor * std::pow(dist_mat(i, j), dimint);
        }
        viol = false;
        k = min_knn;
        n = 1;
        while (!viol) {
            rhg = k / Vols[k];
            dL = std::abs(4 * (rhg * (Vols[k] - Vols[3 * k / 4]) / double(k)));
            if (dL > critV[n]) {
                viol = true;
                continue;
            }
            n += 1;
            k += 4;
            if (k > limit) {
                viol = true;
            }
        }
        Nstar[i] = k - 4;
        if (Nstar[i] < min_knn) {
            Nstar[i] = min_knn;
        }
        Rho_err[i] = std::numeric_limits<double>::max();
        rhg = (float) Nstar[i] / Vols[Nstar[i]]; // Rho without fit.
        Rho[i] = rhg;
        savNstar = Nstar[i];
        Npart = 4;
        fin = Nstar[i] / 2;
        while (Npart <= fin) {
            if (Nstar[i] % Npart < Nstar[i] / 4) {
                if (Nstar[i] % Npart != 0) {
                    Nstar[i] -= Nstar[i] % Npart;
                }
                /// get inv rho of partition
                // FIXME: double free corruption...
                VecDouble x(Npart);
                VecDouble rh(Npart);
                VecDouble rjk(Npart);
                j = Nstar[i] / Npart;
                for (k = 0; k < Npart; ++k) {
                    n = n + j;
                    x[k] = k;
                    if (k == 1) {
                        rh[k] = Vols[n] / a;
                    } else {
                        rh[k] = (Vols[n] - Vols[n - j]) / a;
                    }
                }
                xmean = x.sum() / (double) Npart;
                ymean = rh.sum() / (double) Npart;
                b = 0;
                c = 0;
                for (k = 0; k < Npart; ++k) {
                    a = x[k] - xmean;
                    b += a * (rh[k] - ymean);
                    c += a * a;
                }
                a = b / c;
                //slope = a;
                rjfit = ymean - a * xmean;
                //yintercept = rjfit;
                /// Perform jacknife resampling for estimate the error( it includes the statistical error and curvature error)
                xsum = x.sum();
                ysum = rh.sum();
                x2sum = x.pow(2).sum();
                xysum = (x * rh).sum();

                for (n = 0; n < Npart; ++n) {
                    xmean = (xsum - x[n]) / (double) (Npart - 1);
                    ymean = (ysum - rh[n]) / (double) (Npart - 1);
                    c = x2sum - x[n] * x[n];
                    c -= xmean * xmean * (Npart - 1);
                    b = xysum - x[n] * rh[n];
                    b -= xmean * xmean * (Npart - 1);
                    a = b / c;
                    rjk[n] = ymean - a * xmean;
                }

                rjaver = rjk.sum() / (double) Npart;
                temp_rho = Npart * rjfit - (Npart - 1) * rjaver;
                temp_err = 0.0;
                for (k = 0; k < Npart; ++k) {
                    temp_err += (rjk[k] - rjaver) * (rjk[k] - rjaver);
                }
                temp_err = (Npart - 1) * temp_err / Npart;
                temp_rho = 1. / temp_rho;
                temp_err = std::max(temp_rho / std::sqrt((double) Nstar[i]), temp_rho * temp_rho * std::sqrt(temp_err));
                if (temp_err > Rho_err[i]) {
                    Rho[i] = Rho_err[i] = temp_rho;
                    //partGood = Npart;
                }
            } // end: if(Nstar[i] % Npart < Nstar[i] /4 )
            Npart += 1;
            Nstar[i] = savNstar;
        } // while
    }

    /// Filter with neighbours density (Iterative version)
    filter = false; // all elements = false
    viol = true;
    for (i = 0; i < Nele; ++i) {
        if (Rho[i] < 0 || Rho_err[i] > Rho[i] || Rho[i] > 1E308) {
            filter[i] = true;
        }
    }

    while (viol) {
        viol = false;
        if (!filter[i]) {
            for (i = 0; i < Nele; ++i) {
                a = b = n = 0;
                for (j = 0; j < Nstar[i]; ++j) {
                    a += Rho[Nlist(i, j)];
                    b += Rho[Nlist(i, j)] * Rho[Nlist(i, j)];
                    n += 1;
                }
                if (n > 0) {
                    a /= (double) n; // average
                    if (a * a <= b / (double) n) {
                        b = std::sqrt(b / -a * a);
                        if ((Rho[i] - a) > std::sqrt(b * b + Rho_err[i] * Rho_err[i])) {
                            filter[i] = true;
                            viol = true;
                        }
                    }
                } else {
                    filter[i] = true;
                    viol = true;
                }
            }
        }
    }
}

void UDPClustering::merging() {
    int Nbarr = (Nclus * Nclus - Nclus) / 2; /// number of contacts between clusters
    int i, j, k, n, l, alive, dead;
    double c1, c2, b12;
    VecBool Survive(Nclus);
    bool change;
    VecDouble Barrier(Nbarr), Barrier_err;
    VecInt iBarrier(Nbarr);
    VecInt2d Bcorr(Nbarr, 2);
    VecInt M2O; /// conversion from merged to original cluster number
    VecInt O2M(Nclus); /// Conversion from original cluster number to its equivalent in merged

    n = 0;
    for (i = 0; i < Nclus - 1; ++i) {
        for (j = i + 1; j < Nclus; ++j) {
            n += 1;
            Barrier[n] = Bord(i, j);
            Barrier_err[n] = Bord_err(i, j);
            Bcorr(n, 0) = i;
            Bcorr(n, 1) = j;
        }
    }
    if (Nbarr > 1) {
        //TODO: call HPSORT(Nbarr,Barrier,iBarrier)
        std::sort(Barrier.data(), Barrier.data() + Barrier.size());
    } else {
        iBarrier[0] = 1;
    }

    Survive = true; // set all elements
    change = true;
    Cluster_m = VecInt(Nele); // TODO: realloc?
    while (change) {
        change = false;
        for (n = Nbarr; n > 1; --n) {
            k = iBarrier[n];
            i = Bcorr(k, 0);
            j = Bcorr(k, 1);
            if (Bord(i, j) > 0 && i != j && Survive[i] && Survive[j]) {
                c1 = cent[i] - sensibility * cent_err[i];
                c2 = cent[j] - sensibility * cent_err[j];
                b12 = Bord(i, j) + sensibility + Bord_err(i, j);

                if (c1 < b12 || c2 < b12) {
                    change = true;
                    Bord(i, j) = 0;
                    Bord_err(i, j) = 0;
                    if (c1 > c2) {
                        alive = i;
                        dead = j;
                    } else {
                        alive = j;
                        dead = i;
                    }
                    Bord(alive, alive) = 0;
                    Bord_err(alive, alive) = 0;
                    Survive[dead] = false;
                    for (k = 1; k < Nclus; ++k) {
                        if (Survive[k]) {
                            if (Bord(i, k) > Bord(j, k)) {
                                Bord(alive, k) = Bord(i, k);
                                Bord(k, alive) = Bord(k, i);
                                Bord_err(alive, k) = Bord_err(i, k);
                                Bord_err(k, alive) = Bord_err(k, i);
                            } else {
                                Bord(alive, k) = Bord(j, k);
                                Bord(k, alive) = Bord(k, j);
                                Bord_err(alive, k) = Bord_err(j, k);
                                Bord_err(k, alive) = Bord_err(k, j);
                            }
                        }
                    }
                    for (l = 0; l < Nele; ++l) {
                        if (Cluster_m[l] == dead) {
                            Cluster_m[l] = alive;
                        }
                    }
                    if (n > 1) {
                        for (l = n - 1; l > 1; --l) {
                            k = iBarrier[l];
                            if (Bcorr(k, 0) == dead) Bcorr(k, 0) = alive;
                            if (Bcorr(k, 1) == dead) Bcorr(k, 1) = alive;
                        }
                    }
                }
            }
        }
    }
    /// perform mapping of merged to original clusters
    Nclus_m = 0;
    for (i = 0; i < Nclus; ++i) {
        if (Survive[i]) Nclus_m += 1;
    }
    M2O = VecInt(Nclus_m);
    n = 0;
    O2M = -1;

    /// merged -> orginal
    for (i = 0; i < Nclus; ++i) {
        if (Survive[i]) {
            n += 1;
            M2O[n] = i;
            O2M[i] = n;
        }
    }

    /// get survival characteristics
    if (Nclus_m > 0) {
        for (i = 0; i < Nele; ++i) {
            Cluster_m[i] = O2M[Cluster_m[i]];
        }
    } else {
        Cluster_m = 1;
        throw ONLY_ONE_CLUSTER_AFTER_MERGING;
    }
    printf("end of get densities\n");
}

} // end of namespace udpclust