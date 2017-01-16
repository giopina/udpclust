//
// Created by marscher on 1/16/17.
//

#include "udpclust.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

/**
 * UDP_modules.dp_clustering.get_densities\
    (id_err,dmat,dim,rho,rho_err,filt,Nlist,Nstar)
 * @param id_err
 * @param dist_mat
 * @param Cluster
 * @param Rho
 * @param filter
 * @param dimint
 * @param Nlist
 * @param Nele
 * @param sensitivity
 * @param maxknn
 */
//** todo: replace arrays with basic data types where applicable (id_err, Nele, dimint, maxknn, Nstar)
void get_densities(py::array_t<int> id_err,
                   py::array_t<double> dist_mat,
                   py::array_t<int> Nele,
                   py::array_t<int> dimint,
                   py::array_t<double> Rho,
                   py::array_t<double> Rho_err,
                   py::array_t<bool> filter,
                   py::array_t<int> Nlist,
                   py::array_t<int> Nstar,
                   int maxknn
) {
    auto buf1 = dist_mat.request();
    auto buf3 = Rho.request();

    auto buf4 = filter.request();
    auto buf5 = Nlist.request();
    auto buf6 = Rho_err.request();
    auto buf7 = Nele.request();
    auto buf8 = dimint.request();
    auto buf9 = Nstar.request();

    double *dist_mat_p = (double *) buf1.ptr;
    double *rho_p = (double *) buf3.ptr;
    double *rho_err_p = (double *) buf6.ptr;
    bool *filter_p = (bool *) buf4.ptr;
    int *Nlist_p = (int *) buf5.ptr;
    int *Nstar_p = (int *) buf9.ptr;

    int res = -1;
    /// void get_densities(int* id_err, double* dist_mat, int Nele, int dimint, double* Rho, double* Rho_err,
    //                     bool* filter, int* Nlist, int* Nstar, int maxknn)
    int iNele = ((int *) &buf7.ptr)[0];
    int idimint = ((int *) &buf8.ptr)[0];
    udpclust::get_densities(&res, dist_mat_p, iNele, idimint, rho_p,
                            rho_err_p, filter_p, Nlist_p, Nstar_p, maxknn);

    id_err[0] = res;

}


PYBIND11_PLUGIN(dp_clustering) {
    py::module m("dp_clustering", "brief docstring");
    m.def("dp_advance", &udpclust::dp_advance, "doc");
    m.def("get_densities", &get_densities, "");
    return m.ptr();
}
