//
// Created by marscher on 1/16/17.
//

#include "udpclust.h"
#include <pybind11/pybind11.h>

using py = pybind11;

PYBIND11_PLUGIN(udpclust) {
        py::module m("UDP_modules.dp_clustering", "brief docstring");
        m.def("dp_advance", &udpclust::dp_advance, "doc");
        return m.ptr();
}