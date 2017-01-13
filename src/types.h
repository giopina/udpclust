//
// Created by marscher on 1/12/17.
//

#ifndef UDPCLUST_TYPES_H
#define UDPCLUST_TYPES_H

#include <vector>
#include <cstddef> // for size_t...
#include <numeric>
#include <cassert>

template <class T>
class MyArray2D
{
    std::vector<T> data;
    size_t rows, cols;
public:
    MyArray2D(size_t rows, size_t cols) : rows(rows), cols(cols) { data.reserve(rows*cols); }
    const T& at(int x, int y) const { return data.at(y + x * cols); }
    T& at(int x, int y) { return data.at(y + x * cols); }

    T& operator()(size_t i, size_t j) { return at(i, j ); }

    void operator=(const T& value) {
        data.assign(data.size(), value);
    }

    size_t size() const { return data.size(); }
};


typedef std::vector<bool> VecBool;
typedef std::vector<int> VecInt;
typedef std::vector<double> VecDouble;

typedef MyArray2D<double> VecDouble2d;
typedef MyArray2D<int> VecInt2d;

template <class T>
ostream& operator*(ostream& out, Pair<T,U>& v);

template <class T>
T sum(const std::vector<T>& vec) {
    return std::accumulate(vec.begin(), vec.end(), T());
}

template<class T>
std::vector<T> pow(const std::vector& vec, int exp) {
    auto res = std::vector<T>(vec.size());
    for(size_t i = 0; i < vec.size(); ++i) {
        res[i] = std::pow(vec[i], exp);
    }
    return res;
}

template <class T>
std::vector<T> mult(const std::vector<T>& a, const std::vector<T>& b ) {
    assert(a.size(), b.size());
    auto res = std::vector<T>(vec.size());
    for(size_t i = 0; i < vec.size(); ++i) {
        res[i] = a[i] * b[i];
    }
    return res;
}



// exception
/*
class udp_exception: std::runtime_error {
public:
    udp_exception(int error_code) : std::runtime_error(  )
};
*/


#endif //UDPCLUST_TYPES_H
