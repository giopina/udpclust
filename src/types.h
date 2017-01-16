//
// Created by marscher on 1/12/17.
//

#ifndef UDPCLUST_TYPES_H
#define UDPCLUST_TYPES_H

#include <vector>
#include <cstddef> // for size_t...
#include <numeric>
#include <cassert>
#include <cmath>


template<class T>
class MyArray1D {
    T *ptr;
    size_t N;
    bool own_data;
public:
    MyArray1D() : ptr(0), N(0), own_data(true) {}
    MyArray1D(T *ptr, size_t N) : ptr(ptr), N(N), own_data(false) { if(ptr == nullptr) { ptr = new T[N]; }}
    MyArray1D(size_t N) : ptr(new T[N]), N(N), own_data(true) {}

    virtual ~MyArray1D() { if (own_data) { delete[] ptr; }}

    T &operator[](size_t i) const { return ptr[i]; }

    void operator=(const T &value) {
        for (size_t i = 0; i < N; ++i) {
            ptr[i] = value;
        }
    }

    MyArray1D &operator=(const MyArray1D &other) {
        if (this != &other) { // protect against invalid self-assignment
            // 1: allocate new memory and copy the elements
            T *new_array = new T[other.size()];
            std::copy(other.data(), other.data() + other.size(), new_array);

            // 2: deallocate old memory
            delete[] ptr;

            // 3: assign the new memory to the object
            ptr = new_array;
            N = other.size();
            own_data = other.own_data;
        }
        // by convention, always return *this
        return *this;
    }

    size_t size() const { return N; }

    T *data() { return ptr; };

    const T *data() const { return ptr; }

    MyArray1D operator+(const MyArray1D &o) const {
        auto res = MyArray1D(o.size());
        for (size_t i = 0; i < o.size(); ++i) {
            res[i] = ptr[i] + o[i];
        }
        return res;
    }

    MyArray1D &operator=(T &value) {
        for (size_t i = 0; i < N; ++i) {
            ptr[i] = value;
        }
        return *this;
    }

    MyArray1D operator*(const MyArray1D &o) const {
        auto res = MyArray1D(o.size());
        for (size_t i = 0; i < N; ++i) {
            res[i] = ptr[i] * o[i];
        }
        return res;
    }

    T sum() const {
        T res = 0;
        for (size_t i = 0; i < N; ++i) {
            res += ptr[i];
        }
        return res;
    }

    // TODO: this only valid for T=double?
    MyArray1D pow(int exp) const {
        auto res = *this; // huh?
        for (size_t i = 0; i < size(); ++i) {
            res[i] = std::pow(ptr[i], exp);
        }
        return res;
    }

};

template<class T>
class MyArray2D {
    T *data;
    size_t rows, cols;
    bool own_data;
public:
    MyArray2D() : data(nullptr), own_data(true), rows(0), cols(0) {}
    MyArray2D(T *ptr, size_t rows, size_t cols) : data(ptr), rows(rows), cols(cols), own_data(false) {}
    MyArray2D(size_t rows, size_t cols) : rows(rows), cols(cols), own_data(true) { data = new T[rows * cols]; }

    virtual ~MyArray2D() { if (own_data) { delete[] data; }}

    T &operator()(size_t i, size_t j) const { return data[i + j * cols]; }

    void operator=(const T &value) {
        for (size_t i = 0; i < rows * cols; ++i) {
            data[i] = value;
        }
    }

    size_t size() const { return rows * cols; }

};


typedef MyArray1D<bool> VecBool;
typedef MyArray1D<int> VecInt;
typedef MyArray1D<double> VecDouble;

typedef MyArray2D<double> VecDouble2d;
typedef MyArray2D<int> VecInt2d;

/*
template <class T>
T sum(const std::vector<T>& vec) {
    return std::accumulate(vec.begin(), vec.end(), T());
}*/


#endif //UDPCLUST_TYPES_H
