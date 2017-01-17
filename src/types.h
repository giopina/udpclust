//
// Created by marscher on 1/12/17.
//

#ifndef UDPCLUST_TYPES_H
#define UDPCLUST_TYPES_H

#include <vector>
#include <numeric>
#include <cmath>
#include <sstream>
#include <iostream>
#include <memory>
#include <limits>


using std::size_t;

#undef NDEBUG
#include <cassert>


template<class T, typename D>
class MyArray1D {
    std::unique_ptr<T, D> uptr;
    size_t N;
public:
    MyArray1D() : N(0) {}
    MyArray1D(size_t N) : uptr(new T[N]), N(N) {}
    MyArray1D(T* raw_ptr, size_t N) : uptr(raw_ptr), N(N) {}

    MyArray1D(MyArray1D& other) = delete;
    MyArray1D &operator=(const MyArray1D &other) = delete;
    MyArray1D(MyArray1D&& other) : uptr(std::move(other.uptr)), N(other.N) {
        std::cout << "moved another array1d" << std::endl;
    }

    T &operator[](size_t i) const {
        if (i > N) {
            std::stringstream ss;
            ss << "index out of bounds: " << i << " exceeds " << N;
            throw std::runtime_error(ss.str());
        }
        return uptr.get()[i];
    }
    /// assign a value to the array
    void operator=(const T &value) {
        auto ptr = uptr.get();
        for (size_t i = 0; i < N; ++i) {
            ptr[i] = value;
        }
    }

    MyArray1D& operator=(MyArray1D&& other)
    {
        std::cout << "moved another array1d (assignment)" << std::endl;
        uptr = std::move(other.uptr);
        N = other.N;
        return *this;
    }

    size_t size() const { return N; }

    T *data() { return uptr.get(); };

    const T *data() const { return uptr.get(); }

    MyArray1D operator+(const MyArray1D &o) const {
        auto res = MyArray1D(o.size());
        auto ptr = uptr.get();

        for (size_t i = 0; i < o.size(); ++i) {
            res[i] = ptr[i] + o[i];
        }
        return res;
    }

    MyArray1D &operator=(T &value) {
        auto ptr = uptr.get();
        for (size_t i = 0; i < N; ++i) {
            ptr[i] = value;
        }
        return *this;
    }

    MyArray1D operator*(const MyArray1D &o) const {
        auto res = MyArray1D(o.size());
        auto ptr = uptr.get();
        for (size_t i = 0; i < N; ++i) {
            res[i] = ptr[i] * o[i];
        }
        return res;
    }

    T sum() const {
        T res = 0;
        auto ptr = uptr.get();
        for (size_t i = 0; i < N; ++i) {
            res += ptr[i];
        }
        return res;
    }

    // TODO: this only valid for T=double?
    MyArray1D pow(int exp) const {
        MyArray1D res(N);
        auto ptr = uptr.get();

        for (size_t i = 0; i < size(); ++i) {
            res[i] = std::pow(ptr[i], exp);
        }
        return res;
    }

    T min() const {
        T min = std::numeric_limits<T>::max();
        auto ptr = uptr.get();
        for (size_t i =0; i < size(); ++i) {
            if (ptr[i] < min) min = ptr[i];
        }
        return min;
    }

    T max() const {
        T max = std::numeric_limits<T>::min();
        T* ptr = uptr.get();
        for (size_t i = 0; i < size(); ++i) {
            if (ptr[i] > max) max = ptr[i];
        }
        return max;
    }

};

template<class T, class D>
class MyArray2D {
    std::unique_ptr<T, D> uptr;
    size_t rows, cols;
public:
    MyArray2D() : rows(0), cols(0) {}
    MyArray2D(T* raw_ptr, size_t rows, size_t cols) : uptr(raw_ptr), rows(rows), cols(cols) {}
    MyArray2D(size_t rows, size_t cols) : uptr(new T[rows * cols]), rows(rows), cols(cols) { }

    T &operator()(size_t i, size_t j) const { auto ptr = uptr.get(); return ptr[i + j * cols]; }

    void operator=(const T &value) {
        auto ptr = uptr.get();
        for (size_t i = 0; i < rows * cols; ++i) {
            ptr[i] = value;
        }
    }

    size_t size() const { return rows * cols; }

};

template<typename T>
struct no_delete {
    void operator()(T *__ptr) const {
        std::cout<< "dummy deleter" << std::endl;
    }
};

// internally managed memory arrays
typedef MyArray1D<bool, std::default_delete<bool>> VecBool;
typedef MyArray1D<int, std::default_delete<int>> VecInt;
typedef MyArray1D<double, std::default_delete<double>> VecDouble;

typedef MyArray2D<int, std::default_delete<int>> VecInt2d;
typedef MyArray2D<double, std::default_delete<double>> VecDouble2d;

// externally managed memory arrays
typedef MyArray1D<bool, no_delete<bool>> VecBoolExt;
typedef MyArray1D<int, no_delete<int>> VecIntExt;
typedef MyArray1D<double, no_delete<double>> VecDoubleExt;

typedef MyArray2D<int, no_delete<int>> VecIntExt2d;
typedef MyArray2D<double, no_delete<double>> VecDoubleExt2d;


#endif //UDPCLUST_TYPES_H
