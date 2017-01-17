//
// Created by marscher on 1/12/17.
//

#ifndef UDPCLUST_TYPES_H
#define UDPCLUST_TYPES_H

#include <vector>
//#include <cstddef> // for size_t...
#include <numeric>
#include <cassert>
#include <cmath>
#include <memory>

using std::size_t;


template<class T, typename D>
class MyArray1D {
    // T *ptr;
    std::unique_ptr<T, D> uptr;
    size_t N;
public:
    MyArray1D() : N(0) {}
    MyArray1D(size_t N) : uptr(new T[N]), N(N) {}
    MyArray1D(T* raw_ptr, std::size_t N) : uptr(raw_ptr), N(N) {}

    MyArray1D(const MyArray1D& other) = 0;
    MyArray1D(const MyArray1D&& other) : uptr(std::move(other.uptr)), N(other.N) {
    }

    T &operator[](size_t i) const { return uptr.get()[i]; }

    void operator=(const T &value) {
        auto ptr = uptr.get();
        for (size_t i = 0; i < N; ++i) {
            ptr[i] = value;
        }
    }

    MyArray1D &operator=(const MyArray1D &other) {
        if (this != &other) { // protect against invalid self-assignment
            /*// 1: allocate new memory and copy the elements
            T *new_array = new T[other.size()];
            std::copy(other.data(), other.data() + other.size(), new_array);

            // 2: deallocate old memory
            if (own_data)
                delete[] ptr;

            // 3: assign the new memory to the object
            ptr = new_array;
            N = other.size();
            own_data = true;*/
            //std::swap
        }
        // by convention, always return *this
        return *this;
    }
    MyArray1D& operator=(MyArray1D&& other)
    {
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
        auto res = *this; // make copy
        auto ptr = uptr.get();

        for (size_t i = 0; i < size(); ++i) {
            res[i] = std::pow(ptr[i], exp);
        }
        return res;
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
    void operator()(T *__ptr) const {}
};

// internally managed memory arrays
typedef MyArray1D<bool, std::default_delete<bool[]>> VecBool;
typedef MyArray1D<int, std::default_delete<int[]>> VecInt;
typedef MyArray1D<double, std::default_delete<double[]>> VecDouble;

typedef MyArray2D<int, std::default_delete<int[]>> VecInt2d;
typedef MyArray2D<double, std::default_delete<double[]>> VecDouble2d;

// externally managed memory arrays
typedef MyArray1D<bool, no_delete<bool>> VecBoolExt;
typedef MyArray1D<int, no_delete<int>> VecIntExt;
typedef MyArray1D<double, no_delete<double>> VecDoubleExt;

typedef MyArray2D<int, no_delete<int>> VecIntExt2d;
typedef MyArray2D<double, no_delete<double>> VecDoubleExt2d;


#endif //UDPCLUST_TYPES_H
