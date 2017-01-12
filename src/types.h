//
// Created by marscher on 1/12/17.
//

#ifndef UDPCLUST_TYPES_H
#define UDPCLUST_TYPES_H

#include <vector>
#include <cstddef> // for size_t...

typedef std::vector<bool> VecBool;
typedef std::vector<int> VecInt;
typedef std::vector<double> VecDouble;

template <class T>
class MyArray1d {
    std::vector<T> data;
public:
    MyArray1d(size_t n) { data.reserve(n); }
};

template <class T>
class MyArray2D
{
    std::vector<T> data;
    size_t sizeX, sizeY;
public:
    MyArray2D(size_t rows, size_t cols):sizeX(rows), sizeY(cols) { data.reserve(rows*cols); }
    const T& at(int x, int y) const { return data.at(y + x * sizeY); }
    T& at(int x, int y) { return data.at(y + x * sizeY); }

    T& operator()(size_t i, size_t j) { return at(i, j ); }

    void operator=(const T& value) {
        data.assign(data.size(), value);
    }

// wrap other methods you need of std::vector here
};

typedef MyArray2D<double> VecDouble2d;
typedef MyArray2D<int> VecInt2d;


#endif //UDPCLUST_TYPES_H
