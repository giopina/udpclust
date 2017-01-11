//
// Created by marscher on 1/11/17.
//
#ifndef UDPCLUST_HEAP_SORT_H
#define UDPCLUST_HEAP_SORT_H

namespace {
template<typename T>
void checkRootNode(T *array, size_t root, size_t size) {
    size_t left = 2 * root;
    size_t right = 2 * root + 1;
    if (left < size && array[root] < array[left]) {
        std::swap(array[root], array[left]);
        checkRootNode(array, left, size);
    }
    if (right < size && array[root] < array[right]) {
        std::swap(array[root], array[right]);
        checkRootNode(array, right, size);
    }
}

template<typename T>
void buildHeap(T *array, size_t size) {
    for (size_t i = size / 2; i > 0; --i) {
        checkRootNode(array, i, size);
    }
}
}

namespace heap_sort {

template<typename T>
void sort(T *array, size_t size) {
    buildHeap(array, size);
    while (size > 1) {
        std::swap(array[1], array[size - 1]);
        checkRootNode(array, 1, --size);
    }
}

}

#endif //UDPCLUST_HEAP_SORT_H
