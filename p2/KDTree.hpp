#ifndef BOTTCHER_K_D_TREE_HPP
#define BOTTCHER_K_D_TREE_HPP
#include <iostream>
#include <array>
#include <random>
#include <cmath>
#include <algorithm>
#include <vector>
#include <initializer_list>

namespace cs547 {
    //todo: create helper functions for files, i.e. is valid header for training, is valid query file, etc.
    template<typename T, size_t size>
    struct Point {
        Point() = delete;
        Point(std::array<T, size> coordinates) : coordinates(coordinates) {}
        Point(std::initializer_list<T> list) {
            for(auto &i = 0; i < size; i++) {
                coordinates[i] = list[i];
            }
        }
        friend std::ostream& operator<<(std::ostream& os, const Point& p) {
            os << "(";
            for(size_t i = 0; i < size; i++) {
                if(i > 0)
                    os << ", ";
                os << p.get(i);
            }
            os << ")";
            return os;
        }
        T get(size_t i) const { //want this to be immutable so we don't clown ourselves
            return coordinates[i];
        }
        float distance(const Point& p) { //we are using 32-bit floating point numbers in this assignment
            float distance = 0;
            for(size_t i = 0; i < size; i++) {
                float d = this->get(i) - p.get(i);
                distance += d * d;
            }
            return distance; //return the distance squared
        }
    private:
        std::array<T, size> coordinates;
    };
    //end class Point
    template<typename T, size_t size>
    struct Tree {
        using point = Point<T, size>;
        Tree() {

        }
    };
}
#endif
