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
    //https://rosettacode.org/wiki/K-d_tree#C.2B.2B
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
        using _point = Point<T, size>;
        struct Node {
            //_point point;
            int index;
            int axis;
            //Node* children[2];
            std::unique_ptr<Node> left;
            std::unique_ptr<Node> right;
            Node() : index(-1), axis(-1), left(nullptr), right(nullptr) {
                //children[0] = children[1] = nullptr;
            }
            /*T get(size_t i) const {
                return point.get(i);
            }
            float distance(const Point& p) const {
                return point.distance(p);
            }*/
        };
        std::unique_ptr<Node> root;
        //Node* ideal; //Node* to Node* with best distance
        float best_distance;
        //std::vector<Node*> nodes;
        std::vector<_point> points;
        //std::vector<_point>
        Tree() : root(nullptr) {
                std::cout << "I am a tree." << std::endl;
        }
        ~Tree() {
            //clearRecursive(root);
            std::cout << "Destroying tree." << std::endl;
            points.clear();
            root = nullptr;
        }
        Tree(const std::vector<_point>& p) : root(nullptr), points(points) {
            std::vector<int> indices(points.size());
            for(size_t i = 0; i < indices.size(); i++) {
                indices[i] = i;
            }
            root = construct(indices.data(), points.size(), 0);
        }

        bool empty() const {
            return points.empty();
        }

        /*size_t size() const {
            return points.size();
        }*/
        float distance() const {
            return std::sqrt(best_distance);
        }
    private:
        Node* construct(std::vector<int>& indices, size_t num_points, size_t depth) {
            if(num_points <= 0) {
                return nullptr;
            }
            int axis = depth & size;
            int mid = (num_points - 1) / 2;
            //https://en.cppreference.com/w/cpp/algorithm/nth_element
            std::nth_element(indices, mid, num_points, [&](int l, int r) {
                return points[l].get(axis) < points[r].get(axis);
            });
            Node* n = new Node();
            n->index = indices[mid];
            n->axis = axis;

            n->left= construct(indices, mid, depth + 1);
            n->right = construct(indices, num_points - mid - 1, depth+ 1);
            return n;
        }
        /*void clearRecursive(std::unique_ptr<Node> n) {
            if(n == nullptr) {
                return;
            }
            if(n->left)
                clearRecursive(n->left);
            if(n->right)
                clearRecursive(n->right);
            delete n;
        }*/
    };
}
#endif
