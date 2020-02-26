#include <thread>
#include <atomic>
#include <string>
#include <math.h>
#include <stdlib.h>
//#include <stdio.h>
#include <fstream>
#include "KDTree.hpp"

/*std::vector<unsigned long> validate_training_header(const char* training_file) {
    std::ifstream training(training_file, std::ios_base::binary);
    std::vector<unsigned long> ret;
    if (training) {
        training.seekg(0, training.end);
        size_t size = training.tellg();
        training.seekg(0, training.beg);
        char train_chars[8];
        training.read(train_chars, 8);
        std::cout << train_chars << std::endl;
        //training.read()
        unsigned long id, num_points, dims;
        training.read(id, 8);
        training.read(num_points, 8);
        training.read(dims, 8);
        training.close();
    }
    return ret;
}*/

int main(int argc, char** argv) {
    if (argc < 5) {
        printf("Usage: ./k-nn num_cores training_file query_file result_file\n");
    }
    else {
        size_t num_cores = static_cast<size_t>(atoi(argv[1]));
        const char* training_file = argv[2];
        const char* query_file = argv[3];
        const char* result_file = argv[4];
        std::ifstream training(training_file, std::ios_base::binary);
        std::vector<unsigned long> ret;
        unsigned long id, num_points, dims;
        std::vector<float> points; //
        if (training) {
            training.seekg(0, training.end);
            size_t size = training.tellg();
            training.seekg(0, training.beg);
            char train_chars[8];
            training.read(train_chars, 8);
            std::cout << train_chars << std::endl;
            //training.read()
            //training.read(id, 8);
            //training.read(num_points, 8);
            //training.read(dims, 8);

            //todo: populate a vector of floats from the training file
            training.close();
        }
        cs547::Tree<float, 4> my_tree;

    }
    return 0;
}
