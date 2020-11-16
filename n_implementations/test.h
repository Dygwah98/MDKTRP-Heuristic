#ifndef TEST_H
#define TEST_H

#include<string>
#include<vector>
using namespace std;

class Test {

    public:
        unsigned vehicles;
        string prefix;
        unsigned start_number_file;
        unsigned max_number_file;
        double known_solution;
        unsigned factor_valuations;

        Test(   const unsigned vehicles, 
                string prefix, 
                const unsigned start_number_file, 
                const unsigned max_number_file, 
                const double known_solution, 
                const unsigned factor_valuations) {
            
            this->vehicles = vehicles;
            this->prefix = prefix;
            this->start_number_file = start_number_file;
            this->max_number_file = max_number_file;
            this->known_solution = known_solution;
            this->factor_valuations = factor_valuations;
        }
};

class TestInstances {

    public:
        vector<Test> tests = {
            {5, "Ir", 1, 1, 545, 1},
            {5, "Ir", 2, 2, 832, 1},
            {5, "Ir", 3, 3, 832, 1},
            {10, "Ir", 4, 4, 2082, 1},
            {10, "Ir", 5, 5, 1827, 1},
            {10, "Ir", 6, 6, 1786, 1},
            {20, "Ir", 7, 7, 5424, 1},
            {20, "Ir", 8, 8, 3737, 1},
            {20, "Ir", 9, 9, 3082, 1},
            {20, "Ir", 10, 10, 2786, 1},
            {20, "Ir", 11, 11, 2909, 1},
            {20, "Ir", 12, 12, 3171, 1},
            {25, "Ir", 13, 13, 8288, 1},
            {25, "Ir", 14, 14, 7257, 1},
            {25, "Ir", 15, 15, 8625, 1},
            {25, "Ir", 16, 16, 5265, 1},
            {25, "Ir", 17, 17, 6107, 1},
            {25, "Ir", 18, 18, 5788, 1},
            {35, "p", 1, 1, 660, 1},
            {35, "p", 2, 2, 660, 1},
            {35, "p", 3, 3, 906, 1},
            {35, "p", 4, 4, 1881, 1},
            {35, "p", 5, 5, 1871, 1},
            {35, "p", 6, 6, 1460, 1},
            {35, "p", 7, 7, 1453, 1},
            {35, "p", 8, 8, 0, 15},
            {35, "p", 9, 9, 0, 15},
            {35, "p", 10, 10, 0, 15},
            {35, "p", 11, 11, 0, 15},
            {35, "p", 12, 12, 2769, 1},
            {35, "p", 15, 15, 5618, 10},
            {35, "p", 18, 18, 0, 20},
            {35, "pr", 1, 1, 1167, 1},
            {35, "pr", 2, 2, 2422, 2},
            {35, "pr", 3, 3, 4287, 5},
            {35, "pr", 4, 4, 5582, 10},
            {35, "pr", 5, 5, 6782, 15},
            {35, "pr", 6, 6, 0, 20},
            {35, "pr", 7, 7, 1594, 1},
            {35, "pr", 8, 8, 3817, 5},
            {35, "pr", 9, 9, 5668, 10},
            {35, "pr", 10, 10, 0, 20}
        };
};

#endif