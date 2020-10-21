#ifndef ALG_CONT_H
#define ALG_CONT_H

#include<vector>
#include"../implementations/GA.h"
#include"../implementations/GA_adaptive.h"
#include"../implementations/GA_all_equal.h"

class AlgorithmContainer {

    private:
        using Algorithm = double(*)(const Test&, const Individual&);

        double r;
        unsigned depots;
        unsigned customers;
        double **coordinate_matrix;
        double **distance_matrix;
        const string dir = "./dat/";
        unsigned dpc;
        
        Algorithm algorithm;

    public:
        AlgorithmContainer( Algorithm alg): algorithm(alg) {}

        void execute(const TestInstances& t) {
            
            for(auto& instance : t.tests) {
    
                depots = 0;
                customers = 0;
                unsigned i = instance.start_number_file;  
                string number_instance;
                if (i < 10)
                {
                    number_instance = number_instance + "0" + to_string(i);
                }
                else
                {
                    number_instance = number_instance + to_string(i);
                }
                string file = dir + instance.prefix + number_instance + ".txt";

                std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                
                read_file_ruiz(file, depots, customers, coordinate_matrix);
                dpc = depots + customers;
                distance_matrix = get_distance_matrix(depots, customers, coordinate_matrix);
                
                Individual startingIndividual(instance.vehicles, depots, customers, distance_matrix);
                
                r = (*algorithm)(instance, startingIndividual);

                free_matrix(distance_matrix, depots+customers);
                free_matrix(coordinate_matrix, depots+customers);

                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

                cout << instance.prefix + number_instance << ":" << r << ":" << (double)std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / (double)1000 << "\n";

            }
        }

        double getResult() const { return r; }
};

#endif