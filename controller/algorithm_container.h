#ifndef ALG_CONT_H
#define ALG_CONT_H

#include<vector>
#include"../implementations/GA.h"

class AlgorithmContainer {

    private:
        using Algorithm = double(*)(const Test&, const Individual&, const BaseData&);

        double r;
        unsigned depots;
        unsigned customers;
        double **coordinate_matrix;
        const string dir = "./dat/";
        unsigned dpc;
        
        Algorithm algorithm;
        const BaseData& data;

    public:
        AlgorithmContainer( Algorithm alg, const BaseData& param ): algorithm(alg), data(param) {}

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

                read_file_ruiz(file, depots, customers, coordinate_matrix);
                dpc = depots + customers;
                double **distance_matrix = get_distance_matrix(depots, customers, coordinate_matrix);

                Individual startingIndividual(instance.vehicles, depots, customers, distance_matrix);

                std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                
                r = (*algorithm)(instance, startingIndividual, data);

                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

                cout << instance.prefix + number_instance << ":" << r << ":" << (double)std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / (double)1000 << "\n";

                free_matrix(distance_matrix, depots+customers);
                free_matrix(coordinate_matrix, depots+customers);
            }
        }

        double getResult() const { return r; }
};

#endif