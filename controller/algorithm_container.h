#ifndef ALG_CONT_H
#define ALG_CONT_H

#include<vector>
#include"../implementations/base_algorithm.h"

class AlgorithmContainer {

    private:
        double r;
        BaseAlgorithm* algorithm;

    public:
        AlgorithmContainer( BaseAlgorithm* param ): algorithm(param) {}

        void execute(TestInstances t) {

            for(auto& instance : t.tests) {
    
                unsigned depots = 0;
                unsigned customers = 0;
                double **coordinate_matrix;
                string dir = "./dat/";
    
                for (unsigned i = instance.start_number_file; i <= instance.max_number_file; i++)
                {
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
                    double **distance_matrix = get_distance_matrix(depots, customers, coordinate_matrix);

                    Individual startingIndividual(instance.vehicles, depots, customers, distance_matrix);

                    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                    
                    r = (*algorithm)(instance, startingIndividual);

                    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

                    cout << instance.prefix + number_instance << ":" << r << ":" << (double)std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / (double)1000 << "\n";

                }
            }
        }

        BaseAlgorithm* getImpl() const { return algorithm; }

        double getResult() const { return r; }
};

#endif