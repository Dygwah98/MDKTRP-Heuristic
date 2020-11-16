#ifndef ALG_CONT_H
#define ALG_CONT_H

#include<vector>
#include"../n_implementations/GA.h"
#include"../n_implementations/GA_adaptive.h"
#include"../n_implementations/GA_all_equal.h"

class AlgorithmContainer {

    private:
        using Algorithm = double(*)(const Test&, const Individual&);

        double r;
        unsigned depots;
        unsigned customers;
        double **coordinate_matrix;
        distTable distances;
        const string dir = "./dat/";
        unsigned dpc;
        
        Algorithm algorithm;

    public:
        AlgorithmContainer( Algorithm alg): algorithm(alg) {}

        void execute(const TestInstances& t) {
            
            double mean_gap = 0.0;
            double gap_counter = 0;

            for(auto& instance : t.tests) {
    
                depots = 0;
                customers = 0;
                unsigned i = instance.start_number_file;  
                string number_instance = "";
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
                
                //cout << "file\n";
                read_file_ruiz(file, depots, customers, coordinate_matrix);
                dpc = depots + customers;
             
                distances.clear();

                Individual::setEnv(instance.vehicles, customers, depots, coordinate_matrix, &distances);

                Individual startingIndividual;
                
                r = (*algorithm)(instance, startingIndividual);

                free_matrix(coordinate_matrix, depots+customers);

                Individual::freeEnv();

                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

                cout << instance.prefix + number_instance << " : " << r <<  " : ";
                cout << (double)std::chrono::duration_cast
                        <std::chrono::milliseconds>(end - begin)
                        .count() / (double)1000 << " : ";
#ifdef BASE 
                if(instance.known_solution != 0) {

                    gap_counter += 1;
                    mean_gap += 100.0 * ((floor(r) - instance.known_solution) /instance.known_solution );
                    
                    cout << 100.0 * ((floor(r) - instance.known_solution) /instance.known_solution );
                }
                else
                    cout << 0;
#endif
                cout << "\n";
            }
#ifdef BASE
            cout << "GAP medio : " << mean_gap / gap_counter << "\n";
#endif
        }

        double getResult() const { return r; }
};

#endif