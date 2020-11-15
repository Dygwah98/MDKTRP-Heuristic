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
#ifndef BASE
        double *activation_costs = nullptr;
#endif
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
                

                read_file_ruiz(file, depots, customers, coordinate_matrix);
                dpc = depots + customers;
             
                distances.clear();
#ifndef BASE
                calculate_activation_costs(instance.vehicles);
                Individual::setEnv(instance.vehicles, customers, depots, activation_costs, coordinate_matrix, &distances);
#else
                Individual::setEnv(instance.vehicles, customers, depots, coordinate_matrix, &distances);
#endif
                Individual startingIndividual;
                
                r = (*algorithm)(instance, startingIndividual);

                free_matrix(coordinate_matrix, depots+customers);

                Individual::freeEnv();

#ifndef BASE
                delete[] activation_costs;
#endif

                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

                cout << instance.prefix + number_instance << " : " << r <<  " : ";
                cout << (double)std::chrono::duration_cast
                        <std::chrono::milliseconds>(end - begin)
                        .count() / (double)1000 << " : ";
#ifdef BASE 
                if(instance.known_solution != 0) {

                    gap_counter += 1;
                    mean_gap += 100.0 * ((r - instance.known_solution) /r );
                    
                    cout << 100.0 * ((r - instance.known_solution) /r );
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

        inline unsigned getCustomerIndex(unsigned x, unsigned y) {
            return depots*customers + x*customers + y;
        }

        inline unsigned getDepotIndex(unsigned x, unsigned y) {
            return x*customers + y;
        }

#ifndef BASE
        void calculate_activation_costs(unsigned vehicles) {

            activation_costs = new double[depots];

            //per ogni depot
            for(unsigned i = 0; i < depots; ++i) {

                double mean = 0;

                double max_cost = euclidean_distance(coordinate_matrix[i][0], coordinate_matrix[i][1],
													 coordinate_matrix[depots][0], coordinate_matrix[depots][1]);
                
                //si inserisce la distanza fra il depot e il primo cliente
                unsigned index = getDepotIndex(i, 0);
                distances[index] = max_cost;

                mean += max_cost;
                
                //per ogni cliente dopo il primo
                for(unsigned j = depots + 1; j < depots + customers; ++j) {
                    
                    double cost = euclidean_distance(coordinate_matrix[i][0], coordinate_matrix[i][1],
													 coordinate_matrix[j][0], coordinate_matrix[j][1]);

                    mean += cost;
                    //si inserisce la distanza fra il depot e il cliente
                    index = getDepotIndex(i, j - depots);
                    distances[index] = cost;

                    //se è la distanza massima trovata finora, conservala
                    if(cost > max_cost) {
                        max_cost = cost;
                    }
                }

                //calcolo della media e formula del costo d'attivazione
                mean /= (double)customers;
                double result = (max_cost - mean) * customers;

                //il costo d'attivazione viene salvato
                activation_costs[i] = result/vehicles;
            }

            //cout << "   activation costs:\n";
            //for(unsigned i = 0; i < depots; ++i) {
            //    cout << "       depot: " << i << " " << activation_costs[i] << endl;
            //}
        }
#endif

        double getResult() const { return r; }
};

#endif