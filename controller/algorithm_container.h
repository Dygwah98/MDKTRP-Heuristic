#ifndef ALG_CONT_H
#define ALG_CONT_H

#include<iomanip>
#include<vector>
#include"../implementation/algorithm/GA.h"
#include"../implementation/algorithm/GA_adaptive.h"
#include"../implementation/algorithm/GA_all_equal.h"

//wrapper dell'algoritmo, si occupa di inizializzare i dati ed eseguirlo per ogni istanza
class AlgorithmContainer {

    private:
        using Algorithm = double(*)(const Test&, const Individual&);

        //risultato dell'algoritmo
        double r;
        //numero di depots
        unsigned depots;
        //numero di customers
        unsigned customers;
        //depots + customers
        unsigned dpc;
        //matrice di coordinate [depots+customers][2]
        double **coordinate_matrix;
        //tabella di hash contenente le distanze
        distTable distances;
        const string dir = "./dat/";
        
        Algorithm algorithm;

    public:
        AlgorithmContainer( Algorithm alg): algorithm(alg) {}

        void execute(const TestInstances& t) {

            //per ogni istanza...
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
                
                //prendo in input le coordinate da giusto file
                read_file_ruiz(file, depots, customers, coordinate_matrix);
                dpc = depots + customers;
             
                //resetto le distanze (cambiano da istanza a istanza)
                distances.clear();

                //setto un pò di costanti condivise da tutti gli Individual 
                Individual::setEnv(instance.vehicles, customers, depots, coordinate_matrix, &distances);

                Individual startingIndividual;

                //eseguo l'algoritmo                
                r = (*algorithm)(instance, startingIndividual);

                //libero la memoria dinamica allocata
                free_matrix(coordinate_matrix, depots+customers);
                Individual::freeEnv();

                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

                //stampa dei risultati
                cout << instance.prefix + number_instance << " " << setprecision(2) << fixed << r <<  " ";
                cout << setprecision(2) << fixed
                     << (double)std::chrono::duration_cast
                        <std::chrono::milliseconds>(end - begin)
                        .count() / (double)1000;
                cout << "\n";
            }
        }

        double getResult() const { return r; }
};

#endif