#ifndef ALG_SEL_H
#define ALG_SEL_H

#include<string>
#include<map>
#include<vector>
#include<algorithm>
#include<iostream>
#include<iterator>
#include"algorithm_container.h"
using namespace std;

//si occupa di eseguire l'algoritmo corretto
class AlgorithmSelector {

    private:
        const map<unsigned, AlgorithmContainer > funcs = {
            {0, AlgorithmContainer(GeneticAlgorithm)         },
            {1, AlgorithmContainer(AllEqualGeneticAlgorithm) }
        };

        const map<unsigned, string > names = {
            {0, "GeneticAlgorithm"},
            {1, "AllEqualGeneticAlgorithm"}
        };

        AlgorithmSelector() {};

    public:
        static void execute(const unsigned key) {

            AlgorithmSelector instance;
            //se la key corrisponde a una delle funzioni disponibili...
            if(instance.funcs.find(key) != instance.funcs.end()) {
                
                //...seleziona l'algoritmo e inizializza la lista di istanze da eseguire
                TestInstances test;
                AlgorithmContainer choice = instance.funcs.at(key);
                
                cout << "seed " << seed; 
                cout << endl;
                cout << instance.names.at(key) << endl;
                cout << "Instance cost time" << endl;
                //esegue l'algoritmo scelto
                choice.execute(test);

            }
        }
};

#endif