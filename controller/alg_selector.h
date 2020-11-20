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

class AlgorithmSelector {

    private:
        const map<unsigned, AlgorithmContainer > funcs = {
            {0, AlgorithmContainer(GeneticAlgorithm)         },
            {1, AlgorithmContainer(AllEqualGeneticAlgorithm) },
            {2, AlgorithmContainer(AdaptiveGeneticAlgorithm) }
        };

        const map<unsigned, string > names = {
            {0, "GeneticAlgorithm"},
            {2, "AdaptiveGeneticAlgorithm"},
            {1, "AllEqualGeneticAlgorithm"}
        };

        AlgorithmSelector() {};

    public:
        static void execute(const unsigned key) {

            AlgorithmSelector instance;
            if(instance.funcs.find(key) != instance.funcs.end()) {
                
                TestInstances test;
                AlgorithmContainer choice = instance.funcs.at(key);
                
                cout << "seed " << seed; 
                cout << endl;
                cout << instance.names.at(key) << endl;
                cout << "Instance cost time" << endl;
                choice.execute(test);

            }
        }
};

#endif