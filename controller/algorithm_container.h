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

            //Individual startingIndividual;
            //input goes here

            for(auto& instance : t.tests) {
                //r = (*algorithm)(instance, startingIndividual);
            }
        }

        BaseAlgorithm* getImpl() { return algorithm; }
};

#endif