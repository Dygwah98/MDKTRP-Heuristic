#ifndef ALG_SEL_H
#define ALG_SEL_H

#include<string>
#include<map>
#include<vector>
#include<algorithm>
#include"algorithm_container.h"
#include"../implementations/GA.h"
using namespace std;

class AlgorithmSelector {

    private:
        const map<unsigned, AlgorithmContainer > names = {
            {0, AlgorithmContainer(new GA()) }
        };

        AlgorithmSelector() {};

    public:
        static void execute(const unsigned key) {

            AlgorithmSelector instance;
            if(instance.names.find(key) != instance.names.end()) {
                
                TestInstances test;
                AlgorithmContainer choice = instance.names.at(key);
                choice.execute(test);

            }
        }

        ~AlgorithmSelector() {
            for(auto& el : names)
                delete el.second.getImpl();
        }
};

#endif