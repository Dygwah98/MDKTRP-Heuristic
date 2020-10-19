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
        map<string, AlgorithmContainer > names = {
            {"0", AlgorithmContainer(new GA()) }
        };

        AlgorithmSelector() {};

    public:
        static void execute(string key) {

            AlgorithmSelector instance;
            if(instance.names.find(key) != instance.names.end()) {
                
                TestInstances test;
                AlgorithmContainer choice = instance.names.at(key);
                //execute algorithm.execute(test);

            }
        }

        ~AlgorithmSelector() {
            for(auto& el : names)
                delete el.second.getImpl();
        }
};

#endif