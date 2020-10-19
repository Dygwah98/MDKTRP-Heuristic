#ifndef BA_H
#define BA_H

#include"test.h"
#include"individual.h"
#include"utils.h"

class BaseAlgorithm {

    protected:
        double cost;
        double global_best;

    public:
        virtual double operator()(Test t, Individual ind) = 0;
};

#endif