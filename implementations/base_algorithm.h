#ifndef BA_H
#define BA_H

#include"test.h"
#include"individual.h"
#include"utils.h"

class BaseAlgorithm {

    public:
        virtual double operator()(Test t, Individual ind) = 0;
};

#endif