#include<iostream>
#include"alg_selector.h"

int main(int argc, char const *argv[])
{
     if (argc < 2)
    {
        std::cerr << "At least one parameter for choosing the heuristic" << std::endl;
        return 1;
    }

    TestInstances test;

    AlgorithmSelector.execute(argv[1][0], test);

    return 0;
}
