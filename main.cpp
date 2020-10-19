#include<iostream>
#include"controller/alg_selector.h"

int main(int argc, char const *argv[])
{
     if (argc < 2)
    {
        std::cerr << "At least one parameter for choosing the heuristic" << std::endl;
        return 1;
    }

    AlgorithmSelector::execute( argv[1][0] - '0' );

    return 0;
}
