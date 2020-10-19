#ifndef GA_H
#define GA_H

#include<vector>
#include"base_algorithm.h"

const int population_size = 250;

class GA : public BaseAlgorithm {
    
    private:
        unsigned tries = 1;
        unsigned mu;
        unsigned lambda = 250;
        unsigned max_g;
        unsigned mutator = SWAP2;
        unsigned crossover = TWO_POINT;
        unsigned max_evaluations_GA = 40000000 / 3;

        double body(Test instance, Individual ind) {
            double best_cost = std::numeric_limits<double>::max();

            //vettore che contiene la popolazione
            std::vector<Individual> individuals;
            individuals.assign(population_size, ind);

            Individual best_individual(ind);

            //inizializzazione della popolazione
            for (Individual &i : individuals)
            {
                i.random_inizialize();
                i.calculate_cost();
                if (i.get_cost() < best_cost)
                {
                    best_cost = i.get_cost();

                    best_individual = i;

                    if ((unsigned)best_individual.get_cost() == instance.known_solution)
                    {
                        return best_individual.get_cost();
                    }
                }
            }

            std::uniform_int_distribution<unsigned> random_parent(0, population_size / 2);
            std::uniform_int_distribution<unsigned> random_mut(0, 4);

            std::vector<Individual> new_generation;
            new_generation.assign(population_size, ind);

            std::vector<Individual> *new_generation_p = &new_generation;
            std::vector<Individual> *individuals_p = &individuals;

            unsigned g = 0;
            while (g < max_g)
            {
                //cout<<"G: "<<g<<"\n";
                const std::vector<Individual> &individuals_a = *individuals_p;
                std::vector<Individual> &new_generation_a = *new_generation_p;

                //elitismo al 10%
                const unsigned s = (10 * population_size) / 100;
                unsigned i = 0;
                for (; i < s; i++)
                    new_generation_a[i] = individuals_a[i];

                //dobbiamo inserire il rimanente 90% della popolazione
                const unsigned s2 = (90 * population_size) / 100;
                const unsigned s3 = s2 + s;
                for (; i < s3; i++)
                {
                    const unsigned p1 = random_parent(mt);
                    unsigned p2 = random_parent(mt);
                    while (p1 == p2)
                    {
                        p2 = random_parent(mt);
                    }

                    //Individual child(instance.vehicles, depots, customers, distance_matrix);

                    switch (crossover)
                    {
                    case 0:
                        new_generation_a[i].one_point_cross_over(individuals_a[p1], individuals_a[p2]);
                        break;
                    case 1:
                        new_generation_a[i].two_point_cross_over(individuals_a[p1], individuals_a[p2]);
                        break;
                    case 2:
                        new_generation_a[i].best_order_cross_over(individuals_a[p1], individuals_a[p2], best_individual);
                        break;
                    case 3:
                        new_generation_a[i].position_base_cross_over(individuals_a[p1], individuals_a[p2]);
                        break;
                    case 4:
                        new_generation_a[i].uniform_cross_over(individuals_a[p1], individuals_a[p2]);
                        break;
                    }

                    //mutazione genetica solo con un certo rateo
                    if (random_mut(mt) == 0)
                    {
                        switch (mutator)
                        {
                        case 0:
                            new_generation_a[i].swap2();
                            break;
                        case 1:
                            new_generation_a[i].swap3();
                            break;
                        case 2:
                            new_generation_a[i].scrumble();
                            break;
                        case 3:
                            new_generation_a[i].inversion();
                            break;
                        case 4:
                            new_generation_a[i].insertion();
                            break;
                        }
                    }

                    new_generation_a[i].calculate_cost();

                    if (new_generation_a[i].get_cost() < best_cost)
                    {
                        best_cost = new_generation_a[i].get_cost();
                        best_individual = new_generation_a[i];

                        //std::cout << "Child improved: " << best_cost << "\n";

                        if ((unsigned)best_individual.get_cost() == instance.known_solution)
                        {
                            return best_individual.get_cost();
                        }
                    }
                }

                std::vector<Individual> *const tmp = individuals_p;
                individuals_p = new_generation_p;
                new_generation_p = tmp;

                std::sort((*individuals_p).begin(), (*individuals_p).end());

                g++;
            }

            return best_individual.get_cost();
        }

    public:

    double operator()(Test instance, Individual ind) override {

        return body(instance, ind);
    }
};

#endif