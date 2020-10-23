#ifndef GA_EQ_H
#define GA_EQ_H

#include"GA.h"

class AllEqualGeneticAlgorithmData {

    public:
        const unsigned tries = 1;
        const unsigned population_size = 250;
        const unsigned mutator = SWAP2;
        const unsigned crossover = TWO_POINT;
        const unsigned max_evaluations_GA = 40000000 / 3;
        const unsigned mut_rate = 4;
        const unsigned cross_rate = 4;
        const unsigned mutator_rate = 4;
};

double AllEqualGeneticAlgorithm(const Test& instance, const Individual& ind) {

    const AllEqualGeneticAlgorithmData gdata;
    
    const unsigned max_evaluations = gdata.max_evaluations_GA * instance.factor_valuations;
    const unsigned max_g = (max_evaluations / gdata.tries) / gdata.population_size;
    
    double cost = std::numeric_limits<double>::max();
    double global_best = std::numeric_limits<double>::max();

    for (unsigned i = 0; i < gdata.tries; ++i)
    {
        double best_cost = std::numeric_limits<double>::max();

        //vettore che contiene la popolazione
        std::vector<Individual> individuals;
        individuals.assign(gdata.population_size, ind);

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
                    cout << "known solution ";
                    return best_individual.get_cost();                    
                }
            }
        }

        //std::sort(individuals.begin(), individuals.end());

        const unsigned s = (10 * gdata.population_size) / 100;
        const unsigned s2 = (90 * gdata.population_size) / 100;
        const unsigned s3 = s2 + s;  

        std::uniform_int_distribution<unsigned> random_parent(0, gdata.population_size / 2);
        std::uniform_int_distribution<unsigned> random_mut(0, gdata.mut_rate);
        std::uniform_int_distribution<unsigned> random_cross(0, gdata.cross_rate);
        std::uniform_int_distribution<unsigned> random_mutator(0, gdata.mutator_rate);

        std::vector<Individual> new_generation;
        new_generation.assign(gdata.population_size, ind);

        std::vector<Individual> *new_generation_p = &new_generation;
        std::vector<Individual> *individuals_p = &individuals;

        unsigned g = 0;
        
        unsigned p1;
        unsigned p2; 
        while (g < max_g)
        {
            //cout<<"G: "<<g<<"\n";
            const std::vector<Individual> &individuals_a = *individuals_p;
            std::vector<Individual> &new_generation_a = *new_generation_p;

            //elitismo al 10%
            unsigned i = 0;
            for (; i < s; ++i)
                new_generation_a[i] = individuals_a[i];

            //dobbiamo inserire il rimanente 90% della popolazione
            for (; i < s3; ++i)
            {
                p1 = random_parent(mt);
                p2 = random_parent(mt);
                while (p1 == p2)
                {
                    p2 = random_parent(mt);
                }

                //Individual child(instance.vehicles, depots, customers, distance_matrix);

                const unsigned cross = random_cross(mt);
                switch (cross)
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
                    const unsigned mut = random_mutator(mt);
                    switch (mut)
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
                        cout << "known solution ";
                        return best_individual.get_cost();
                    }
                }
            }

            std::vector<Individual> *const tmp = individuals_p;
            individuals_p = new_generation_p;
            new_generation_p = tmp;

            std::sort((*individuals_p).begin(), (*individuals_p).end());

            ++g;
        }

        cost = best_individual.get_cost();
        
         if (cost < global_best)
        {
            global_best = cost;
            if ((unsigned)global_best == instance.known_solution)
            {
                cout << "known_solution ";
                return global_best;
            }
        }
    }

    return cost;
}

#endif