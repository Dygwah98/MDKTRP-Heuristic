#ifndef GA_H
#define GA_H

#include"individual/misc.h"

class GeneticAlgorithmData {

    public:
        const unsigned tries = 1;
        const unsigned population_size = 500;
        const unsigned mutator = SCRAMBLE;
        const unsigned crossover = TWO_POINT;
        const unsigned max_evaluations_GA = 40000000 / 3;
        const unsigned mut_rate = 2;
};

double GeneticAlgorithm(const Test& instance, const Individual& ind) {

    //variabili per misurare le prestazioni
    
    Analyzer printStats; 
    Timer& initialize = printStats.initialize;
    Timer& crossover = printStats.crossover;
    Timer& mutation = printStats.mutation;
    Timer& improvement = printStats.improvement;
    Timer& repair = printStats.repair;
    Timer& costs = printStats.costs;

    //schema genetico

    const GeneticAlgorithmData gdata;
    
    const unsigned max_evaluations = gdata.max_evaluations_GA * instance.factor_valuations;
    const unsigned max_g = (max_evaluations / gdata.tries) / gdata.population_size;
    const unsigned popsize = gdata.population_size;
    
    double cost = std::numeric_limits<double>::max();
    double global_best = std::numeric_limits<double>::max();

    Individual best_individual2(ind);
    
    initialize.random_initialize(best_individual2);

    //cout << "       parameters set\n";
    const unsigned max_tries = gdata.tries;

    std::vector<unsigned> I;
    for(unsigned i = 0; i < popsize; ++i)
        I.push_back(i);

    for (unsigned tries_i = 0; tries_i < max_tries; ++tries_i)
    {
        double best_cost = std::numeric_limits<double>::max();

        //vettore che contiene la popolazione
        std::vector<Individual> individuals_original;
        individuals_original.assign(popsize, ind);

        Individual best_individual(best_individual2);
        
        double mean_cost = 0;
        //inizializzazione della popolazione
        for (Individual &i : individuals_original)
        {
            initialize.random_initialize(i);

            costs.calculate_cost(i);

            double cost = i.get_cost(); 
            if (cost < best_cost)
            {
                best_cost = cost;
                best_individual = i;
                best_individual2 = i;

                if (cost == instance.known_solution)
                {
                    printStats();
                    cout << "known solution ";
                    //best_individual.print_tour();
                    return cost;                    
                }
            }

            mean_cost += cost;
        }
        //cout << "       population initialized\n";

        mean_cost /= (double)popsize;

        const unsigned s  = (10 * popsize) / 100;
        
        std::uniform_int_distribution<unsigned> random_parent(1, popsize/2 - 2);
        std::uniform_int_distribution<unsigned> random_mut(0, gdata.mut_rate);
        std::uniform_int_distribution<unsigned> random_choice(0, 1);
        std::uniform_int_distribution<unsigned> substitution(0, 10);

        std::vector<Individual> new_generation_original;
        new_generation_original.assign(popsize, ind);

        keySort sorts(individuals_original, new_generation_original);

        std::sort(I.begin(), I.end(), sorts);
        sorts._swap();

        unsigned g = 0;
        
        unsigned p1;
        unsigned p2; 

        std::vector<Individual> *individuals_ptr = &individuals_original;
        std::vector<Individual> *new_generation_ptr = &new_generation_original; 

        //cout << "       starting evaluations\n";
        while (g < max_g)
        {
            double new_mean_cost = 0;
            //cout<<"G: "<<g<<"\n";

            const std::vector<Individual> &individuals = *individuals_ptr;
            std::vector<Individual> &new_generation = *new_generation_ptr;

            //elitismo al 10%
            unsigned i = 0;
            for (; i < s; ++i) {
                new_generation[ I[i] ] = individuals[ I[i] ];
                
                new_mean_cost += new_generation[i].get_cost();
            }

            //cout << "       iteration " << g << ": elite population processed\n";

            //dobbiamo inserire il rimanente 90% della popolazione
            for (; i < popsize; ++i)
            {
                        
                p1 = random_parent(mt);
                std::uniform_int_distribution<unsigned> left(0, p1-1);
                std::uniform_int_distribution<unsigned> right(p1 + 1, popsize - 1);
                if(random_choice(mt))
                    p2 = left(mt);
                else
                    p2 = right(mt);

                switch (gdata.crossover)
                {
                case 0:

                    crossover.one_point_crossover(new_generation[ I[i] ], individuals[ I[p1] ], individuals[ I[p2] ]);
                    break;
                case 1:

                    crossover.two_point_crossover(new_generation[ I[i] ], individuals[ I[p1] ], individuals[ I[p2] ]);
                    break;
                case 2:

                    crossover.best_order_crossover(new_generation[ I[i] ], individuals[ I[p1] ], individuals[ I[p2] ], best_individual);
                    break;
                case 3:

                    crossover.position_based_crossover(new_generation[ I[i] ], individuals[ I[p1] ], individuals[ I[p2] ]);
                    break;
                case 4:

                    crossover.uniform_crossover(new_generation[ I[i] ], individuals[ I[p1] ], individuals[ I[p2] ]);
                    break;
                }

                //mutazione genetica solo con un certo rateo
                if (random_mut(mt) == 0)
                {
                    switch (gdata.mutator)
                    {
                    case 0:

                        mutation.swap2(new_generation[ I[i] ]);
                        break;
                    case 1:

                        mutation.swap3(new_generation[ I[i] ]);
                        break;
                    case 2:

                        mutation.scramble(new_generation[ I[i] ]);

                        break;
                    case 3:

                        mutation.inversion(new_generation[ I[i] ]);
                        break;
                    case 4:

                        mutation.insertion(new_generation[ I[i] ]);
                        break;
                    }
                }

                repair.repair(new_generation[ I[i] ]);

                improvement.improvement_algorithm(new_generation[ I[i] ]);

                costs.calculate_cost( new_generation[ I[i] ] );
                new_mean_cost += new_generation[i].get_cost();

                if (new_generation[ I[i] ].get_cost() < best_cost)
                {
                    best_cost = new_generation[ I[i] ].get_cost();
                    best_individual = new_generation[ I[i] ];
                    best_individual2 = new_generation[ I[i] ];

                    if ((unsigned)best_individual.get_cost() == instance.known_solution)
                    {
                        printStats();
                        cout << "known solution ";
                        //best_individual.print_tour();
                        return best_individual.get_cost();
                    }
                    //std::cout << "Child improved: " << best_cost << "\n";
                }

                new_mean_cost /= (double)popsize;
                mean_cost = new_mean_cost;
            }

            //cout << "       iteration " << g << ": entire population processed\n";

            std::vector<Individual> *const temp = individuals_ptr;
            individuals_ptr = new_generation_ptr;
            new_generation_ptr = temp;

            
            std::sort(I.begin(), I.end(), sorts);
            sorts._swap();
            
            //cout << "       iteration " << g << ": entire population sorted\n";
            ++g;
        }

        cost = best_individual.get_cost();
        if (cost < global_best)
        {
            best_individual2 = best_individual;
            global_best = cost;    

            if(cost == instance.known_solution)
            {            
                printStats();
                cout << "known solution ";
                //best_individual.print_tour();
                return cost;
            }
        } 
        
    }

    printStats();
    best_individual2.print_tour();
    return cost;
}

#endif