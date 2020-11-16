#ifndef GA_H
#define GA_H

#include"individual/misc.h"

struct GeneticAlgorithmData {

    static constexpr unsigned tries = 1;
    static constexpr unsigned population_size = 500;
    static constexpr unsigned mutator = SCRAMBLE;
    static constexpr unsigned crossover = TWO_POINT;
#ifdef TIMELIMIT
    static constexpr double timelimit = 10.0;
#else 
    static constexpr unsigned max_evaluations_GA = 40000000 / 3; 
#endif
    static constexpr unsigned mut_rate = 2;
};

double GeneticAlgorithm(const Test& instance, const Individual& ind) {

    //variabili per misurare le prestazioni
#ifdef PRINT    
    Analyzer printStats; 
    Timer& initialize = printStats.initialize;
    Timer& crossover = printStats.crossover;
    Timer& mutation = printStats.mutation;
    Timer& improvement = printStats.improvement;
    Timer& repair = printStats.repair;
    Timer& costs = printStats.costs;
#endif
    //schema genetico
#ifdef TIMELIMIT    
    double timelimit = 
        (GeneticAlgorithmData::timelimit / GeneticAlgorithmData::tries) 
      * ( (instance.factor_valuations) );
    
    if(timelimit > 7200)
        timelimit = 7200;
#else 
    const unsigned max_evaluations = GeneticAlgorithmData::max_evaluations_GA * instance.factor_valuations;
    const unsigned max_g = (max_evaluations / GeneticAlgorithmData::tries) / GeneticAlgorithmData::population_size;
#endif


    const unsigned popsize = GeneticAlgorithmData::population_size;
    
    double cost = std::numeric_limits<double>::max();
    double global_best = std::numeric_limits<double>::max();

    Individual best_individual(ind);
    Individual spare_son(ind);
#ifdef PRINT    
    initialize.measure_time(best_individual, &Individual::random_initialize);
#else
    best_individual.random_initialize();
#endif
    //cout << "       parameters set\n";
    const unsigned max_tries = GeneticAlgorithmData::tries;

    std::vector<unsigned> I;
    for(unsigned i = 0; i < popsize; ++i)
        I.push_back(i);

    for (unsigned tries_i = 0; tries_i < max_tries; ++tries_i)
    {
        unsigned repeated = 0;
#ifdef TIMELIMIT
        std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
#endif
        double best_cost = std::numeric_limits<double>::max();

        //vettore che contiene la popolazione
        std::vector<Individual> individuals_original;
        individuals_original.assign(popsize, ind);
        
        double mean_cost = 0;
        //inizializzazione della popolazione
        for (Individual &i : individuals_original)
        {
#ifdef PRINT
            initialize.measure_time(i, &Individual::random_initialize);

            costs.measure_time(i, &Individual::calculate_cost);
#else
            i.random_initialize();
            i.calculate_cost();
#endif
            double cost = i.get_cost(); 
            if (cost < best_cost)
            {
                best_cost = cost;
                best_individual = i;

                if (cost == instance.known_solution)
                {
#ifdef PRINT
                    printStats();
                    best_individual.print_tour();
#endif
                    cout << "known solution ";
                    return cost;                    
                }
            }

            mean_cost += cost;
        }
        //cout << "       population initialized\n";

        mean_cost /= (double)popsize;

        const unsigned s  = (10 * popsize) / 100;
        
        std::uniform_int_distribution<unsigned> random_parent(1, popsize/2 - 2);
        std::uniform_int_distribution<unsigned> random_mut(0, GeneticAlgorithmData::mut_rate);
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
        double time_elapsed = 0;
#ifdef TIMELIMIT
        while(true) {
#else 
        while (g < max_g) {
#endif
            double new_mean_cost = 0;
            ++repeated;
            //cout<<"G: "<<g<<"\n";

            const std::vector<Individual> &individuals = *individuals_ptr;
            std::vector<Individual> &new_generation = *new_generation_ptr;

            //elitismo al 10%
            unsigned i = 0;
            if(repeated == 500) {
                for (; i < s/2; ++i) {
                    new_generation[ I[i] ] = individuals[ I[i] ];
#ifdef PRINT
                    initialize.measure_time(new_generation[ I[i] ], &Individual::random_restart);
                    repair.measure_time(new_generation[ I[i] ], &Individual::repair);
#else
                    new_generation[ I[i] ].random_restart();
                    new_generation[ I[i] ].repair();
#endif
                    new_mean_cost += new_generation[i].get_cost();
                }

                for (; i < s; ++i) {
                    new_generation[ I[i] ] = individuals[ I[i] ];
                    new_mean_cost += new_generation[i].get_cost();
                }

                repeated = 0;

            } else {
                
                for (; i < s; ++i) {
                    new_generation[ I[i] ] = individuals[ I[i] ];
                    new_mean_cost += new_generation[i].get_cost();
                }
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
#ifdef PRINT
                switch (GeneticAlgorithmData::crossover)
                {
                case 0:

                    spare_son = crossover.measure_time(new_generation[ I[i] ], &Individual::one_point_cross_over, 
                                            individuals[ I[p1] ], individuals[ I[p2] ]);
                    break;
                case 1:

                    spare_son = crossover.measure_time(new_generation[ I[i] ], &Individual::two_point_cross_over,
                                            individuals[ I[p1] ], individuals[ I[p2] ]);
                    break;
                case 2:

                    spare_son = crossover.measure_time(new_generation[ I[i] ], &Individual::best_order_cross_over,
                                            individuals[ I[p1] ], individuals[ I[p2] ], best_individual);
                    break;
                case 3:

                    spare_son = crossover.measure_time(new_generation[ I[i] ], &Individual::position_based_cross_over,
                                            individuals[ I[p1] ], individuals[ I[p2] ]);
                    break;
                case 4:

                    spare_son = crossover.measure_time(new_generation[ I[i] ], &Individual::uniform_cross_over, 
                                            individuals[ I[p1] ], individuals[ I[p2] ]);
                    break;
                }

                repair.measure_time(new_generation[ I[i] ], &Individual::repair);
                repair.measure_time(spare_son, &Individual::repair);

                costs.measure_time(new_generation[ I[i] ], &Individual::calculate_cost);
                costs.measure_time(spare_son, &Individual::calculate_cost);

                if(spare_son.get_cost() < new_generation[ I[i] ].get_cost())
                    new_generation[ I[i] ] = spare_son;

                //mutazione genetica solo con un certo rateo
                if (random_mut(mt) == 0)
                {
                    switch (GeneticAlgorithmData::mutator)
                    {
                    case 0:

                        mutation.measure_time(new_generation[ I[i] ], &Individual::swap2);
                        mutation.measure_time(spare_son, &Individual::swap2);
                        break;
                    case 1:

                        mutation.measure_time(new_generation[ I[i] ], &Individual::swap3);
                        mutation.measure_time(spare_son, &Individual::swap3);
                        break;
                    case 2:

                        mutation.measure_time(new_generation[ I[i] ], &Individual::scramble);
                        mutation.measure_time(spare_son, &Individual::scramble);
                        break;
                    case 3:

                        mutation.measure_time(new_generation[ I[i] ], &Individual::inversion);
                        mutation.measure_time(spare_son, &Individual::inversion);
                        break;
                    case 4:

                        mutation.measure_time(new_generation[ I[i] ], &Individual::insertion);
                        mutation.measure_time(spare_son, &Individual::insertion);
                        break;
                    }
                }

                improvement.measure_time(new_generation[ I[i] ], &Individual::improvement_algorithm);
                improvement.measure_time(spare_son, &Individual::improvement_algorithm);

                costs.measure_time( new_generation[ I[i] ], &Individual::calculate_cost );
                costs.measure_time( spare_son, &Individual::calculate_cost );
#else
                switch (GeneticAlgorithmData::crossover)
                {
                case 0:

                    spare_son = new_generation[ I[i] ].one_point_cross_over(individuals[ I[p1] ], individuals[ I[p2] ]);
                    break;
                case 1:

                    spare_son = new_generation[ I[i] ].two_point_cross_over(individuals[ I[p1] ], individuals[ I[p2] ]);
                    break;
                case 2:

                    spare_son = new_generation[ I[i] ].best_order_cross_over(individuals[ I[p1] ], individuals[ I[p2] ], best_individual);
                    break;
                case 3:

                    spare_son = new_generation[ I[i] ].position_based_cross_over(individuals[ I[p1] ], individuals[ I[p2] ]);
                    break;
                case 4:

                    spare_son =  new_generation[ I[i] ].uniform_cross_over(individuals[ I[p1] ], individuals[ I[p2] ]);
                    break;
                }

                new_generation[ I[i] ].repair();
                spare_son.repair();

                //mutazione genetica solo con un certo rateo
                if (random_mut(mt) == 0)
                {
                    switch (GeneticAlgorithmData::mutator)
                    {
                    case 0:

                        new_generation[ I[i] ].swap2();
                        spare_son.swap2();
                        break;
                    case 1:

                        new_generation[ I[i] ].swap3();
                        spare_son.swap3();
                        break;
                    case 2:
                        
                        new_generation[ I[i] ].scramble();
                        spare_son.scramble();
                        break;
                    case 3:

                        new_generation[ I[i] ].inversion();
                        spare_son.inversion();
                        break;
                    case 4:

                        new_generation[ I[i] ].insertion();
                        spare_son.insertion();
                        break;
                    }
                }

                new_generation[ I[i] ].improvement_algorithm();
                spare_son.improvement_algorithm();

                new_generation[ I[i] ].calculate_cost();
                spare_son.calculate_cost();
#endif

                if(spare_son.get_cost() < new_generation[ I[i] ].get_cost()) {
                    new_generation[ I[i] ] = spare_son;
                }

                new_mean_cost += new_generation[i].get_cost();

                if (new_generation[ I[i] ].get_cost() < best_cost)
                {
                    best_cost = new_generation[ I[i] ].get_cost();
                    best_individual = new_generation[ I[i] ];
                    repeated = 0;

                    if ((unsigned)best_individual.get_cost() == instance.known_solution)
                    {
#ifdef PRINT
                        printStats();
                        best_individual.print_tour();
#endif
                        cout << "known solution ";
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
#ifdef TIMELIMIT            
            std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

            time_elapsed = (double)std::chrono::duration_cast
                <std::chrono::nanoseconds>(end_time - start_time)
                .count() / (double)1000000000;
            
            if(time_elapsed > timelimit)
                break;
#else
            ++g;
#endif
        }

        cost = best_individual.get_cost();
        if (cost < global_best)
        {
            global_best = cost;    

            if(cost == instance.known_solution)
            {  
#ifdef PRINT          
                printStats();
                best_individual.print_tour();
#endif
                cout << "known solution ";
                //best_individual.print_tour();
                return cost;
            }
        } 
        
    }
#ifdef PRINT
    
    printStats();
    best_individual.print_tour();
#endif
    return cost;
}

#endif