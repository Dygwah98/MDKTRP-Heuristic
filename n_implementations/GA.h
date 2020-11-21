#ifndef GA_H
#define GA_H

#include"individual/misc.h"

struct GeneticAlgorithmData {

    static constexpr unsigned tries = 1;
    static constexpr unsigned population_size = 25;
    static constexpr unsigned mutator = SCRAMBLE;
    static constexpr unsigned crossover = TWO_POINT;
#ifdef TIMELIMIT
    static constexpr double timelimit = 30.0;
#endif
    static constexpr unsigned max_evaluations_GA = 10000; 
    //va letto: 1/mut_rate prob di mutazione
    static constexpr unsigned mut_rate = 6; 
    //numero di iterazioni senza miglioramenti prima di attivare random retart/hypermutation
    static constexpr unsigned mut_update_window = 100;
    static constexpr double elite_ratio = 0.4;
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
    Timer& metrics = printStats.metrics;
    Timer& sort_time = printStats.sort;
#endif
    //schema genetico
#ifdef TIMELIMIT    
    double timelimit = 
        (GeneticAlgorithmData::timelimit / GeneticAlgorithmData::tries) 
      * ( (instance.factor_valuations) );
    
    if(timelimit > 7200)
        timelimit = 7200;
#endif
    const unsigned max_evaluations = GeneticAlgorithmData::max_evaluations_GA * instance.factor_valuations;
    const unsigned max_g = (max_evaluations / GeneticAlgorithmData::tries) / GeneticAlgorithmData::population_size;

    const unsigned original_mut_rate = GeneticAlgorithmData::mut_rate;
    unsigned mut_rate = GeneticAlgorithmData::mut_rate;

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

    std::vector<unsigned> D;
    for(unsigned i = 0; i < popsize; ++i)
        D.push_back(i);

    const unsigned mut_update_window = GeneticAlgorithmData::mut_update_window;
    
    const unsigned max_repeated = (mut_rate * mut_update_window)*2;

    for (unsigned tries_i = 0; tries_i < max_tries; ++tries_i)
    {
        unsigned repeated = 0;

        std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

        double best_cost = std::numeric_limits<double>::max();

        //vettore che contiene la popolazione
        std::vector<Individual> individuals_original;
        individuals_original.assign(popsize, ind);
        
        double mean_cost = 0;
        //inizializzazione della popolazione
        for (unsigned pos = 0; pos < popsize; ++pos)
        {
            auto& i = individuals_original[pos];
#ifdef PRINT
            initialize.measure_time(i, &Individual::random_initialize);
            costs.measure_time(i, &Individual::calculate_cost);
#else
            i.random_initialize();
            i.calculate_cost();
#endif
            cost = i.get_cost(); 
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
                    //cout << "known solution ";
                    return cost;                    
                }
            }

            mean_cost += cost;
        }
        //cout << "       population initialized\n";

        mean_cost /= (double)popsize;

        for(unsigned pos = 0; pos < popsize; ++pos) {
            auto& i = individuals_original[pos];
#ifdef PRINT
            metrics.measure_time(i, &Individual::calculate_diversity_ratio, pos, individuals_original);
#else
            i.calculate_diversity_ratio(pos, individuals_original);
#endif
            i.set_normalized(best_cost);
        }

        const double elite_ratio = GeneticAlgorithmData::elite_ratio;
        const unsigned s  = elite_ratio * (double)popsize;
        
        std::uniform_int_distribution<unsigned> random_parent(1, popsize/2 - 2);
        std::uniform_int_distribution<unsigned> random_mut(1, mut_rate);
        std::uniform_int_distribution<unsigned> random_choice(0, 1);

        std::vector<Individual> new_generation_original;
        new_generation_original.assign(popsize, ind);

        keySort sorts(individuals_original, new_generation_original);
        keyDiversitySort sorts2(individuals_original, new_generation_original, 1.0 - elite_ratio);

#ifdef PRINT
        sort_time.measure_time(D, sorts2);
#else
        std::sort(D.begin(), D.end(), sorts2);
        sorts2._swap();
#endif
        unsigned g = 0;
        
        unsigned p1;
        unsigned p2; 

        std::vector<Individual> *individuals_ptr = &individuals_original;
        std::vector<Individual> *new_generation_ptr = &new_generation_original; 

        //cout << "       starting evaluations\n";
        double time_elapsed = 0;

        while (g < max_g && time_elapsed < 7200) {

            double new_mean_cost = 0;
            ++repeated;
            //cout<<"G: "<<g << "| ";

            const std::vector<Individual> &individuals = *individuals_ptr;
            std::vector<Individual> &new_generation = *new_generation_ptr;

            
            unsigned i = 0;
            if(repeated % mut_update_window == 0) {
                
                if(mut_rate-1 >= 1) {
                    --mut_rate;
                    random_mut = std::uniform_int_distribution<unsigned>(1, mut_rate);
                }

                //cout << mut_update_window << " elite ";
                for (; i < s; ++i) {
                    new_generation[ D[i] ] = individuals[ D[i] ];
                    new_mean_cost += new_generation[ D[i] ].get_cost();
                }

                for(; i < popsize; ++i) {
#ifdef PRINT
                    initialize.measure_time(new_generation[ D[i] ], &Individual::random_restart);
                    costs.measure_time(new_generation[ D[i] ], &Individual::calculate_cost);
#else
                    new_generation[ D[i] ].random_restart();
                    new_generation[ D[i] ].calculate_cost();
#endif
                    new_mean_cost += new_generation[ D[i] ].get_cost();

                    if (best_cost - new_generation[ D[i] ].get_cost() > std::numeric_limits<double>::epsilon())
                    {
                        best_cost = new_generation[ D[i] ].get_cost();
                        best_individual = new_generation[ D[i] ];
                        repeated = 0;
                        if(mut_rate < 6) {
                            ++mut_rate;
                            random_mut = std::uniform_int_distribution<unsigned>(1, mut_rate);
                        }

                        if ((unsigned)best_individual.get_cost() == instance.known_solution)
                        {
#ifdef PRINT
                            printStats();
                            best_individual.print_tour();
#endif
                            //cout << "known solution ";
                            return best_individual.get_cost();
                        }
                        //std::cout << "Child improved: " << best_cost << "\n";
                    }
                }

                if(repeated == max_repeated) {
                    repeated = 0;
                }

            } else {

                //elitismo al 10%
                for (; i < s; ++i) {
                    new_generation[ D[i] ] = individuals[ D[i] ];
                    new_mean_cost += new_generation[ D[i] ].get_cost();
                }
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

                        spare_son = crossover.measure_time(new_generation[ D[i] ], &Individual::one_point_cross_over, 
                                                individuals[ D[p1] ], individuals[ D[p2] ]);
                        break;
                    case 1:

                        spare_son = crossover.measure_time(new_generation[ D[i] ], &Individual::two_point_cross_over,
                                                individuals[ D[p1] ], individuals[ D[p2] ]);
                        break;
                    case 2:

                        spare_son = crossover.measure_time(new_generation[ D[i] ], &Individual::uniform_cross_over, 
                                                individuals[ D[p1] ], individuals[ D[p2] ]);
                        break;
                    default:
                        break;
                    }

                    //mutazione genetica solo con un certo rateo
                    if (random_mut(mt) == 1)
                    {
                        switch (GeneticAlgorithmData::mutator)
                        {
                        case 0:

                            mutation.measure_time(new_generation[ D[i] ], &Individual::swap3);
                            mutation.measure_time(spare_son, &Individual::swap3);
                            break;
                        case 1:

                            mutation.measure_time(new_generation[ D[i] ], &Individual::swap5);
                            mutation.measure_time(spare_son, &Individual::swap5);
                            break;
                        case 2:

                            mutation.measure_time(new_generation[ D[i] ], &Individual::scramble);
                            mutation.measure_time(spare_son, &Individual::scramble);
                            break;
                        case 3:

                            mutation.measure_time(new_generation[ D[i] ], &Individual::inversion);
                            mutation.measure_time(spare_son, &Individual::inversion);
                            break;
                        default:
                            break;
                        }
                    }

                    repair.measure_time(new_generation[ D[i] ], &Individual::repair);
                    repair.measure_time(spare_son, &Individual::repair);

                    improvement.measure_time(new_generation[ D[i] ], &Individual::local_search);
                    improvement.measure_time(spare_son, &Individual::local_search);

                    costs.measure_time( new_generation[ D[i] ], &Individual::calculate_cost );
                    costs.measure_time( spare_son, &Individual::calculate_cost );
        #else
                    switch (GeneticAlgorithmData::crossover)
                    {
                    case 0:
                        spare_son = new_generation[ D[i] ].one_point_cross_over(individuals[ D[p1] ], individuals[ D[p2] ]);
                        break;
                    case 1:
                        spare_son = new_generation[ D[i] ].two_point_cross_over(individuals[ D[p1] ], individuals[ D[p2] ]);
                        break;
                    case 2:
                        spare_son =  new_generation[ D[i] ].uniform_cross_over(individuals[ D[p1] ], individuals[ D[p2] ]);
                        break;
                    default:
                        break;
                    }

                    //mutazione genetica solo con un certo rateo
                    if (random_mut(mt) == 1)
                    {
                        switch (GeneticAlgorithmData::mutator)
                        {
                        case 0:

                            new_generation[ D[i] ].swap3();
                            spare_son.swap3();
                            break;
                        case 1:

                            new_generation[ D[i] ].swap5();
                            spare_son.swap5();
                            break;
                        case 2:
                            
                            new_generation[ D[i] ].scramble();
                            spare_son.scramble();
                            break;
                        case 3:

                            new_generation[ D[i] ].inversion();
                            spare_son.inversion();
                            break;
                        default:
                            break;
                        }
                    }
                    
                    new_generation[ D[i] ].repair();
                    spare_son.repair();

                    new_generation[ D[i] ].local_search();
                    spare_son.local_search();

                    new_generation[ D[i] ].calculate_cost();
                    spare_son.calculate_cost();
        #endif

                    if(new_generation[ D[i] ].get_cost() - spare_son.get_cost() > std::numeric_limits<double>::epsilon()) {
                        new_generation[ D[i] ] = spare_son;
                    }

                    new_mean_cost += new_generation[ D[i] ].get_cost();

                    if (new_generation[ D[i] ].is_feasible() 
                        && (best_cost - new_generation[ D[i] ].get_cost() > std::numeric_limits<double>::epsilon()))
                    {
                        best_cost = new_generation[ D[i] ].get_cost();
                        best_individual = new_generation[ D[i] ];
                        repeated = 0;
                        if(mut_rate < original_mut_rate) {
                            ++mut_rate;
                            random_mut = std::uniform_int_distribution<unsigned>(1, mut_rate);
                        }

                        if ((unsigned)best_individual.get_cost() == instance.known_solution)
                        {
#ifdef PRINT
                            printStats();
                            best_individual.print_tour();
#endif
                            //cout << "known solution ";
                            return best_individual.get_cost();
                        }
                        //std::cout << "Child improved: " << best_cost << "\n";
                    }
                }
                
                new_mean_cost /= (double)popsize;
                mean_cost = new_mean_cost;
                
                unsigned j = 0;

                for(; j < popsize; ++j) {
#ifdef PRINT
                    metrics.measure_time(new_generation[ D[j] ], &Individual::calculate_diversity_ratio, D[j], new_generation);
#else
                    new_generation[ D[j] ].calculate_diversity_ratio(D[j], new_generation);
#endif
                    new_generation[ D[j] ].set_normalized(best_cost);
                }

            }
            //cout << "\n";
            //cout << "       iteration " << g << ": elite population processed\n";

            //cout << "       iteration " << g << ": entire population processed\n";

            std::vector<Individual> *const temp = individuals_ptr;
            individuals_ptr = new_generation_ptr;
            new_generation_ptr = temp;
#ifdef PRINT
            sort_time.measure_time(D, sorts2);
#else
            std::sort(D.begin(), D.end(), sorts2);
            sorts2._swap();
#endif           
            std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

            time_elapsed = (double)std::chrono::duration_cast
                        <std::chrono::nanoseconds>(end_time - start_time)
                        .count() / (double)1000000000;
#ifdef TIMELIMIT             
            //cout << time_elapsed << "\n";
            if(time_elapsed > timelimit)
                break;
#endif

            ++g;

        }

        cost = best_individual.get_cost();
        if (global_best - cost > std::numeric_limits<double>::epsilon())
        {
            global_best = cost;    

            if(cost == instance.known_solution)
            {  
#ifdef PRINT          
                printStats();
                best_individual.print_tour();
#endif
                //cout << "known solution ";
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