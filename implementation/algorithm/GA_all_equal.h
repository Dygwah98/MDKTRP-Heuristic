#ifndef GA_EQ_H
#define GA_EQ_H

#include"../individual/misc.h"

//parametri dell'algoritmo
struct AllEqualGeneticAlgorithmData {

    //numero di esecuzioni dell'algoritmo
    static constexpr unsigned tries = 1;
    //cardinalità della popolazione
    static constexpr unsigned population_size = 25;
    // numero di mutators e crossover (può essere incrementato di uno per aggiungere la possibilità che non venga eseguito)
    static constexpr unsigned mutator_choice = 3;
    static constexpr unsigned crossover_choice = 2;
#ifdef TIMELIMIT
    //timelimit minimo (a seconda dell'istanza, può essere fino a 20 volte più grande)
    static constexpr double timelimit = 30.0;
#endif
    //numero massimo di iterazioni
    static constexpr unsigned max_evaluations_GA = 10000; 
    //probabilità di mutazione, va letto: 1/mut_rate
    static constexpr unsigned mut_rate = 6; 
    //numero di iterazioni senza miglioramenti prima di attivare random retart/hypermutation
    static constexpr unsigned mut_update_window = 100;
    //frazione di individui mantenuti durante l'esecuzione
    static constexpr double elite_ratio = 0.4;

    static void print() {
        cout << "   Tries: " << AllEqualGeneticAlgorithmData::tries << endl;
        cout << "   Population size: " << AllEqualGeneticAlgorithmData::population_size << endl;
        cout << "   Different mutators: " << AllEqualGeneticAlgorithmData::mutator_choice << endl;
        cout << "   Different crossovers: " << AllEqualGeneticAlgorithmData::crossover_choice << endl;
#ifdef TIMELIMIT
        cout << "   Timelimit: " << AllEqualGeneticAlgorithmData::timelimit << endl;
#endif
        cout << "   Max iterations: " << AllEqualGeneticAlgorithmData::max_evaluations_GA << endl;
        cout << "   Mutation rate: " << AllEqualGeneticAlgorithmData::mut_rate << endl;
        cout << "   Iterations to restart: " << AllEqualGeneticAlgorithmData::mut_update_window << endl;
        cout << "   Elite ratio: " << AllEqualGeneticAlgorithmData::elite_ratio << endl;
    }
};

double AllEqualGeneticAlgorithm(const Test& instance, const Individual& ind) {
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
    
    const unsigned max_tries = GeneticAlgorithmData::tries;

    //vettore d'appoggio, contenente gli indici degli individui
    std::vector<unsigned> D;
    for(unsigned i = 0; i < popsize; ++i)
        D.push_back(i);

    const unsigned mut_update_window = GeneticAlgorithmData::mut_update_window;
    
    const unsigned max_repeated = (mut_rate * mut_update_window)*2;

#ifndef TIMELIMIT
    double best_time = 7200;
#else
    double best_time = timelimit;
#endif
    unsigned best_iteration = max_evaluations;

    for (unsigned tries_i = 1; tries_i <= max_tries; ++tries_i)
    {
        unsigned repeated = 0;

        std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

        double best_cost = std::numeric_limits<double>::max();

        //vettore che contiene la popolazione
        std::vector<Individual> individuals_original;
        individuals_original.assign(popsize, ind);
        
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
                std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

                best_time = (double)std::chrono::duration_cast
                            <std::chrono::nanoseconds>(end_time - start_time)
                            .count() / (double)1000000000;
                
                best_iteration = tries_i;

                if (cost == instance.known_solution)
                {
#ifdef PRINT
                    printStats();
                    best_individual.print_tour();
#endif
                    cout << best_time << " " << best_iteration << " ";
                    return cost;                    
                }
            }

        }
    

        for(unsigned pos = 0; pos < popsize; ++pos) {
            auto& i = individuals_original[pos];
#ifdef PRINT
            metrics.measure_time(i, &Individual::calculate_diversity_ratio, pos, D, individuals_original);
#else
            i.calculate_diversity_ratio(pos, D, individuals_original);
#endif
            i.set_normalized(best_cost);
        }

        const double elite_ratio = GeneticAlgorithmData::elite_ratio;
        const unsigned s  = elite_ratio * (double)popsize;
        
        std::uniform_int_distribution<unsigned> random_parent(1, popsize/2 - 2);
        std::uniform_int_distribution<unsigned> random_mut(1, mut_rate);
        std::uniform_int_distribution<unsigned> random_choice(0, 1);
        std::uniform_int_distribution<unsigned> random_cross(0, AllEqualGeneticAlgorithmData::crossover_choice);
        std::uniform_int_distribution<unsigned> random_mutator(0, AllEqualGeneticAlgorithmData::mutator_choice);

        std::vector<Individual> new_generation_original;
        new_generation_original.assign(popsize, ind);

        //keySort sorts(individuals_original, new_generation_original);

        //si occupa di ordinare il vector di indici in base al rank degli individui
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

        double time_elapsed = 0;

        while (g < max_g && time_elapsed < 7200) {

            ++repeated;
        
            const std::vector<Individual> &individuals = *individuals_ptr;
            std::vector<Individual> &new_generation = *new_generation_ptr;
            
            unsigned i = 0;
            //restart procedure
            if(repeated % mut_update_window == 0) {
                
                if(mut_rate-1 >= 1) {
                    --mut_rate;
                    random_mut = std::uniform_int_distribution<unsigned>(1, mut_rate);
                }

                for (; i < s; ++i) {
                    new_generation[ D[i] ] = individuals[ D[i] ];
                }

                for(; i < popsize; ++i) {
#ifdef PRINT
                    initialize.measure_time(new_generation[ D[i] ], &Individual::random_restart);
                    costs.measure_time(new_generation[ D[i] ], &Individual::calculate_cost);
#else
                    new_generation[ D[i] ].random_restart();
                    new_generation[ D[i] ].calculate_cost();
#endif

                    if (best_cost - new_generation[ D[i] ].get_cost() > std::numeric_limits<double>::epsilon())
                    {
                        best_cost = new_generation[ D[i] ].get_cost();
                        best_individual = new_generation[ D[i] ];

                        std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

                        best_time = (double)std::chrono::duration_cast
                                    <std::chrono::nanoseconds>(end_time - start_time)
                                    .count() / (double)1000000000;
                        
                        best_iteration = tries_i*g;

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
                            cout << best_time << " " << best_iteration << " ";
                            return best_individual.get_cost();
                        }
                    }
                }

                if(repeated == max_repeated) {
                    repeated = 0;
                }

            } else {

                //elitismo
                for (; i < s; ++i) {
                    
                    new_generation[ D[i] ] = individuals[ D[i] ];
                }

                //dobbiamo inserire la popolazione rimanente
                for (; i < popsize; ++i) {

                    //scelta dei parent per il crossover
                    //il primo è in range [0, popsize/2]
                    p1 = random_parent(mt);
                    //il secondo è in range [0, popsize-1]
                    std::uniform_int_distribution<unsigned> left(0, p1-1);
                    std::uniform_int_distribution<unsigned> right(p1 + 1, popsize - 1);
                    if(random_choice(mt))
                        p2 = left(mt);
                    else
                        p2 = right(mt);
        
#ifdef PRINT
                    switch (random_cross(mt))
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

                    costs.measure_time( new_generation[ D[i] ], &Individual::calculate_cost );
                    costs.measure_time( spare_son, &Individual::calculate_cost );

                    if(new_generation[ D[i] ].is_feasible() && (best_cost - new_generation[ D[i] ].get_cost() > std::numeric_limits<double>::epsilon())
                    || spare_son.is_feasible() && (best_cost - spare_son.get_cost() > std::numeric_limits<double>::epsilon())) {

                        if(new_generation[ D[i] ].get_cost() - spare_son.get_cost() > std::numeric_limits<double>::epsilon())
                            new_generation[ D[i] ] = spare_son;

                        best_cost = new_generation[ D[i] ].get_cost();
                        best_individual = new_generation[ D[i] ];

                        std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

                        best_time = (double)std::chrono::duration_cast
                                    <std::chrono::nanoseconds>(end_time - start_time)
                                    .count() / (double)1000000000;
                        
                        best_iteration = tries_i*g;

                        repeated = 0;
                        if(mut_rate < original_mut_rate) {
                            ++mut_rate;
                            random_mut = std::uniform_int_distribution<unsigned>(1, mut_rate);
                        }

                        if ((unsigned)best_individual.get_cost() == instance.known_solution)
                        {

                            printStats();
                            best_individual.print_tour(); 
                            cout << best_time << " " << best_iteration << " " << endl;
                            return best_individual.get_cost();
                        }
                    } else {
                        
                        //mutazione genetica solo con un certo rateo
                        if (random_mut(mt) == 1)
                        {
                            switch (random_mutator(mt))
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

                        costs.measure_time( new_generation[ D[i] ], &Individual::calculate_cost );
                        costs.measure_time( spare_son, &Individual::calculate_cost );

                        if(new_generation[ D[i] ].is_feasible() && (best_cost - new_generation[ D[i] ].get_cost() > std::numeric_limits<double>::epsilon())
                        || spare_son.is_feasible() && (best_cost - spare_son.get_cost() > std::numeric_limits<double>::epsilon())) {

                            if(new_generation[ D[i] ].get_cost() - spare_son.get_cost() > std::numeric_limits<double>::epsilon())
                                new_generation[ D[i] ] = spare_son;

                            best_cost = new_generation[ D[i] ].get_cost();
                            best_individual = new_generation[ D[i] ];

                            std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

                            best_time = (double)std::chrono::duration_cast
                                        <std::chrono::nanoseconds>(end_time - start_time)
                                        .count() / (double)1000000000;
                            
                            best_iteration = tries_i*g;

                            repeated = 0;
                            if(mut_rate < original_mut_rate) {
                                ++mut_rate;
                                random_mut = std::uniform_int_distribution<unsigned>(1, mut_rate);
                            }

                            if ((unsigned)best_individual.get_cost() == instance.known_solution)
                            {

                                printStats();
                                best_individual.print_tour();
                                cout << best_time << " " << best_iteration << " ";
                                return best_individual.get_cost();
                            }
                        } else {
                            repair.measure_time(new_generation[ D[i] ], &Individual::repair);
                            repair.measure_time(spare_son, &Individual::repair);

                            improvement.measure_time(new_generation[ D[i] ], &Individual::local_search);
                            improvement.measure_time(spare_son, &Individual::local_search);

                            costs.measure_time( new_generation[ D[i] ], &Individual::calculate_cost );
                            costs.measure_time( spare_son, &Individual::calculate_cost );

                            if(new_generation[ D[i] ].is_feasible() && (best_cost - new_generation[ D[i] ].get_cost() > std::numeric_limits<double>::epsilon())
                            || spare_son.is_feasible() && (best_cost - spare_son.get_cost() > std::numeric_limits<double>::epsilon()))
                            {
                                //conserviamo solo il figlio migliore
                                if(new_generation[ D[i] ].get_cost() - spare_son.get_cost() > std::numeric_limits<double>::epsilon()) {
                                    new_generation[ D[i] ] = spare_son;
                                }

                                best_cost = new_generation[ D[i] ].get_cost();
                                best_individual = new_generation[ D[i] ];

                                std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

                                best_time = (double)std::chrono::duration_cast
                                            <std::chrono::nanoseconds>(end_time - start_time)
                                            .count() / (double)1000000000;
                                
                                best_iteration = tries_i*g;

                                repeated = 0;
                                if(mut_rate < original_mut_rate) {
                                    ++mut_rate;
                                    random_mut = std::uniform_int_distribution<unsigned>(1, mut_rate);
                                }

                                if ((unsigned)best_individual.get_cost() == instance.known_solution)
                                {
                                    printStats();
                                    best_individual.print_tour(); 
                                    cout << best_time << " " << best_iteration << " ";
                                    return best_individual.get_cost();
                                }
                            }
                        }                        
                    }                   
#else
                    switch (random_cross(mt))
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

                    new_generation[ D[i] ].calculate_cost();
                    spare_son.calculate_cost();

                    if(new_generation[ D[i] ].is_feasible() && (best_cost - new_generation[ D[i] ].get_cost() > std::numeric_limits<double>::epsilon())
                    || spare_son.is_feasible() && (best_cost - spare_son.get_cost() > std::numeric_limits<double>::epsilon())) {

                        if(new_generation[ D[i] ].get_cost() - spare_son.get_cost() > std::numeric_limits<double>::epsilon())
                            new_generation[ D[i] ] = spare_son;

                        best_cost = new_generation[ D[i] ].get_cost();
                        best_individual = new_generation[ D[i] ];

                        std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

                        best_time = (double)std::chrono::duration_cast
                                    <std::chrono::nanoseconds>(end_time - start_time)
                                    .count() / (double)1000000000;
                        
                        best_iteration = tries_i*g;

                        repeated = 0;
                        if(mut_rate < original_mut_rate) {
                            ++mut_rate;
                            random_mut = std::uniform_int_distribution<unsigned>(1, mut_rate);
                        }

                        if ((unsigned)best_individual.get_cost() == instance.known_solution)
                        {
                            cout << best_time << " " << best_iteration << " ";
                            return best_individual.get_cost();
                        }
                    } 
                    else 
                    {
                        //mutazione genetica solo con un certo rateo
                        if (random_mut(mt) == 1)
                        {
                            switch (random_mutator(mt))
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
                        
                        new_generation[ D[i] ].calculate_cost();
                        spare_son.calculate_cost();

                        if(new_generation[ D[i] ].is_feasible() && (best_cost - new_generation[ D[i] ].get_cost() > std::numeric_limits<double>::epsilon())
                        || spare_son.is_feasible() && (best_cost - spare_son.get_cost() > std::numeric_limits<double>::epsilon())) {

                            if(new_generation[ D[i] ].get_cost() - spare_son.get_cost() > std::numeric_limits<double>::epsilon())
                                new_generation[ D[i] ] = spare_son;

                            best_cost = new_generation[ D[i] ].get_cost();
                            best_individual = new_generation[ D[i] ];

                            std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

                            best_time = (double)std::chrono::duration_cast
                                        <std::chrono::nanoseconds>(end_time - start_time)
                                        .count() / (double)1000000000;
                            
                            best_iteration = tries_i*g;

                            repeated = 0;
                            if(mut_rate < original_mut_rate) {
                                ++mut_rate;
                                random_mut = std::uniform_int_distribution<unsigned>(1, mut_rate);
                            }

                            if ((unsigned)best_individual.get_cost() == instance.known_solution)
                            {
                                cout << best_time << " " << best_iteration << " ";
                                return best_individual.get_cost();
                            }
                        } else {
                            new_generation[ D[i] ].repair();
                            spare_son.repair();

                            new_generation[ D[i] ].local_search();
                            spare_son.local_search();

                            new_generation[ D[i] ].calculate_cost();
                            spare_son.calculate_cost();

                            //conserviamo solo il figlio migliore
                            if(new_generation[ D[i] ].get_cost() - spare_son.get_cost() > std::numeric_limits<double>::epsilon()) {
                                new_generation[ D[i] ] = spare_son;
                            }

                            if (new_generation[ D[i] ].is_feasible() 
                                && (best_cost - new_generation[ D[i] ].get_cost() > std::numeric_limits<double>::epsilon()))
                            {
                                best_cost = new_generation[ D[i] ].get_cost();
                                best_individual = new_generation[ D[i] ];

                                std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

                                best_time = (double)std::chrono::duration_cast
                                            <std::chrono::nanoseconds>(end_time - start_time)
                                            .count() / (double)1000000000;
                                
                                best_iteration = tries_i*g;

                                repeated = 0;
                                if(mut_rate < original_mut_rate) {
                                    ++mut_rate;
                                    random_mut = std::uniform_int_distribution<unsigned>(1, mut_rate);
                                }

                                if ((unsigned)best_individual.get_cost() == instance.known_solution)
                                {
                                    cout << best_time << " " << best_iteration << " ";
                                    return best_individual.get_cost();
                                }
                            }
                        }  
                    }
#endif
                }
                
                unsigned j = 0;
                for(; j < popsize; ++j) {
#ifdef PRINT
                    metrics.measure_time(new_generation[ D[j] ], &Individual::calculate_diversity_ratio, j, D, new_generation);
#else
                    new_generation[ D[j] ].calculate_diversity_ratio(j, D, new_generation);
#endif
                    new_generation[ D[j] ].set_normalized(best_cost);
                }

            }

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
                cout << best_time << " " << best_iteration << " ";
                return cost;
            }
        } 
        
    }
#ifdef PRINT
    
    printStats();
    best_individual.print_tour();
#endif
    cout << best_time << " " << best_iteration << " ";
    return cost;
}

#endif