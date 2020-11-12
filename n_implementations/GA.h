#ifndef GA_H
#define GA_H

#include<vector>
#include<utility>
#include"test.h"
#include"utils.h"
#include"individual/individual.h"


class GeneticAlgorithmData {

    public:
        const unsigned tries = 1;
        const unsigned population_size = 250;
        const unsigned mutator = SWAP2;
        const unsigned crossover = TWO_POINT;
        const unsigned max_evaluations_GA = 40000000 / 3;
        const unsigned mut_rate = 4;
};

double GeneticAlgorithm(const Test& instance, const Individual& ind) {

    //variabili per misurare le prestazioni
    
    //random_initialize()
    double random_initialize_time = 0;
    double random_initialize_calls = 0;
    std::chrono::steady_clock::time_point random_initialize_begin;
    std::chrono::steady_clock::time_point random_initialize_end;

    //repair()
    double repair_time = 0;
    double repair_calls = 0;
    std::chrono::steady_clock::time_point repair_begin;
    std::chrono::steady_clock::time_point repair_end;

    //improvement_algorithm()
    double improvement_algorithm_time = 0;
    double improvement_algorithm_calls = 0;
    std::chrono::steady_clock::time_point improvement_algorithm_begin;
    std::chrono::steady_clock::time_point improvement_algorithm_end;

    //crossover()
    double crossover_time = 0;
    double crossover_calls = 0;
    std::chrono::steady_clock::time_point crossover_begin;
    std::chrono::steady_clock::time_point crossover_end;

    //mutation()
    double mutation_time = 0;
    double mutation_calls = 0;
    std::chrono::steady_clock::time_point mutation_begin;
    std::chrono::steady_clock::time_point mutation_end;


    const GeneticAlgorithmData gdata;
    
    const unsigned max_evaluations = gdata.max_evaluations_GA * instance.factor_valuations;
    const unsigned max_g = (max_evaluations / gdata.tries) / gdata.population_size;
    
    double cost = std::numeric_limits<double>::max();
    double global_best = std::numeric_limits<double>::max();

    Individual best_individual2(ind);
    random_initialize_begin = std::chrono::steady_clock::now();
    
    best_individual2.random_initialize();
    
    random_initialize_end = std::chrono::steady_clock::now();

    random_initialize_time += (double)std::chrono::duration_cast
            <std::chrono::nanoseconds>(random_initialize_end - random_initialize_begin)
            .count() / (double)1000000000;
    random_initialize_calls += 1;

    //cout << "       parameters set\n";
    for (unsigned tries_i = 0; tries_i < gdata.tries; ++tries_i)
    {
        double best_cost = std::numeric_limits<double>::max();

        //vettore che contiene la popolazione
        std::vector<Individual> individuals_original;
        individuals_original.assign(gdata.population_size, ind);

        Individual best_individual(best_individual2);
        
        //inizializzazione della popolazione
        for (Individual &i : individuals_original)
        {
            random_initialize_begin = std::chrono::steady_clock::now();
            
            i.random_initialize();
            
            random_initialize_end = std::chrono::steady_clock::now();
            random_initialize_time += (double)std::chrono::duration_cast
                <std::chrono::nanoseconds>(random_initialize_end - random_initialize_begin)
                .count() / (double)1000000000;
            random_initialize_calls += 1;

            repair_begin =std::chrono::steady_clock::now();
            i.repair();
            repair_end = std::chrono::steady_clock::now();

            repair_time += (double)std::chrono::duration_cast
                <std::chrono::nanoseconds>(repair_end - repair_begin)
                .count() / (double)1000000000;
            repair_calls +=1;

            i.calculate_cost();

            if (i.get_cost() < best_cost)
            {
                best_cost = i.get_cost();
                best_individual = i;
                best_individual2 = i;

                if ((unsigned)best_individual.get_cost() == instance.known_solution)
                {
                    cout <<"    AUDITING:\n";
                    cout <<"    (function): (total time)\n";
                    cout <<"    random_initalize: " << random_initialize_time << endl;
                    cout <<"    repair: " << repair_time << endl;
                    cout <<"    improvement algorithm: " << improvement_algorithm_time << endl;
                    cout <<"    crossover: "  << crossover_time << endl;
                    cout <<"    mutation: " << mutation_time << endl;
                    cout << "known solution ";
                    //best_individual.print_tour();
                    return best_individual.get_cost();                    
                }
            }

            
        }
        //cout << "       population initialized\n";

        std::sort(individuals_original.begin(), individuals_original.end());

        const unsigned s  = (unsigned)ceil( (10.0 * gdata.population_size) / 100.0);
        const unsigned s3 = gdata.population_size;  

        std::uniform_int_distribution<unsigned> random_parent(1, gdata.population_size / 2);
        std::uniform_int_distribution<unsigned> random_mut(0, gdata.mut_rate);
        std::uniform_int_distribution<unsigned> random_choice(0, 1);

        std::vector<Individual> new_generation_original;
        new_generation_original.assign(gdata.population_size, ind);

        unsigned g = 0;
        
        unsigned p1;
        unsigned p2; 

        std::vector<Individual> *individuals_ptr = &individuals_original;
        std::vector<Individual> *new_generation_ptr = &new_generation_original; 

        //cout << "       starting evaluations\n";
        while (g < max_g)
        {
            //cout<<"G: "<<g<<"\n";

            const std::vector<Individual> &individuals = *individuals_ptr;
            std::vector<Individual> &new_generation = *new_generation_ptr;

            //elitismo al 10%
            unsigned i = 0;
            for (; i < s; ++i) {
                new_generation[i] = individuals[i];

                improvement_algorithm_begin = std::chrono::steady_clock::now();

                new_generation[i].improvement_algorithm();

                improvement_algorithm_end = std::chrono::steady_clock::now();

                improvement_algorithm_time += (double)std::chrono::duration_cast
                    <std::chrono::nanoseconds>(improvement_algorithm_end - improvement_algorithm_begin)
                    .count() / (double)1000000000;

                improvement_algorithm_calls += 1;

                repair_begin =std::chrono::steady_clock::now();
                new_generation[i].repair();
                repair_end = std::chrono::steady_clock::now();

                repair_time += (double)std::chrono::duration_cast
                    <std::chrono::nanoseconds>(repair_end - repair_begin)
                    .count() / (double)1000000000;
                repair_calls +=1;
            }

            //cout << "       iteration " << g << ": elite population processed\n";

            //dobbiamo inserire il rimanente 90% della popolazione
            for (; i < s3; ++i)
            {
                p1 = random_parent(mt);
                std::uniform_int_distribution<unsigned> left(0, p1-1);
                std::uniform_int_distribution<unsigned> right(p1 + 1, (gdata.population_size / 2) + 1);
                if(random_choice(mt))
                    p2 = left(mt);
                else
                    p2 = right(mt);

                //cout << "parents chosen: " << &individuals[p1] << " " << &individuals[p2] << endl;
                //Individual child(instance.vehicles, depots, customers, distance_matrix);

                switch (gdata.crossover)
                {
                case 0:

                    crossover_begin = std::chrono::steady_clock::now(); 
                    new_generation[i].one_point_cross_over(individuals[p1], individuals[p2]);
                    crossover_end = std::chrono::steady_clock::now();

                    crossover_time += (double)std::chrono::duration_cast
                        <std::chrono::nanoseconds>(crossover_end - crossover_begin)
                        .count() / (double)1000000000;

                    crossover_calls += 1;
                    
                    repair_begin = std::chrono::steady_clock::now();
                    new_generation[i].repair();
                    repair_end = std::chrono::steady_clock::now();

                    repair_time += (double)std::chrono::duration_cast
                        <std::chrono::nanoseconds>(repair_end - repair_begin)
                        .count() / (double)1000000000;
                    repair_calls +=1;
                    break;
                case 1:
                    crossover_begin = std::chrono::steady_clock::now(); 
                    new_generation[i].two_point_cross_over(individuals[p1], individuals[p2]);
                    crossover_end = std::chrono::steady_clock::now();

                    crossover_time += (double)std::chrono::duration_cast
                        <std::chrono::nanoseconds>(crossover_end - crossover_begin)
                        .count() / (double)1000000000;

                    crossover_calls += 1;

                    repair_begin =std::chrono::steady_clock::now();
                    new_generation[i].repair();
                    repair_end = std::chrono::steady_clock::now();

                    repair_time += (double)std::chrono::duration_cast
                        <std::chrono::nanoseconds>(repair_end - repair_begin)
                        .count() / (double)1000000000;
                    repair_calls +=1;
                    break;
                case 2:
                    crossover_begin = std::chrono::steady_clock::now(); 
                    new_generation[i].best_order_cross_over(individuals[p1], individuals[p2], best_individual);
                    crossover_end = std::chrono::steady_clock::now();

                    crossover_time += (double)std::chrono::duration_cast
                        <std::chrono::nanoseconds>(crossover_end - crossover_begin)
                        .count() / (double)1000000000;

                    crossover_calls += 1;

                    repair_begin =std::chrono::steady_clock::now();
                    new_generation[i].repair();
                    repair_end = std::chrono::steady_clock::now();

                    repair_time += (double)std::chrono::duration_cast
                        <std::chrono::nanoseconds>(repair_end - repair_begin)
                        .count() / (double)1000000000;
                    repair_calls +=1;
                    break;
                case 3:
                    crossover_begin = std::chrono::steady_clock::now();    
                    new_generation[i].position_based_cross_over(individuals[p1], individuals[p2]);
                    crossover_end = std::chrono::steady_clock::now();

                    crossover_time += (double)std::chrono::duration_cast
                        <std::chrono::nanoseconds>(crossover_end - crossover_begin)
                        .count() / (double)1000000000;

                    crossover_calls += 1;

                    repair_begin =std::chrono::steady_clock::now();
                    new_generation[i].repair();
                    repair_end = std::chrono::steady_clock::now();

                    repair_time += (double)std::chrono::duration_cast
                        <std::chrono::nanoseconds>(repair_end - repair_begin)
                        .count() / (double)1000000000;
                    repair_calls +=1;
                    break;
                case 4:
                    crossover_begin = std::chrono::steady_clock::now();    
                    new_generation[i].uniform_cross_over(individuals[p1], individuals[p2]);
                    crossover_end = std::chrono::steady_clock::now();

                    crossover_time += (double)std::chrono::duration_cast
                        <std::chrono::nanoseconds>(crossover_end - crossover_begin)
                        .count() / (double)1000000000;

                    crossover_calls += 1;

                    repair_begin =std::chrono::steady_clock::now();
                    new_generation[i].repair();
                    repair_end = std::chrono::steady_clock::now();

                    repair_time += (double)std::chrono::duration_cast
                        <std::chrono::nanoseconds>(repair_end - repair_begin)
                        .count() / (double)1000000000;
                    repair_calls +=1;
                    break;
                }

                //cout << "       iteration " << g << ": individual " << &(new_generation[i]) <<": crossover executed\n";

                //mutazione genetica solo con un certo rateo
                if (random_mut(mt) == 0)
                {
                    switch (gdata.mutator)
                    {
                    case 0:
                        mutation_begin = std::chrono::steady_clock::now(); 
                        new_generation[i].swap2();
                        mutation_end = std::chrono::steady_clock::now();

                        mutation_time += (double)std::chrono::duration_cast
                            <std::chrono::nanoseconds>(mutation_end - mutation_begin)
                            .count() / (double)1000000000;
                        
                        mutation_calls += 1;

                        repair_begin = std::chrono::steady_clock::now();
                        new_generation[i].repair();
                        repair_end = std::chrono::steady_clock::now();

                        repair_time += (double)std::chrono::duration_cast
                            <std::chrono::nanoseconds>(repair_end - repair_begin)
                            .count() / (double)1000000000;
                        repair_calls +=1;
                        break;
                    case 1:
                        mutation_begin = std::chrono::steady_clock::now(); 
                        new_generation[i].swap3();
                        mutation_end = std::chrono::steady_clock::now();

                        mutation_time += (double)std::chrono::duration_cast
                            <std::chrono::nanoseconds>(mutation_end - mutation_begin)
                            .count() / (double)1000000000;
                        
                        mutation_calls += 1;

                        repair_begin =std::chrono::steady_clock::now();
                        new_generation[i].repair();
                        repair_end = std::chrono::steady_clock::now();

                        repair_time += (double)std::chrono::duration_cast
                            <std::chrono::nanoseconds>(repair_end - repair_begin)
                            .count() / (double)1000000000;
                        repair_calls +=1;
                        break;
                    case 2:
                        mutation_begin = std::chrono::steady_clock::now(); 
                        new_generation[i].scramble();
                        mutation_end = std::chrono::steady_clock::now();

                        mutation_time += (double)std::chrono::duration_cast
                            <std::chrono::nanoseconds>(mutation_end - mutation_begin)
                            .count() / (double)1000000000;
                        
                        mutation_calls += 1;

                        repair_begin =std::chrono::steady_clock::now();
                        new_generation[i].repair();
                        repair_end = std::chrono::steady_clock::now();

                        repair_time += (double)std::chrono::duration_cast
                            <std::chrono::nanoseconds>(repair_end - repair_begin)
                            .count() / (double)1000000000;
                        repair_calls +=1;
                        break;
                    case 3:
                        mutation_begin = std::chrono::steady_clock::now(); 
                        new_generation[i].inversion();
                        mutation_end = std::chrono::steady_clock::now();

                        mutation_time += (double)std::chrono::duration_cast
                            <std::chrono::nanoseconds>(mutation_end - mutation_begin)
                            .count() / (double)1000000000;
                        
                        mutation_calls += 1;

                        repair_begin =std::chrono::steady_clock::now();
                        new_generation[i].repair();
                        repair_end = std::chrono::steady_clock::now();

                        repair_time += (double)std::chrono::duration_cast
                            <std::chrono::nanoseconds>(repair_end - repair_begin)
                            .count() / (double)1000000000;
                        repair_calls +=1;
                        break;
                    case 4:
                        mutation_begin = std::chrono::steady_clock::now(); 
                        new_generation[i].insertion();
                        mutation_end = std::chrono::steady_clock::now();

                        mutation_time += (double)std::chrono::duration_cast
                            <std::chrono::nanoseconds>(mutation_end - mutation_begin)
                            .count() / (double)1000000000;
                        
                        mutation_calls += 1;
                        
                        repair_begin =std::chrono::steady_clock::now();
                        new_generation[i].repair();
                        repair_end = std::chrono::steady_clock::now();

                        repair_time += (double)std::chrono::duration_cast
                            <std::chrono::nanoseconds>(repair_end - repair_begin)
                            .count() / (double)1000000000;
                        repair_calls +=1;
                        break;
                    }
                }

                //cout << "       iteration " << g << ": individual " << &(new_generation[i]) <<": mutation executed\n";

                
                improvement_algorithm_begin = std::chrono::steady_clock::now();

                new_generation[i].improvement_algorithm();

                improvement_algorithm_end = std::chrono::steady_clock::now();

                improvement_algorithm_time += (double)std::chrono::duration_cast
                    <std::chrono::nanoseconds>(improvement_algorithm_end - improvement_algorithm_begin)
                    .count() / (double)1000000000;

                improvement_algorithm_calls += 1;

                repair_begin =std::chrono::steady_clock::now();
                new_generation[i].repair();
                repair_end = std::chrono::steady_clock::now();

                repair_time += (double)std::chrono::duration_cast
                    <std::chrono::nanoseconds>(repair_end - repair_begin)
                    .count() / (double)1000000000;
                repair_calls +=1;

                new_generation[i].calculate_cost();

                if (new_generation[i].get_cost() < best_cost)
                {
                    best_cost = new_generation[i].get_cost();
                    best_individual = new_generation[i];
                    best_individual2 = new_generation[i];
                    if ((unsigned)best_individual.get_cost() == instance.known_solution)
                    {
                        cout <<"    AUDITING:\n";
                        cout <<"    (function): (total time)\n";
                        cout <<"    random_initalize: " << random_initialize_time << endl;
                        cout <<"    repair: " << repair_time << endl;
                        cout <<"    improvement algorithm: " << improvement_algorithm_time << endl;
                        cout <<"    crossover: "  << crossover_time << endl;
                        cout <<"    mutation: " << mutation_time << endl;
                        cout << "known solution ";
                        //best_individual.print_tour();
                        return best_individual.get_cost();
                    }
                    //std::cout << "Child improved: " << best_cost << "\n";
                }
            }

            //cout << "       iteration " << g << ": entire population processed\n";

            std::vector<Individual> *const temp = individuals_ptr;
            individuals_ptr = new_generation_ptr;
            new_generation_ptr = temp;

            std::sort(new_generation.begin(), new_generation.end());
            
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
                cout <<"    AUDITING:\n";
                cout <<"    (function): (total time)\n";
                cout <<"    random_initalize: " << random_initialize_time << endl;
                cout <<"    repair: " << repair_time << endl;
                cout <<"    improvement algorithm: " << improvement_algorithm_time << endl;
                cout <<"    crossover: "  << crossover_time << endl;
                cout <<"    mutation: " << mutation_time << endl;
                cout << "known solution ";
                //best_individual.print_tour();
                return cost;
            }
        } 
        
    }

    cout <<"    AUDITING:\n";
    cout <<"    (function): (total time)\n";
    cout <<"    random_initalize: " << random_initialize_time << endl;
    cout <<"    repair: " << repair_time << endl;
    cout <<"    improvement algorithm: " << improvement_algorithm_time << endl;
    cout <<"    crossover: "  << crossover_time << endl;
    cout <<"    mutation: " << mutation_time << endl;
    best_individual2.print_tour();
    return cost;
}

#endif