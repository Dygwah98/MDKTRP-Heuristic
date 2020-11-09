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

    const GeneticAlgorithmData gdata;
    
    const unsigned max_evaluations = gdata.max_evaluations_GA * instance.factor_valuations;
    const unsigned max_g = (max_evaluations / gdata.tries) / gdata.population_size;
    
    double cost = std::numeric_limits<double>::max();
    double global_best = std::numeric_limits<double>::max();

    Individual best_individual2(ind);
    best_individual2.random_initialize();
    //cout << "       parameters set\n";
    for (unsigned i = 0; i < gdata.tries; ++i)
    {
        double best_cost = std::numeric_limits<double>::max();

        //vettore che contiene la popolazione
        std::vector<Individual>* individuals_ptr = new std::vector<Individual>(gdata.population_size, ind);

        Individual best_individual(best_individual2);
        
        //inizializzazione della popolazione
        for (Individual &i : (*individuals_ptr) )
        {
            i.random_initialize();
            if (i.get_cost() < best_cost)
            {
                best_cost = i.get_cost();
                best_individual = i;
                best_individual2 = i;

                if ((unsigned)best_individual.get_cost() == instance.known_solution)
                {
                    cout << "known solution ";
                    delete individuals_ptr;
                    //best_individual.print_tour();
                    return best_individual.get_cost();                    
                }
            }
        }
        //cout << "       population initialized\n";

        //std::sort(individuals.begin(), individuals.end());

        const unsigned s  = (unsigned)floor( (10.0 * gdata.population_size) / 100.0);
        const unsigned s3 = gdata.population_size;  

        std::uniform_int_distribution<unsigned> random_parent(1, gdata.population_size / 2);
        std::uniform_int_distribution<unsigned> random_mut(0, gdata.mut_rate);
        std::uniform_int_distribution<unsigned> random_choice(0, 1);

        std::vector<Individual>* new_generation_ptr = new std::vector<Individual>(gdata.population_size, ind);

        unsigned g = 0;
        
        unsigned p1;
        unsigned p2; 

        //cout << "       starting evaluations\n";
        while (g < max_g)
        {
            //cout<<"G: "<<g<<"\n";
            const std::vector<Individual>& individuals = *individuals_ptr;
            std::vector<Individual>& new_generation = *new_generation_ptr;
            //elitismo al 10%
            unsigned i = 0;
            for (; i < s; ++i) {
                new_generation[i] = individuals[i];
                //new_generation[i].improvement_algorithm();
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
                    new_generation[i].one_point_cross_over(individuals[p1], individuals[p2]);
                    break;
                case 1:
                    new_generation[i].two_point_cross_over(individuals[p1], individuals[p2]);
                    break;
                case 2:
                    new_generation[i].best_order_cross_over(individuals[p1], individuals[p2], best_individual);
                    break;
                case 3:
                    new_generation[i].position_based_cross_over(individuals[p1], individuals[p2]);
                    break;
                case 4:
                    new_generation[i].uniform_cross_over(individuals[p1], individuals[p2]);
                    break;
                }

                //cout << "       iteration " << g << ": individual " << &(new_generation[i]) <<": crossover executed\n";

                //mutazione genetica solo con un certo rateo
                if (random_mut(mt) == 0)
                {
                    switch (gdata.mutator)
                    {
                    case 0:
                        new_generation[i].swap2();
                        break;
                    case 1:
                        new_generation[i].swap3();
                        break;
                    case 2:
                        new_generation[i].scramble();
                        break;
                    case 3:
                        new_generation[i].inversion();
                        break;
                    case 4:
                        new_generation[i].insertion();
                        break;
                    }
                }

                //cout << "       iteration " << g << ": individual " << &(new_generation[i]) <<": mutation executed\n";

                new_generation[i].calculate_cost();

                if (new_generation[i].get_cost() < best_cost)
                {
                    best_cost = new_generation[i].get_cost();
                    best_individual = new_generation[i];
                    best_individual2 = new_generation[i];
                    if ((unsigned)best_individual.get_cost() == instance.known_solution)
                    {
                        cout << "known solution ";

                        delete individuals_ptr;
                        delete new_generation_ptr;
                        //best_individual.print_tour();
                        return cost;
                    }
                    //std::cout << "Child improved: " << best_cost << "\n";
                }
            }

            //cout << "       iteration " << g << ": entire population processed\n";

            std::swap(individuals_ptr, new_generation_ptr);
            //for(auto& i : individuals)
            //    i.print_tour();

            std::sort(individuals_ptr->begin(), individuals_ptr->end());
            
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
                cout << "known solution ";
                delete individuals_ptr;
                delete new_generation_ptr;
                //best_individual.print_tour();
                return cost;
            }
        } 
    
        delete individuals_ptr;
        delete new_generation_ptr;
    }

    //best_individual2.print_tour();
                
    return cost;
}

#endif