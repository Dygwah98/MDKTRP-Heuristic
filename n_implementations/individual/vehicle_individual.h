#ifndef V_IND_H
#define V_IND_H

#include"tour_individual.h"

class VehicleIndividual {

    private:
        //cardinalità veicoli, customers, depots
        const unsigned vehicles;
        const unsigned customers;
        const unsigned depots;
        //customers * depot
        const unsigned N;
        //costo
        double cost;
        //giant tour
        unsigned *const tours;
        //associa ad ogni inizio subtour il suo depot
        map<unsigned, unsigned> tours_start;
    #ifndef BASE    
        //costi di attivazione per ogni depot
        double *const activation_costs;
    #endif
        //matrice delle coordinate (per calcolare le distanze)
        double**const coordinate_matrix;
        //tabella di hash contente le distanze node <-> customer
        distTable& distance_table;
        //generatori di numeri casuali
        std::uniform_int_distribution<unsigned> random_cell;
		std::uniform_int_distribution<unsigned> random_depot;
        std::uniform_int_distribution<unsigned> random_bit;

    public:
        VehicleIndividual(const unsigned vehicles, const unsigned customers, const unsigned depots, 
                        unsigned*const tours, double*const activation_costs, double**const coordinate_matrix, distTable& distances):
            vehicles(vehicles),
            customers(customers),
            depots(depots),
            N(customers*depots),
            cost(0),
            tours(tours),
            tours_start(),
            activation_costs(activation_costs),
            coordinate_matrix(coordinate_matrix),
            distance_table(distances),
            random_cell(0, customers-1),
            random_depot(0, depots-1),
            random_bit(0, 1)
        {

        }

        VehicleIndividual(const VehicleIndividual& v);
        VehicleIndividual& operator=(const VehicleIndividual& v);
        VehicleIndividual(const TourIndividual& t);

        inline unsigned getCustomerIndex(unsigned x, unsigned y) {
            return N + x*customers + y;
        }

        inline unsigned getDepotIndex(unsigned x, unsigned y) {
            return x*customers + y;
        }

};

#endif