#ifndef T_IND_H
#define T_IND_H

#include"individual.h"

class TourIndividual {

    private:
        const unsigned customers;
        const unsigned depots;
        //customers * depot
        const unsigned N;
        //costo
        double cost;
        //giant tour
        unsigned* tours;
        //matrice delle coordinate (per calcolare le distanze)
        double**const coordinate_matrix;
        //tabella di hash contente le distanze node <-> customer
        distTable& distance_table;
        //generatori di numeri casuali
        std::uniform_int_distribution<unsigned> random_cell;
		std::uniform_int_distribution<unsigned> random_depot;
        std::uniform_int_distribution<unsigned> random_bit;
    
    public:
        TourIndividual(const unsigned customers, const unsigned depots, double**const coordinate_matrix, distTable& distance_table):
            customers(customers),
            depots(depots),
            N(customers * depots),
            cost(0),
            coordinate_matrix(coordinate_matrix),
            distance_table(distance_table),
            random_cell(0, customers-1),
            random_depot(0, depots-1),
            random_bit(0, 1)
        {
            tours = new unsigned[customers];
        }

        TourIndividual(const TourIndividual& t);
        TourIndividual& operator=(const TourIndividual& t);
        TourIndividual(const Individual& o);

        ~TourIndividual() {
            delete[] tours;
        }

        inline unsigned getCustomerIndex(unsigned x, unsigned y) {
            return N + x*customers + y;
        }

        inline unsigned getDepotIndex(unsigned x, unsigned y) {
            return x*customers + y;
        }

};

#endif