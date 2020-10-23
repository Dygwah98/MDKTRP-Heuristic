#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <iostream>
#include <algorithm>
#include <random>
#include <set>
#include <chrono>
#include <unordered_map>
using namespace std;

std::uint_least64_t seed(11363585407706497491U);     
std::mt19937 mt(seed);

class Individual {

    using distTable = unordered_map<unsigned, double>;

    public:

        Individual(const unsigned v, const unsigned d, const unsigned c, distTable& ct, distTable& dt, const double *const ac):
            vehicles( v ), 
            depots( d ), 
            customers( c ), 
            N( d+c ), 
            vehicles_n( v-1 ), 
            cost( 0 ), 
            tours( new unsigned[c] ),
	        vehicles_positions( new unsigned[vehicles] ),
            depot_vehicles( new unsigned[vehicles] ),
            activation_costs( ac ),
            customer_table( ct ),
            depot_table( dt ),
            random_cell( 0, d+c-1 ),
            random_depot( 0, d-1 ),
            random_bit( 0, 1 )
        {}

        Individual(const Individual& o):
            vehicles( o.vehicles ), 
            depots( o.depots ), 
            customers( o.customers ), 
            N( o.N ), 
            vehicles_n( o.vehicles_n ), 
            cost( 0 ), 
            tours( new unsigned[o.customers] ),
	        vehicles_positions( new unsigned[o.vehicles] ),
            depot_vehicles( new unsigned[o.vehicles] ),
            activation_costs( o.activation_costs ),
            customer_table( o.customer_table ),
            depot_table( o.depot_table ),
            random_cell( 0, o.N - 1 ),
            random_depot( 0, o.depots - 1 ),
            random_bit( 0, 1 )
        {
            for(unsigned i = 0; i < customers; ++i)
                tours[i] = o.tours[i];

            for(unsigned i = 0; i < vehicles; ++i)
                vehicles_positions[i] = o.vehicles_positions[i];

            for(unsigned i = 0; i < vehicles; ++i)
                depot_vehicles[i] = o.depot_vehicles[i];
        }

        void initialize();
        void education();
    	void swap2();
		void swap3();
		void inversion();
		void scrumble();
		void insertion();
		void insertion_repeated();
		double calculate_tour_latency(const unsigned );
		void one_point_cross_over(const Individual&, const Individual&);
		void two_point_cross_over(const Individual&, const Individual&);
		void best_order_cross_over(const Individual&, const Individual&);
		void position_base_cross_over(const Individual&, const Individual&, const Individual&);
		void uniform_cross_over(const Individual&, const Individual&);

        //da valutare: solo mosse ammissibili?
        //bool sanity_check()

        void calculate_cost() {

            const unsigned *const tours = this->tours;
            const unsigned *const pv = this->vehicles_positions;
            const double *const ac = this->activation_costs;
            const unsigned *const dv = this->depot_vehicles;
            distTable& dt = this->customer_table;

            double sum = 0;
            for(unsigned v = 0; v < vehicles_n; ++v) {
                
                sum += ac[ dv[v] ];

                unsigned first = pv[v];
                const unsigned last = pv[v+1]-1;
                unsigned len = last - first;
                for(; first < last; ++first) {
                    
                    const unsigned index = getIndex( tours[first], tours[first+1] );
                    if(dt.find(index) != dt.end()) {
                            sum += dt.at(index) * len;
                    } else {
                        //calcola distanza e inseriscila nella table
                    }

                    --len;
                }
            }

            sum += ac[ dv[vehicles_n] ];

            unsigned first = pv[vehicles_n];
            const unsigned last = pv[vehicles_n+1] - 1;
            unsigned len = last - first;
            for(; first < last; ++first) {
                
                const unsigned index = getIndex( tours[first], tours[first+1] );
                if(dt.find(index) != dt.end()) {
                        sum += dt.at(index) * len;
                } else {
                    //calcola distanza e inseriscila nella table
                }

                --len;
            }

            cost = sum;
        }

        inline double get_cost() const {

            return cost;
        }

        void print_tour_matrix() const {

            cout << "depot : vehicle : tour\n";
            unsigned j;
            for(unsigned i = 0; i < vehicles_n; ++i) {
                
                cout << depot_vehicles[i] << " : ";
                cout << i << " : ";
                
                j = vehicles_positions[i];
                const unsigned max = vehicles_positions[i+1];
                for(; j < max; ++j) {
                    cout << tours[j] << " ";
                }

                cout << endl;
            }

            cout << depot_vehicles[vehicles_n] << " : ";
            cout << vehicles_n << " : ";
                
            j = vehicles_positions[vehicles_n];
            const unsigned max = customers;
            for(; j < max; ++j) {
                cout << tours[j] << " ";
            }

            cout << endl;
        }

        Individual& operator=(const Individual& m) {

            if(this == &m)
                return *this;

            cost = m.cost;

            unsigned *const tours_t = this->tours;
            const unsigned *const tours_m = m.tours;
            for(unsigned i = 0; i < customers; ++i) {

                tours_t[i] = m.tours[i];
            }

            unsigned *const vehicles_t = this->vehicles_positions;
            const unsigned *const vehicles_m = m.vehicles_positions;
            for(unsigned i = 0; i < vehicles; ++i) {

                vehicles_t[i] = vehicles_m[i];
            }

            unsigned *const depots_t = this->depot_vehicles;
            const unsigned *const depots_m = m.depot_vehicles;
            for(unsigned i = 0; i < vehicles; ++i) {

                depots_t[i] = depots_m[i];
            }

            return *this;
        }

        bool operator<(const Individual& o) const {

            return cost < o.cost;
        }

        bool operator==(const Individual& o) const {

            if(this == &o)
                return true;

            const unsigned *const tours_t = this->tours;
            const unsigned *const tours_o = o.tours;
            for(unsigned i = 0; i < customers; ++i) {
                
                if(tours_t[i] != tours_o[i]) {
                    return false;
                }
            }

            const unsigned *const vehicles_t = this->vehicles_positions;
            const unsigned *const vehicles_o = o.vehicles_positions;
            for(unsigned i = 0; i < vehicles; ++i) {
                
                if(vehicles_t[i] != vehicles_o[i]) {
                    return false;
                }
            }

            const unsigned *const depots_t = this->depot_vehicles;
            const unsigned *const depots_o = o.depot_vehicles;
            for(unsigned i = 0; i < vehicles; ++i) {
                
                if(depots_t[i] != depots_o[i]) {
                    return false;
                }
            }

            return true;
        }

        ~Individual() {
            delete[] tours;
            delete[] vehicles_positions;
            delete[] depot_vehicles;
        }
        
    private:
        //cardinalità veicoli, customers, depots
        const unsigned vehicles;
        const unsigned customers;
        const unsigned depots;
        //customers + depot
        const unsigned N;
        //vehicles - 1
        const unsigned vehicles_n;
        //costo
        double cost;
        //giant tour
        unsigned* tours;
        //posizioni di ogni veicolo
        unsigned* vehicles_positions;
        //per ogni veicolo, indica il depot di appartenenza
        unsigned* depot_vehicles;
        //costi di attivazione per ogni depot
        const double *const activation_costs;
        //tabella di hash contente le distanze depot <-> customer
        distTable& depot_table;
        //tabella di hash contente le distanze customer <-> customer
        distTable& customer_table;
        //generatori di numeri casuali
        std::uniform_int_distribution<unsigned> random_cell;
		std::uniform_int_distribution<unsigned> random_depot;
        std::uniform_int_distribution<unsigned> random_bit;
        
        inline unsigned getIndex(unsigned x, unsigned y) {
            return x*customers + y;
        }
};

#endif