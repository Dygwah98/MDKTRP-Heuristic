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
using distTable = unordered_map<unsigned, double>;

class Individual;

class IndividualSet {

    friend Individual;
    public:
        IndividualSet(unsigned pop_size, Individual seed, distTable& ct, distTable& dt, const double *const ac):
            individuals(pop_size, seed),
            activation_costs( ac ),
            customer_table( ct ),
            depot_table( dt ),
            random_cell( 0, seed.N-1 ),
            random_depot( 0, seed.depots-1 ),
            random_bit( 0, 1 ) 
        {}
    
    private:
        vector<Individual> individuals;
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
};

class Individual {

    friend IndividualSet;
    public:
        Individual(const unsigned v, const unsigned d, const unsigned c, IndividualSet* parent = nullptr):
            vehicles( v ), 
            depots( d ), 
            customers( c ), 
            N( d+c ), 
            vehicles_n( v-1 ), 
            cost( 0 ), 
            tours( new unsigned[c] ),
	        vehicles_positions( new unsigned[vehicles] ),
            depot_vehicles( new unsigned[vehicles] ),
            parent(parent)
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
            parent(o.parent)
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
		
        void position_based_cross_over(const Individual&p1, const Individual&p2) {

            std::uniform_int_distribution<unsigned> random_n(0, customers-1);
            const Individual*const parents[2] = {&p1, &p2};
            
            const unsigned n_positions = random_n(mt);
            unsigned p = parent->random_bit(mt);
            unsigned index = 0;
            unsigned ind = 0;
            unsigned pt = 0;
            
            bool *customers_not_inserted = new bool[customers];
            for(unsigned i = 0; i < customers; ++i) {
                
                customers_not_inserted[i] = true;
            }

            unsigned *const positions = new unsigned[customers];
            for(unsigned i = 0; i < customers; ++i) {
                
                positions[i] = i;
            }

            std::shuffle(positions, positions + customers, mt);

            for(unsigned i = 0; i < n_positions; ++i) {
                
                index = positions[i];
                
                tours[ index ] = parents[p]->tours[ index ];
                customers_not_inserted[ index ] = false;
            }

            p = !p;
            for(unsigned i = n_positions; i < customers; ++i) {

                index = positions[i];
                pt = parents[ p ]->tours[ index ];
                
                if(customers_not_inserted[ index ]) {
                
                    tours[ index ] = pt;
                    customers_not_inserted[ index ] = false;
                
                } else {
                
                    while(ind < customers) {
                
                        pt = parents[ p ]->tours[ ind ];

                        if(customers_not_inserted[ ind ]) {
                
                            tours[ ind ] = pt;
                            customers_not_inserted[ ind ] = false;
                            break;
                        }
                
                        ++ind;
                    }
                }
            }

            delete[] positions;
            delete[] customers_not_inserted;
        }
		
        void uniform_cross_over(const Individual &p1, const Individual &p2) {

            auto& random_bit = parent->random_bit;
            unsigned index[2] = {0, 0};
            const Individual*const parents[2] = {&p1, &p2};
            
            unsigned pt = 0;
            unsigned pv = 0;

            bool *customers_not_inserted = new bool[customers];
            for(unsigned i = 0; i < customers; ++i) {
            
                customers_not_inserted[i] = true;
            }

            for(unsigned i = 0; i < customers; ++i) {

                const unsigned rand = random_bit(mt);
                pt = parents[rand]->tours[i];
            
                if(customers_not_inserted[i]) {
            
                    tours[i] = pt;
                    customers_not_inserted[i] = false;
            
                } else {
            
                    while(index[rand] < customers) {
            
                        pt = parents[rand]->tours[ index[rand] ];
            
                        if(customers_not_inserted[ index[rand] ]) {
            
                            tours[ index[rand] ] = pt;
                            customers_not_inserted[ index[rand] ] = false;
                            break;
                        }
                        ++index[rand];
                    }
                }
            }

            unsigned rand2 = random_bit(mt);
            
            vehicles_positions[0] = parents[rand2]->vehicles_positions[0];
            depot_vehicles[0] = parents[rand2]->depot_vehicles[0];

            for(unsigned i = 1; i < vehicles; ++i) {

                const unsigned rand3 = random_bit(mt);
                pv = parents[rand3]->vehicles_positions[i];

                while(pv <= vehicles_positions[i-1]) {
                    ++pv;
                }

                vehicles_positions[i] = pv;
                depot_vehicles[i] = parents[rand3]->depot_vehicles[i];
            }
        
            delete[] customers_not_inserted;
        }

        //da valutare: solo mosse ammissibili?
        //bool sanity_check()

        void calculate_cost() {

            const unsigned *const tours = this->tours;
            const unsigned *const pv = this->vehicles_positions;
            const double *const ac = this->parent->activation_costs;
            const unsigned *const dv = this->depot_vehicles;
            distTable& dt = this->parent->customer_table;

            double sum = 0;
            for(unsigned v = 0; v < vehicles_n; ++v) {
                
                unsigned first = pv[v];
                const unsigned last = pv[v+1]-1;
                unsigned len = last - first;
                
                sum += ac[ dv[v] ];

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

        inline double get_cost() const {

            return cost;
        }

        inline void set_parent(IndividualSet& set) {
            parent = &set;
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
        //container di individuals
        IndividualSet* parent;

        inline unsigned getIndex(unsigned x, unsigned y) {
            return x*customers + y;
        }
};

#endif