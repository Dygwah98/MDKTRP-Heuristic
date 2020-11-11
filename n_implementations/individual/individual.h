#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <iostream>
#include <algorithm>
#include <random>
#include <set>
#include <chrono>
#include <map>
#include <unordered_map>
#include "../utils.h"
using namespace std;

std::uint_least64_t seed(11363585407706497491U);     
std::mt19937 mt(seed);

using distTable = unordered_map<unsigned, double>;

class Individual {

    public:
        Individual(const unsigned v, const unsigned d, const unsigned c, distTable& dt, 
    #ifndef BASE    
        double *const ac, 
    #endif    
        double **const _coordinate_matrix):
            vehicles( v ), 
            depots( d ), 
            customers( c ), 
            N( d*c ), 
            cost( 0 ), 
            tours( new unsigned[c] ),
	        tours_start(),
            coordinate_matrix(_coordinate_matrix),
        #ifndef BASE
            activation_costs( ac ),
        #endif
            distance_table( dt ),
            random_cell( 0, c-1 ),
            random_depot( 0, d-1 ),
            random_bit( 0, 1 )
        {
            //cout << "           default constructor called\n";
            random_initialize();
        }

        Individual(const Individual& o):
            vehicles( o.vehicles ), 
            depots( o.depots ), 
            customers( o.customers ), 
            N( o.N ), 
            cost( o.cost ), 
            tours( new unsigned[o.customers] ),
            tours_start(o.tours_start),
            coordinate_matrix(o.coordinate_matrix),
        #ifndef BASE
            activation_costs( o.activation_costs ),
        #endif
            distance_table( o.distance_table ),
            random_cell( 0, o.customers-1 ),
            random_depot( 0, o.depots-1 ),
            random_bit( 0, 1 )
        {
            //cout << "           copy constructor called\n";

            for(unsigned i = 0; i < customers; ++i)
                tours[i] = o.tours[i];
        }

        Individual& operator=(const Individual& m) {

            if(this == &m)
                return *this;

            unsigned *const tours_t = this->tours;
            const unsigned *const tours_m = m.tours;
            
            for(unsigned i = 0; i < customers; ++i) {

                tours_t[i] = m.tours[i];
            }

            tours_start = m.tours_start;

            random_cell = m.random_cell;
            random_depot = m.random_depot;
            random_bit = m.random_bit;

            if(m.cost != 0) {
                cost = m.cost;
            }
            else
                calculate_cost();

            return *this;
        }

        ~Individual() {
            delete[] tours;
        }
    	
        void swap2()   {
/*
            auto& random_cell = this->random_cell;

            unsigned fst = random_cell(mt);
            
            std::uniform_int_distribution<unsigned> left(-1, fst-1);
            std::uniform_int_distribution<unsigned> right(fst+1, customers);
            
            unsigned snd;
            if(random_bit(mt)) {
                snd = left(mt);
                if(snd == -1) ++snd;
            } else {
                snd = right(mt);
                if(snd == customers) --snd;
            } 

            const unsigned node = tours[fst];
            tours[fst] = tours[snd];
            tours[snd] = node;
*/

            if(random_bit(mt)) {
                //interroute swap

                const auto& tours_start = this->tours_start;

                std::uniform_int_distribution<unsigned> route(0, tours_start.size()-1);

                auto it = next(tours_start.begin(), route(mt));
                const unsigned route_c = it->first; 
                ++it;
                const unsigned route_e = it == tours_start.end() ? customers : it->first;

                std::uniform_int_distribution<unsigned> choice(route_c, route_e-1);
                const unsigned fst = choice(mt);
                const unsigned snd = choice(mt);

                const unsigned node = tours[fst];
                tours[fst] = tours[snd];
                tours[snd] = node;

            } else {

                auto& random_cell = this->random_cell;

                const unsigned fst = random_cell(mt);
                const unsigned snd = random_cell(mt);

                const unsigned node = tours[fst];
                tours[fst] = tours[snd];
                tours[snd] = node;
            }

            //cout << "           individual " << this << ": swap2 executed\n";
        }
		
        void swap3()   {

            auto& random_cell = this->random_cell;

            unsigned fst = random_cell(mt);
            unsigned snd = random_cell(mt);

            while(snd == fst) {
                snd = random_cell(mt);
            }

            unsigned trd = random_cell(mt);

            while(trd == fst || trd == snd) {
                trd = random_cell(mt);
            }

            const unsigned node = tours[fst];
            tours[fst] = tours[snd];
            tours[snd] = tours[trd];
            tours[trd] = node;
        }
		
        void inversion()   {

            auto& random_cell = this->random_cell;

            unsigned start = random_cell(mt);
            unsigned end = random_cell(mt);

            while(end == start) {
                end = random_cell(mt);
            }

            if(end < start)
                std::swap(start, end);

            unsigned *const tours = this->tours;
            std::reverse(tours + start, tours + end);

        }
		
        void scramble()   {

            auto& random_cell = this->random_cell;

            unsigned start = random_cell(mt);
            unsigned end = random_cell(mt);

            while(end == start) {
                end = random_cell(mt);
            }

            if(end < start)
                std::swap(start, end);

            unsigned *const tours = this->tours;
            std::shuffle(tours + start, tours + end, mt);

        }
		
        void insertion()   {

            auto& dt = this->tours_start;
            const unsigned position = this->random_cell(mt);
            const unsigned node = tours[position];

            unsigned scroll1 = position;
            unsigned pos = node;
            while(dt.find(pos) == dt.end()) {
                pos = tours[--scroll1];
            }

            unsigned scroll2 = position;
            unsigned end = node;
            while(scroll2 < customers && dt.find(end) == dt.end()) {
                end = tours[++scroll2];
            }
            
            //viene cercato un depot migliore
            if( random_bit(mt) ) {

                unsigned best_depot = dt.at(pos);
                double min_latency = calculate_tour_cost(scroll1, scroll2);
                bool improved = false;
                for(unsigned i = 0; i < depots; ++i) {
                    
                    dt[pos] = i;
                    
                    double latency = calculate_tour_cost(scroll1, scroll2);
                    if(latency < min_latency) {

                        best_depot = i;
                        min_latency = latency;
                        improved = true;
                    }
                }

                if(!improved) {
                    dt[pos] = best_depot;
                }
            
            } 
            //viene cercata una posizione migliore
            else {

                unsigned best_position = scroll1;
                double min_latency = calculate_tour_cost(scroll1, scroll2);
                for(unsigned swap = scroll1; swap < scroll2; ++swap) {

                    tours[position] = tours[swap];
                    tours[swap] = node;

                    double latency = calculate_tour_cost(scroll1, scroll2);
                    if(latency < min_latency) {
                        min_latency = latency;
                        best_position = swap;
                    } else {
                        tours[swap] = tours[position];
                        tours[position] = node;
                    }
                }
            }

        }
		
        void one_point_cross_over(const Individual& p1, const Individual& p2)   {
            
            unsigned *const tours = this->tours;
		    const unsigned *const tours_p1 = p1.tours;
		    const unsigned *const tours_p2 = p2.tours;
            auto& random_cell = this->random_cell;      

            const unsigned cutting_point = random_cell(mt);

            unsigned index_child = 0;
            for(; index_child < cutting_point; ++index_child) {
                tours[index_child] = tours_p1[index_child];
            }

            for(; index_child < customers; ++index_child) {
                tours[index_child] = tours_p2[index_child];
            }

            auto& vmap = this->tours_start;
            auto& vmap_p1 = p1.tours_start;
            auto& vmap_p2 = p2.tours_start;
            const unsigned bound = min(vmap_p1.size(), vmap_p2.size());
            std::uniform_int_distribution<unsigned> random_vehicle(0, bound);
            const unsigned vcutting_point = random_vehicle(mt);

            vmap.clear();
            
            auto p1b = vmap_p1.begin();
            auto vcp_1 = next(p1b, vcutting_point); 
            vmap.insert(p1b, vcp_1);
           
            auto p2b = next(vmap_p2.begin(), vcutting_point);
            auto vcp_2 = vmap_p2.end();
            vmap.insert(p2b, vcp_2);
        
        }
		
        void two_point_cross_over(const Individual& p1, const Individual& p2)   {
/*       
            unsigned *const tours = this->tours;
            const unsigned *const tours_p1 = p1.tours;
            const unsigned *const tours_p2 = p2.tours;
            auto& random_cell = this->random_cell;

            unsigned cutting_point1 = random_cell(mt);
            unsigned cutting_point2 = random_cell(mt);

            if(cutting_point1 == customers - 1) {
                while(cutting_point2 == cutting_point1) {
                    cutting_point2 = random_cell(mt);
                }
            }

            unsigned first_cutting_point = std::min(cutting_point1, cutting_point2);
            unsigned last_cutting_point = std::max(cutting_point1, cutting_point2);
        
            unsigned i = 0;
            for(; i < first_cutting_point; ++i) {
                
                tours[i] = tours_p2[i];
            }

            
            for(; i < last_cutting_point; ++i) {
                
                tours[i] = tours_p1[i];
                
            }

            for(; i < customers; ++i) {
                
                tours[i] = tours_p2[i];
            }

            auto& vp = this->tours_start;
            auto& vp_p1 = p1.tours_start;
            auto& vp_p2 = p2.tours_start;
            const unsigned bound = min(vp_p1.size(), vp_p2.size());
            std::uniform_int_distribution<unsigned> random_v(0, bound);

            unsigned vcutting_point1 = random_v(mt);

            if(vcutting_point1 == bound) {
                while(vcutting_point1 == bound) {
                   vcutting_point1 = random_v(mt);
                }
            }

            auto& first = (vp_p1.size() == bound) ? vp_p2 : vp_p1;
            auto& second = (vp_p1.size() == bound) ? vp_p1 : vp_p2;

            vp.clear();

            //copio il primo segmento di veicoli
            auto fb = first.begin();
            auto cp1 = next(fb, vcutting_point1);
            vp.insert(fb, cp1);
            //copio il secondo segmento di veicoli
            auto sb = next(second.begin(), vcutting_point1);
            auto cp2 = second.end();
            vp.insert(sb, cp2);
            //copio il terzo segmento di veicoli
            auto fb2 = next(fb, second.size());
            auto fe = first.end();
            vp.insert(fb2, fe);
*/

            auto& ts = this->tours_start;
            map<unsigned, unsigned>().swap(ts);

            ts.insert(p1.tours_start.begin(), p1.tours_start.end());
            ts.insert(p2.tours_start.begin(), p2.tours_start.end());

            const unsigned size = ts.size();
            if(size >= 3) {

                const unsigned bound = size > 3 ? size - 3 : 1;
                std::uniform_int_distribution<unsigned> cpoint1(0, bound);
                std::uniform_int_distribution<unsigned> cpoint2(bound + 1, size-1);

                const unsigned first_cutting_point = next(ts.begin(), cpoint1(mt))->first;
                const unsigned second_cutting_point = next(ts.begin(), cpoint2(mt))->first;

                unsigned *const tours = this->tours;
                const unsigned *const toursp1 = p1.tours;
                const unsigned *const toursp2 = p2.tours;
                unsigned i = 0;
                for(; i < first_cutting_point; ++i) {
                    tours[i] = toursp1[i];
                }
                for(; i < second_cutting_point; ++i) {
                    tours[i] = toursp2[i];
                }
                for(; i < customers; ++i) {
                    tours[i] = toursp1[i];
                }

            } 
            
        }
		
        void best_order_cross_over(const Individual&p1, const Individual&p2, const Individual& best)   {

            std::uniform_int_distribution<unsigned> len_random_cutting_point(1, customers / 3);
		    std::uniform_int_distribution<unsigned> random_value_sequence(0, 2);
		    const Individual*const parents[3] = {&p1, &p2, &best};
            unsigned len_1 = len_random_cutting_point(mt);
		    unsigned cutting_point_1 = 0 + len_1;
            const unsigned value_sequence = random_value_sequence(mt);
            unsigned last_cutting_point = cutting_point_1;
		    unsigned len_before_end = customers - last_cutting_point;

            bool *customers_in_sequence = new bool[customers];
            for (unsigned i = depots; i < customers; ++i)
            {
                customers_in_sequence[i] = false;
            }
            
            if (value_sequence == 0)
            {
                //copia interamente la sequenza dal genitore principale
                for (unsigned i = 0; i < cutting_point_1; ++i)
                {
                    const unsigned node = p1.tours[i];
                    tours[i] = node;
                }
                
            }
            else
            {
                //copia la sequenza dal genitore principale ma con l'ordine dettato da un altro genitore
                const unsigned *other_tours = p2.tours;
                if (value_sequence == 2)
                {
                    other_tours = best.tours;
                }

                for (unsigned i = 0; i < cutting_point_1; ++i)
                {
                    customers_in_sequence[p1.tours[i]] = true;
                }

                for (unsigned index_child = 0; index_child < cutting_point_1; ++index_child)
                {
                    const unsigned node = other_tours[index_child];
                    
                    if (customers_in_sequence[node])
                    {
                        tours[index_child] = node;
                        customers_in_sequence[node] = false;
                    }
                    
                }
            }

            const Individual& p = *(parents[value_sequence]);
            auto p1b = p.tours_start.begin();
            auto p1e = p1b;
            while(p1e != p.tours_start.end() && p1e->first < cutting_point_1) {
                ++p1e;
            }
            tours_start.clear();
            tours_start.insert(p1b, p1e);

            while (len_before_end > customers / 3)
            {
                unsigned next_cutting_point = last_cutting_point + len_random_cutting_point(mt);

                //scegli una modalità
                const unsigned value_sequence = random_value_sequence(mt);
                if (value_sequence == 0)
                {
                    //copia interamente la sequenza dal genitore principale
                    for (unsigned i = last_cutting_point; i < next_cutting_point; ++i)
                    {
                        const unsigned node = p1.tours[i];
                        tours[i] = node;
                    }
                }
                else
                {
                    //copia la sequenza dal genitore principale ma con l'ordine dettato da un altro genitore
                    const unsigned *other_tours = p2.tours;
                    if (value_sequence == 2)
                    {
                        other_tours = best.tours;
                    }

                    //conta il numero di depot presenti nella sequenza
                    for (unsigned i = last_cutting_point; i < next_cutting_point; ++i)
                    {
                        const unsigned node = p1.tours[i];
                        customers_in_sequence[node] = true;
                        
                    }

                    
                    for (unsigned index_child = last_cutting_point; index_child < next_cutting_point; ++index_child)
                    {
                        const unsigned node = other_tours[index_child];
                        //voglio inserire un customer
                        if (customers_in_sequence[node])
                        {
                            tours[index_child] = node;
                            customers_in_sequence[node] = false;
                        }
                    
                    }
                }

                auto p2e = p1e;
                while(p2e != p.tours_start.end() && p2e->first < next_cutting_point) {
                    ++p2e;
                }
                tours_start.insert(p1e, p2e);

                p1e = p2e;

                last_cutting_point = next_cutting_point;
                len_before_end = customers - last_cutting_point;
            }

            const unsigned end_sequence_value = random_value_sequence(mt);
            if (end_sequence_value == 0)
            {
                //copia interamente la sequenza dal genitore principale
                for (unsigned i = last_cutting_point; i < customers; ++i)
                {
                    const unsigned node = p1.tours[i];
                    tours[i] = node;
                }
            }
            else
            {
                //copia la sequenza dal genitore principale ma con l'ordine dettato da un altro genitore
                const unsigned *other_tours = p2.tours;
                if (value_sequence == 2)
                {
                    other_tours = best.tours;
                }

                //conta il numero di depot presenti nella sequenza
                unsigned vehicles_insertable = 0;
                for (unsigned i = last_cutting_point; i < customers; ++i)
                {
                    const unsigned node = p1.tours[i];
                    customers_in_sequence[node] = true;
                }

                
                for (unsigned index_child = last_cutting_point; index_child < customers; ++index_child)
                {
                    const unsigned node = other_tours[index_child];
                    
                    //voglio inserire un customer
                    if (customers_in_sequence[node])
                    {
                        tours[index_child] = node;
                
                        customers_in_sequence[node] = false;
                    }
                
                }
            }

            tours_start.insert(p1e, p.tours_start.end());
            
		    delete[] customers_in_sequence;

        }
		
        void position_based_cross_over(const Individual&p1, const Individual&p2)   {

            unsigned *const tours = this->tours;
            const Individual *const parents[2] = {&p1, &p2};
            auto& random_n = this->random_cell;
            auto& random_bit = this->random_bit;

            const unsigned n_positions = random_n(mt);
            unsigned p = random_bit(mt);

            set<unsigned> swap_pos;
            for(unsigned i = 0; i < n_positions; ++i) {
                swap_pos.insert(random_n(mt));
            }

            for(unsigned i = 0; i < customers; ++i) {
                tours[i] = parents[p]->tours[i];
            }

            p = (p + 1) % 2;
            for(auto& swap_index : swap_pos) {
                
                bool swap = false;
                const unsigned swap_2 = parents[p]->tours[ swap_index ];
                
                unsigned swap_1 = tours[ swap_index ];
                tours[ swap_index ] = swap_2;
                
                for(unsigned j = 0; !swap && j < swap_index; ++j) {
                    if(tours[j] == swap_2 ) {
                        tours[j] = swap_1;
                        swap = true;
                    }
                }
                
                for(unsigned j = swap_index + 1; !swap && j < customers; ++j) {
                    if(tours[j] == swap_2 ) {
                        tours[j] = swap_1;
                        swap = true;
                    }
                }
            }

            unsigned vp = random_bit(mt);
            unsigned nvp = (vp + 1) % 2;
            const unsigned vdim = min(parents[nvp]->tours_start.size(), parents[p]->tours_start.size());
            auto random_v = std::uniform_int_distribution<unsigned>(0, vdim);
            
            const unsigned vn_positions = random_v(mt);
            
            set<unsigned> vswap_pos;
            for(unsigned i = 0; i < vn_positions; ++i) {
                vswap_pos.insert(random_v(mt));
            }

            tours_start.clear();
            tours_start.insert(parents[vp]->tours_start.begin(), parents[vp]->tours_start.end());

            for(auto& i : vswap_pos) {
                tours_start[i] = parents[vp]->tours_start.at(i);
            }

        }
		
        void uniform_cross_over(const Individual &p1, const Individual &p2)   {

            auto& random_bit = this->random_bit;
            unsigned index[2] = {0, 0};
            const Individual*const parents[2] = {&p1, &p2};

            for(unsigned i = 0; i < customers; ++i) {
                
                tours[i] = parents[random_bit(mt)]->tours[i];
                
            }

            const unsigned vdim = min(p1.tours_start.size(), p2.tours_start.size());

            tours_start.clear();
            for(unsigned i = 0; i < vdim; ++i) {

                const unsigned rand3 = random_bit(mt);
                const auto it = parents[rand3]->tours_start.begin();
                tours_start.insert( *(next(it, i)) );
            }
        
        }

        //per ora ignoro il problema dello splitting
        void random_initialize()   {
            
            auto& random_depot = this->random_depot;

            //inserisco i customers in maniera sequenziale
            for(unsigned i = 0; i < customers; ++i) {
                tours[i] = i;
            }

            //randomizzo le posizioni di ogni customer
            std::shuffle(tours, tours + customers, mt);

            //stabilisco il numero di veicoli da inserire
            std::uniform_int_distribution<unsigned> r_vehicles(1, vehicles);
            const unsigned vn = r_vehicles(mt);

            //inserisco in maniera ordinata i veicoli, selezionando posizioni randomiche per i depot
            unsigned start = 0;
            for(unsigned i = 0; i < vn - 1; ++i) {
                std::uniform_int_distribution<unsigned> r_nextstop(start, customers-vn-2);

                const unsigned end = r_nextstop(mt);
                tours_start[start] = random_depot(mt);

                start = end;
            }
            tours_start[start] = random_depot(mt);

            //cout << "           individual " << this << " random tour:\n";
            //print_tour();

            improvement_algorithm();
            calculate_cost();

            //cout << "           individual " << this << " randomly initialized\n";
        }
        
        void improvement_algorithm()   {

/*  prima parte: neighboorhood search su giant tour */    
            unsigned *const tours = this->tours;
            auto& dt = this->distance_table;
            //il neighboor è definito qui come:
            //un tour ottenibile tramite un 2-swap di nodi nell'originale
            bool improved = false;
            for(unsigned i = 0; !improved && i+1 < customers; ++i) {
                const unsigned node = tours[i];
                const unsigned next = tours[i+1];
                for(unsigned j = i+1; !improved && j+1 < customers; ++j) {
                    
                    const unsigned jnode = tours[j];
                    const unsigned jnext = tours[j+1];

                    double variation;

                    const unsigned d_node = node + depots;
                    const unsigned d_next = next + depots;
                    const unsigned d_jnode = jnode + depots;
                    const unsigned d_jnext = jnext + depots;
                    
                    const unsigned index_i = getCustomerIndex(node, next);
                    if(dt.find(index_i) == dt.end()) {
                        dt[index_i] = euclidean_distance(coordinate_matrix[d_node][0], coordinate_matrix[d_node][1],
                                                        coordinate_matrix[d_next][0], coordinate_matrix[d_next][1]);
                    }
                    variation += dt.at(index_i);

                    const unsigned index_j = getCustomerIndex(jnode, jnext);
                    if(dt.find(index_j) == dt.end()) {
                        dt[index_j] = euclidean_distance(coordinate_matrix[d_jnode][0], coordinate_matrix[d_jnode][1],
                                                        coordinate_matrix[d_jnext][0], coordinate_matrix[d_jnext][1]);
                    }
                    variation += dt.at(index_j);

                    const unsigned n_index_i = getCustomerIndex(jnode, next);
                    if(dt.find(n_index_i) == dt.end()) {
                        dt[n_index_i] = euclidean_distance(coordinate_matrix[d_jnode][0], coordinate_matrix[d_jnode][1],
                                                        coordinate_matrix[d_next][0], coordinate_matrix[d_next][1]);
                    }
                    variation -= dt.at(n_index_i);

                    const unsigned n_index_j = getCustomerIndex(node, jnext);
                    if(dt.find(n_index_j) == dt.end()) {
                        dt[n_index_j] = euclidean_distance(coordinate_matrix[d_node][0], coordinate_matrix[d_node][1],
                                                        coordinate_matrix[d_jnext][0], coordinate_matrix[d_jnext][1]);
                    }
                    variation -= dt.at(n_index_j);

                    //se viene trovato un neighboor migliore la search termina
                    if(variation < -0.9) {
                        //improved = true;

                        tours[i] = tours[j];
                        tours[j] = node;  
                    }
                }
            }

/* seconda parte: ottimizzazione del costo del depot (sia attivazione che primo arco) */
        #ifndef BASE
            const auto& ac = this->activation_costs;
        #endif
            auto& ts = this->tours_start;

            for(auto& it : ts) {
                
                const unsigned depot = it.second;
                for(unsigned i = 0; i < depots; ++i) {
        #ifndef BASE        
                    //se il depot ha un costo di attivazione minore
                    if(ac[i] < ac[depot]) {
                        //si somma il costo di attivazione col costo dell'arco depot->primo customer del subtour
                        const unsigned first = tours[it.first];

                        const double nval = ac[i] + dt.at( getDepotIndex(i, first) );
                        const double oval = ac[depot] + dt.at(  getDepotIndex(depot, first) );

                        //se il costo complessivo è minore, si effettua lo swap di depot
                        if(nval < oval) {
                            it.second = i;
                            break;
                        }
                    }
        #else
                    //si somma il costo di attivazione col costo dell'arco depot->primo customer del subtour
                    const unsigned first = tours[it.first] + depots;
                    const unsigned dindex1 = getDepotIndex(i, first - depots );
                    if(dt.find(dindex1) == dt.end()) {
                        dt[dindex1] = euclidean_distance(coordinate_matrix[i][0], coordinate_matrix[i][1],
                                                        coordinate_matrix[first][0], coordinate_matrix[ first ][1]);
                    }

                    const unsigned dindex2 = getDepotIndex(depot, first - depots );
                    if(dt.find(dindex2) == dt.end()) {
                        dt[dindex2] = euclidean_distance(coordinate_matrix[depot][0], coordinate_matrix[depot][1],
                                                        coordinate_matrix[ first ][0], coordinate_matrix[ first ][1]);
                    }
                    const double nval = dt.at(dindex1);
                    const double oval = dt.at(dindex2);

                    //se il costo complessivo è minore, si effettua lo swap di depot
                    if(nval < oval) {
                        it.second = i;
                        break;
                    }
                    
        #endif
                }
            }

            //SPLITTING ALGORITHM TIME

        }

        void repair()   {

            //cout << "PREREPAIR\n";
            //print_tour();

            //array di interi per contare il numero di occorrenze di ogni nodo
            unsigned *inserted = new unsigned[customers];
            for(unsigned i = 0; i < customers; ++i) {
                inserted[i] = 0;
            } 

            //conto il numero di occorrenze di ogni nodo
            for(unsigned i = 0; i < customers; ++i) {
                inserted[ tours[i] ] += 1; 
            }

            //salvo il primo nodo non inserito nel tour
            unsigned pos = 0;
            for(; pos < customers; ++pos) {
                if( inserted[pos] == 0 ) {
                    break;
                }
            }

            //riparo il tour
            for(unsigned si = 0; si < customers; ++si) {
                //se tutti i nodi sono stati inseriti, chiusura anticipata
                if(pos >= customers) {
                        break;
                }
                
                //sostituisce un nodo ripetuto con un nodo non inserito
                if( inserted[ tours[si] ] > 1) {

                    //decrementa le occorrenze del nodo
                    inserted[ tours[si] ] -= 1;
                    tours[si] = pos;
                    //incrementa le occorrenze del nodo appena inserito
                    inserted[ pos ] += 1;

                    //trova il successivo nodo non inserito
                    while(inserted[pos] > 0 && pos < customers) {
                        ++pos;
                    }
                }
            }

            delete[] inserted;
            
            //usando una map per rappresentare i veicoli, l'ordine è già mantenuto
            //bisogna limitare il numero di veicoli

            auto& ts = this->tours_start;
            const unsigned size = ts.size();
            if(size > vehicles) {

                unsigned *pos = new unsigned[size];
                unsigned *len = new unsigned[size];

                unsigned k = 0;
                for(auto it = ts.begin(); it != ts.end(); ++it) {
                    
                    auto next_it = it;
                    ++next_it;
                    const unsigned start = it->first;
                    const unsigned end = next_it == ts.end() ? customers : next_it->first;

                    pos[k] = start;
                    len[k] = end - start;

                    ++k;
                }

                int i, j;
                for (i = 1; i < size; i++) {
                        unsigned tmp = len[i];
                        unsigned tmp2 = pos[i];
                        for (j = i; j >= 1 && tmp < len[j-1]; j--) {
                            len[j] = len[j-1];
                            pos[j] = pos[j-1];
                        }
                        len[j] = tmp;
                        pos[j] = tmp2;
                }

                const unsigned excess = size - vehicles;
                for(unsigned l = 0; l < excess; ++l) {
                    ts.erase(pos[l]);
                }

                //se nessun veicolo parte dal primo subtour, si sposta il primo veicolo
                if(ts.find(0) == ts.end()) {
                    const unsigned depot_temp = ts.begin()->second;
                    ts.erase(ts.begin()->first);
                    ts[0] = depot_temp;
                }

                delete[] pos;
                delete[] len;

            }


            //cout << "           individual " << this << " repaired\n";
        }

        double calculate_tour_cost(const unsigned start_pos, const unsigned end_pos)   {

            auto& ts = this->tours_start;
            auto& dc = this->distance_table;
        #ifndef BASE
            const double *const ac = this->activation_costs;
        #endif

            unsigned start = start_pos;
            const unsigned end = end_pos;
            const unsigned depot = ts.at(start_pos);
            unsigned len = end - start;

            double sum = 0;
        #ifndef BASE
            sum += ac[ depot ];
        #endif
            const unsigned tsindex = tours[start] + depots;
            const unsigned di = getDepotIndex(depot, tsindex - depots);
            if(dc.find(di) == dc.end()) {
               dc[di] = euclidean_distance(coordinate_matrix[depot][0], coordinate_matrix[depot][1],
                                            coordinate_matrix[ tsindex ][0], coordinate_matrix[ tsindex ][1]);
            }
            sum += dc.at(di)*len;
            --len;

            //cout << "       depot activation + first arc costs: " << sum << endl;

            for(; start + 1 < end; ++start) {
                
                const unsigned c1 = tours[start] + depots;
                const unsigned c2 = tours[start+1] + depots;
                const unsigned ci = getCustomerIndex(c1 - depots, c2 - depots);
                if(dc.find(ci) == dc.end()) {
                    dc[ci] = euclidean_distance(coordinate_matrix[c1][0], coordinate_matrix[c1][1],
                                                coordinate_matrix[c2][0], coordinate_matrix[c2][1]);
                }
                sum += dc.at(ci)*len;
                --len;

                //cout << "       sum: " << sum << endl;
            }

            if(len < 0) 
                cout << "ERROR LEN\n";

            //cout << "           individual " << this << " subtour: " << start_pos << "to: " << end_pos << "calculated\n";

            return sum;
 
        }

        void calculate_cost()   {

            //print_tour();

            double sum = 0;
            auto& ts = this->tours_start;
            
            auto it = tours_start.begin();
            auto end = tours_start.end();
            auto next_it = next(it, 1);

            while(next_it != end) {
                
                sum += calculate_tour_cost(it->first, next_it->first);
                //cout << "sum : " << sum << endl;

                ++it;
                ++next_it;
            }

            if(next_it == end) {

                sum += calculate_tour_cost(it->first, customers);
                

                //cout << "sum : " << sum << endl;
            }

            double oldcost = cost;
            cost = sum;

            //cout << "           individual " << this << " cost updated from: " << oldcost <<" to: " << cost <<"\n";
        }
 
        void print_tour() const   {

            cout << "\ngiant tour of individual: " << this << "\n";

            for(unsigned i = 0; i < customers; ++i)
                cout << tours[i] << " ";
            cout << "\n";

            cout << "{depot : vehicle : tour}\n";
            
            unsigned i = 0;
            unsigned j;
            const auto end = tours_start.end();
            for(auto it = tours_start.begin(); it != end; ++it) {

                cout << it->second << " : ";
                cout << i << " : ";

                j = it->first;
                auto it2(it);
                const unsigned max = (++it2 == end) ? customers : it2->first;
                for(; j < max; ++j) {
                    cout << tours[j] << " ";
                }

                cout << endl;
                ++i;
            }
            cout << endl;
        }

        inline double get_cost() const {

            return cost;
        }

        bool operator<(const Individual& o) const   {

            //cout << "           individual " << this << " operator< called versus " << &o << "\n";
            //print_tour();
            //cout << cost << "\n";
            //o.print_tour();
            //cout << o.cost << "\n";

            if(this == &o) {
                return true;
            }

            return (int)floor(cost) < (int)floor(o.cost);
        }

        bool operator==(const Individual& o) const   {
            
            //cout << "           individual " << this << " operator== called\n";

            if(this == &o)
                return true;

            return (int)floor(cost) == (int)floor(o.cost);
/*
            const unsigned *const tours_t = this->tours;
            const unsigned *const tours_o = o.tours;
            
            for(unsigned i = 0; i < customers; ++i) {
                
                if(tours_t[i] != tours_o[i]) {
                    return false;
                }
            }

            if(tours_start.size() != o.tours_start.size())
                return false;

            auto it = tours_start.begin();
            auto ot = o.tours_start.begin();
            const auto iend = tours_start.end();
            const auto oend = o.tours_start.end();
            
            while(it != iend && ot != oend) {
                if((*it) != (*ot))
                    return false;
                
                ++it;
                ++ot;
            }
*/
            return true;
        }

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
        unsigned* tours;
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
        
        inline unsigned getCustomerIndex(unsigned x, unsigned y) {
            return N + x*customers + y;
        }

        inline unsigned getDepotIndex(unsigned x, unsigned y) {
            return x*customers + y;
        }
};

#endif