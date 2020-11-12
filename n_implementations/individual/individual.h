#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <iostream>
#include <algorithm>
#include <random>
#include <set>
#include <chrono>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <deque>
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
            needs_repair(true),
            improved_called(false)
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
            needs_repair( o.needs_repair ),
            improved_called( o.improved_called )
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

            if(std::uniform_int_distribution<unsigned>(0,1)(mt)) {
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

                std::uniform_int_distribution<unsigned> random_cell(0, customers-1);

                const unsigned fst = random_cell(mt);
                const unsigned snd = random_cell(mt);

                const unsigned node = tours[fst];
                tours[fst] = tours[snd];
                tours[snd] = node;
            }

            improved_called = false;

            //cout << "           individual " << this << ": swap2 executed\n";
        }
		
        void swap3()   {

            std::uniform_int_distribution<unsigned> random_cell(0, customers-1);

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

            
            improved_called = false;
        }
		
        void inversion()   {

            std::uniform_int_distribution<unsigned> random_cell(0, customers-1);

            unsigned start = random_cell(mt);
            unsigned end = random_cell(mt);

            while(end == start) {
                end = random_cell(mt);
            }

            if(end < start)
                std::swap(start, end);

            unsigned *const tours = this->tours;
            std::reverse(tours + start, tours + end);

            
            improved_called = false;

        }
		
        void scramble()   {

            std::uniform_int_distribution<unsigned> random_cell(0, customers-1);

            unsigned start = random_cell(mt);
            unsigned end = random_cell(mt);

            while(end == start) {
                end = random_cell(mt);
            }

            if(end < start)
                std::swap(start, end);

            unsigned *const tours = this->tours;
            std::shuffle(tours + start, tours + end, mt);

            
            improved_called = false;

        }
		
        void insertion()   {

            auto& dt = this->tours_start;
            const unsigned position = std::uniform_int_distribution<unsigned>(0, customers-1)(mt);
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
            if( std::uniform_int_distribution<unsigned> (0, 1)(mt) ) {

                unsigned best_depot = dt.at(pos);
                double min_latency = calculate_tour_cost(scroll1, scroll2, true);
                bool improved = false;
                for(unsigned i = 0; i < depots; ++i) {
                    
                    dt[pos] = i;
                    
                    double latency = calculate_tour_cost(scroll1, scroll2, true);
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
                double min_latency = calculate_tour_cost(scroll1, scroll2, true);
                for(unsigned swap = scroll1; swap < scroll2; ++swap) {

                    tours[position] = tours[swap];
                    tours[swap] = node;

                    double latency = calculate_tour_cost(scroll1, scroll2, true);
                    if(latency < min_latency) {
                        min_latency = latency;
                        best_position = swap;
                    } else {
                        tours[swap] = tours[position];
                        tours[position] = node;
                    }
                }
            }

            needs_repair = true;
            improved_called = false;
        }
		
        void one_point_cross_over(const Individual& p1, const Individual& p2)   {
            
            unsigned *const tours = this->tours;
		    const unsigned *const tours_p1 = p1.tours;
		    const unsigned *const tours_p2 = p2.tours;
            std::uniform_int_distribution<unsigned> random_cell(0, customers-1);      

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
        
            needs_repair = true;
            improved_called = false;
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
      
            auto& this_ts = this->tours_start;
            this_ts.clear();
            auto& p1_ts = p1.tours_start;
            auto& p2_ts = p2.tours_start;
            const auto& p2_ts_end = p2_ts.end();
            const auto& p1_ts_end = p1_ts.end();

            unsigned i = 0;
            for(; i < first_cutting_point; ++i) {
                
                tours[i] = tours_p2[i];
           
                if(p2_ts.find(i) != p2_ts_end) {
                    this_ts[i] = p2_ts.at(i);
                }

            }

            
            for(; i < last_cutting_point; ++i) {
                
                tours[i] = tours_p1[i];
           
                if(p1_ts.find(i) != p1_ts_end) {
                    this_ts[i] = p1_ts.at(i);
                }

            }

            for(; i < customers; ++i) {
                
                tours[i] = tours_p2[i];
          
                if(p2_ts.find(i) != p2_ts_end) {
                    this_ts[i] = p2_ts.at(i);
                }

            }
*/

            auto& ts = this->tours_start;
            map<unsigned, unsigned>().swap(ts);

            ts.insert(p1.tours_start.begin(), p1.tours_start.end());
            ts.insert(p2.tours_start.begin(), p2.tours_start.end());

            if(ts.size() > vehicles) {
                auto it = next(ts.begin(), vehicles);
                ts.erase(it, ts.end());
            }

            const unsigned size = ts.size();
            if(size >= 3) {

                const unsigned bound = size > 3 ? size - 3 : 1;
                std::uniform_int_distribution<unsigned> cpoint1(1, bound);
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

            needs_repair = true;
            improved_called = false;

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

            needs_repair = true;
            improved_called = false;

        }
		
        void position_based_cross_over(const Individual&p1, const Individual&p2)   {

            unsigned *const tours = this->tours;
            const Individual *const parents[2] = {&p1, &p2};
            std::uniform_int_distribution<unsigned> random_n(0,customers-1);
            std::uniform_int_distribution<unsigned> random_bit(0,1);

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

            needs_repair = true;
            improved_called = false;
        }
		
        void uniform_cross_over(const Individual &p1, const Individual &p2)   {

            std::uniform_int_distribution<unsigned> random_bit(0,1);
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
        
            needs_repair = true;
            improved_called = false;

        }

        //per ora ignoro il problema dello splitting
        void random_initialize()   {
            
            //euristica costruttiva: si parte da un nodo random e da lì si prende ogni volta la scelta migliore

            std::uniform_int_distribution<unsigned>random_depot(0,depots-1);

            auto& dt = this->distance_table;

            tours[0] = std::uniform_int_distribution<unsigned>(0,customers-1)(mt);
            tours_start[0] = random_depot(mt);

            unordered_set<unsigned> visited;
            visited.insert(tours[0]);

            auto end = visited.end();
            unsigned next = 0;
            double cost;
            for(unsigned i = 0; i < customers-1; ++i) {
                
                double best_cost = std::numeric_limits<double>::max();
                
                for(unsigned j = 0; j < customers; ++j) {
                    
                    if(tours[i] != j && visited.find(i) == end) {

                        const unsigned index_i = getCustomerIndex(tours[i], j);
                        if(dt.find(index_i) == dt.end()) {
                            dt[index_i] = euclidean_distance(coordinate_matrix[tours[i] + depots][0], coordinate_matrix[tours[i] + depots][1],
                                                            coordinate_matrix[j + depots][0], coordinate_matrix[j + depots][1]);
                        }
                        
                        cost = dt.at(index_i);
                        if(cost < best_cost) {
                            next = j;
                            best_cost = cost;
                        }
                    
                    }
                }
    
                tours[i+1] = next;
                visited.insert(next);
                end = visited.end();
            }

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

            needs_repair = true;
            improved_called = false;

            improvement_algorithm();
            //cout << "           individual " << this << " randomly initialized\n";
        }
        
        void improvement_algorithm()   {
        
            if(!improved_called) {

                improved_called = true;

    //variant prima parte: splitting algorithm (basato sul Bellman-Ford Shortest Path algorithm)
    /*
        vector<double>distances(customers+1, std::numeric_limits<double>::max());
        vector<bool> inqueue(customers+1, false);
        std::vector<unsigned> predecessor(customers+1, customers);
        deque<int> q;

        auto& dt = this->distance_table;
    #ifndef BASE
        const double *const ac = this->activation_costs;
        const double ac_base = ac[0];
    #endif

        distances[0] = 0;
        q.push_back(1);
        inqueue[1] = true;
        while (!q.empty()) {
            if(distances[q.back()] < distances[q.front()]) {
                int u = q.back();
                q.pop_back();
                q.push_front(u);
            }
            int v = q.front();
            q.pop_front();
            inqueue[v] = false;

            const unsigned dindex = getDepotIndex(0, v-1);
            if(dt.find(dindex) == dt.end()) {
                dt[dindex] = euclidean_distance(coordinate_matrix[0][0], 
                                                coordinate_matrix[0][1],
                                                coordinate_matrix[tours[v-1] + depots][0], 
                                                coordinate_matrix[tours[v-1] + depots][1]);
            }

            double cost = calculate_tour_cost(v-1, v, false);

            cost += dt.at(dindex);

        #ifndef BASE
            if( distances[v-1] + cost + ac_base < distances[v] ) {
        #else 
            if( distances[v-1] + cost < distances[v] ) {
        #endif
                distances[v] = distances[v-1] + cost;
                predecessor[v] = 0;
                if (!inqueue[v]) {
                    q.push_back(v);
                    inqueue[v] = true;
                }
            }

            for (unsigned j = v+1; j < customers; ++j) {
                
                double scost = calculate_tour_cost(v-1, j, false);

        #ifdef BASE
                if( distances[v-1]*(j-v+1) + scost < distances[j] ) {
        #else
                if( distances[v-1]*(j-v-1) + scost + ac_base < distances[j] ) {
        #endif 
                    distances[j] = distances[v-1]*(j-v-1) + scost;
                    predecessor[j] = v-1;
                    if (!inqueue[j]) {
                        q.push_front(j);
                        inqueue[j] = true;
                    }
                }
            }
        }
    */
    // prima parte: splitting algorithm (basato sul Bellman-Ford Shortest Path algorithm)

                std::vector<double> distances(customers+1, std::numeric_limits<double>::max());
                std::vector<unsigned> predecessor(customers+1, customers);

                distances[0] = 0.0;

                auto& dt = this->distance_table;
            #ifndef BASE
                const double *const ac = this->activation_costs;
                const double ac_base = ac[0];
            #endif

                for( unsigned i = 1; i < customers; ++i ) {
                    bool improved = false;

                    const unsigned dindex = getDepotIndex(0, i-1);
                    if(dt.find(dindex) == dt.end()) {
                        dt[dindex] = euclidean_distance(coordinate_matrix[0][0], 
                                                        coordinate_matrix[0][1],
                                                        coordinate_matrix[tours[i-1] + depots][0], 
                                                        coordinate_matrix[tours[i-1] + depots][1]);
                    }

                    double cost = calculate_tour_cost(i-1, i, false);

                    cost += dt.at(dindex);

                #ifndef BASE
                    if( distances[i-1] + cost + ac_base < distances[i] ) {
                #else 
                    if( distances[i-1] + cost < distances[i] ) {
                #endif
                        distances[i] = distances[i-1] + cost;
                        predecessor[i] = 0;
                        improved = true;
                    } 

                    for( unsigned j = i+1; j < customers; ++j ) {
                        
                        double scost = calculate_tour_cost(i-1, j, false);

                #ifdef BASE
                        if( distances[i-1] + scost < distances[j] ) {
                #else
                        if( distances[i-1]*(j-i-1) + scost + ac_base < distances[j] ) {
                #endif 
                            distances[j] = distances[i-1]*(j-i-1) + scost;
                            predecessor[j] = i-1;
                            improved = true;
                        }
                                
                    }

                    if(!improved)
                        break;
                }

                //traduzione del risultato dell'algoritmo in termini di depots
                auto& ts = this->tours_start;
                std::uniform_int_distribution<unsigned> random_depot(0,depots-1);
            #ifndef BASE
                ts.clear();
                ts[0] = random_depot(mt);
            #endif

                for(unsigned i = 1; i < customers; ++i) {
                    if(predecessor[i] == 0 && ts.size() < vehicles) {
                        tours_start[ i-1 ] = random_depot(mt);
                    }
                }


    /* seconda parte: ottimizzazione del costo del depot (sia attivazione che primo arco) */
                
                for(auto& it : ts) {
                    
                    const unsigned depot = it.second;
                    for(unsigned i = 0; i < depots; ++i) {
            #ifndef BASE        
                        
                        //si somma il costo di attivazione col costo dell'arco depot->primo customer del subtour
                        const unsigned first = tours[it.first];

                        const double nval = ac[i] + dt.at( getDepotIndex(i, first) );
                        const double oval = ac[depot] + dt.at(  getDepotIndex(depot, first) );

                        //se il costo complessivo è minore, si effettua lo swap di depot
                        if(nval < oval) {
                            it.second = i;
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
                        }
                        
            #endif
                    }
                }

            }
        }

        void repair()   {
            
            if(needs_repair) {

                needs_repair = false;
                improved_called = false;
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

                if(tours_start.size() > vehicles) {
    /*            
                    std::uniform_int_distribution<unsigned> vc(0, vehicles-1);
                    while(tours_start.size() > vehicles) {
                        tours_start.erase( next(tours_start.begin(), vc(mt)) );
                    }

                    if(tours_start.find(0) == tours_start.end()) {
                        tours_start.erase(tours_start.begin());
                        tours_start[0] = random_depot(mt);
                    }
    */

                    auto& ts = this->tours_start;
                    const unsigned size = ts.size();
                    if(size > vehicles) {

                        vector<unsigned> pos(size, 0);
                        vector<unsigned> len(size, 0);

                        while(ts.size() > vehicles) {
                            
                            unsigned k = 0;
                            for(auto it = ++ts.begin(); it != ts.end(); ++it) {
                                
                                auto next_it = it;
                                ++next_it;

                                auto prev_it = it;
                                --prev_it;

                                const unsigned start = it->first;
                                const unsigned end = next_it == ts.end() ? customers : next_it->first;

                                pos[k] = start;
                                len[k] = end - prev_it->first;

                                ++k;
                            }

                            for(; k < size; ++k) {
                                pos[k] = customers;
                                len[k] = customers;
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
                            
                            ts.erase(pos[0]);
                            len[0] = customers;
                            pos[0] = customers;
                        }

                    }        
                    
                    if(tours_start.find(0) == tours_start.end()) {
                        tours_start.erase(tours_start.begin());
                        tours_start[0] = std::uniform_int_distribution<unsigned>(0,depots-1)(mt);
                    }
                }
    /*              
                    vector<unsigned> dist(tours_start.size()-1, customers);
                    vector<unsigned> veic(tours_start.size()-1, 0);

                    unsigned index = 0;                
                    for(auto it = tours_start.begin(); next(it,2) != tours_start.end(); ++it) {
                        
                        auto next = it;
                        ++(++next);

                        dist[index] = next->first - it->first;
                        veic[index] = index;
                        
                        ++index;
                    }

                    int i, j;
                    for (i = 1; i < tours_start.size()-1; i++) {
                            unsigned tmp = dist[i];
                            unsigned tmp2 = veic[i];
                            for (j = i; j >= 1 && tmp < dist[j-1]; j--) {
                                dist[j] = dist[j-1];
                                veic[j] = veic[j-1];
                            }
                            dist[j] = tmp;
                            veic[j] = tmp2;
                    }

                    unsigned p = 0;
                    while(tours_start.size() > vehicles) {
                        tours_start.erase(veic[p]);
                    }

                    if(tours_start.find(0) == tours_start.end()) {
                        tours_start[0] = random_depot(mt);
                    }
                }
    */
                //cout << "           individual " << this << " repaired\n";
            }
        }

        double calculate_tour_cost(const unsigned start_pos, const unsigned end_pos, bool consider_depot)   {

            auto& ts = this->tours_start;
            auto& dc = this->distance_table;
        #ifndef BASE
            const double *const ac = this->activation_costs;
        #endif

            unsigned start = start_pos;
            const unsigned end = end_pos;
            unsigned len = end - start;

            double sum = 0;

            if(consider_depot) {
                
                const unsigned depot = ts.at(start_pos);
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
                
            }

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
                
                sum += calculate_tour_cost(it->first, next_it->first, true);
                //cout << "sum : " << sum << endl;

                ++it;
                ++next_it;
            }

            if(next_it == end) {

                sum += calculate_tour_cost(it->first, customers, true);
                

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
        //booleani d'utilità
        bool needs_repair;
        bool improved_called;
        
        inline unsigned getCustomerIndex(unsigned x, unsigned y) {
            return N + x*customers + y;
        }

        inline unsigned getDepotIndex(unsigned x, unsigned y) {
            return x*customers + y;
        }
};

#endif