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
		
        void two_point_cross_over(const Individual& p1, const Individual& p2) {

            const unsigned size = std::min(p1.tours_start.size(), p2.tours_start.size());

            unsigned first_cutting_point;
            unsigned second_cutting_point;
            if(size >= 3) {

                const unsigned bound = size / 2;
                std::uniform_int_distribution<unsigned> cpoint1(1, bound);
                std::uniform_int_distribution<unsigned> cpoint2(bound + 1, size-1);

                first_cutting_point = next(p1.tours_start.begin(), cpoint1(mt))->first;
                second_cutting_point = next(p2.tours_start.begin(), cpoint2(mt))->first;

            } else {
                
                std::uniform_int_distribution<unsigned> crossing_point1(1, customers-2);

                first_cutting_point = crossing_point1(mt);
                
                std::uniform_int_distribution<unsigned> crossing_point2(first_cutting_point, customers-1);

                second_cutting_point = crossing_point2(mt);

            }

            unsigned *const tours = this->tours;
            const unsigned *const toursp1 = p1.tours;
            const unsigned *const toursp2 = p2.tours;

            unsigned i = 0;
            for(; i < first_cutting_point; ++i) {
                
                const unsigned node = toursp1[i];
                tours[i] = node;

            }

            for(; i < second_cutting_point; ++i) {
                
                const unsigned node = toursp2[i];
                tours[i] = node;
            }
            
            for(; i < customers; ++i) {
                
                const unsigned node = toursp1[i];
                tours[i] = node;
            }

            auto& ts = this->tours_start;

            const auto& tsp1 = p1.tours_start;
            auto itp1 = tsp1.begin();
            const auto endp1 = tsp1.end();
            while(itp1 != endp1) {
                
                if(itp1->first < first_cutting_point || itp1->first > second_cutting_point)
                    ts.insert(*itp1);
                ++itp1;
            }

            const auto& tsp2 = p2.tours_start;
            auto itp2 = tsp2.begin();
            const auto endp2 = tsp2.end();
            while(itp2 != endp2
                && (itp2->first > first_cutting_point && itp2->first < second_cutting_point)) {
                
                ts.insert(*itp2);
                ++itp2;
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

        void random_initialize()   {
            
            //euristica costruttiva: si parte da un nodo random e da lì si prende ogni volta la scelta migliore

            std::uniform_int_distribution<unsigned>random_depot(0,depots-1);

            auto& dt = this->distance_table;

            tours[0] = std::uniform_int_distribution<unsigned>(0,customers-1)(mt);
            tours_start[0] = optimizeDepot(0);

            unordered_set<unsigned> visited;
            visited.insert(tours[0]);

            auto end = visited.end();
            unsigned next = 0;
            double cost;
            for(unsigned i = 0; i < customers-1; ++i) {
                
                double best_cost = std::numeric_limits<double>::max();
                
                for(unsigned j = 0; j < customers; ++j) {
                    
                    if(visited.find(j) == end) {

                        cost = getCustomerCost(tours[i], j);
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

            repair();
            splitting_algorithm();
            repair();

            //cout << "           individual " << this << " randomly initialized\n";
        }
        
        void improvement_algorithm()   {
        
        //cout << "(";
            if(!improved_called) {

                improved_called = true;

        // prima parte: variante del k-opt con nodi invece di archi

                    auto& dt = this->distance_table;
                    auto& ts = this->tours_start;

                    const unsigned customers_n = customers_n;
                    for(unsigned not_used = 0; not_used < customers_n; ++not_used) {
                        
                        double old_cost = 0;
                        double new_cost = 0;

                        const unsigned in_broken = not_used+1;
                        const unsigned in_joined = not_used+2;

                        const unsigned chain = tours[not_used];
                        const unsigned candidate = tours[in_broken];
                        const unsigned candidate2 = tours[in_joined];

                        if(ts.find(in_broken) == ts.end()) {

                            old_cost += getCustomerCost(chain, candidate);
                            new_cost += getCustomerCost(chain, candidate);

                        } else {

                            const unsigned depot = ts.at(in_broken);
                            old_cost += getDepotCost(depot, candidate);
                            new_cost += getDepotCost(depot, candidate2);
                        }

                        if(old_cost - new_cost > 0) {

                            if(ts.find(in_broken) != ts.end()) {
                                const unsigned depot = ts.at(in_broken);
                                ts.erase(in_broken);
                                ts[in_joined] = depot;
                            }

                            if(ts.find(in_joined) != ts.end()) {
                                const unsigned depot = ts.at(in_joined);
                                ts.erase(in_joined);
                                ts[in_broken] = depot;
                            }

                            const unsigned node = tours[in_broken]; 

                            tours[in_broken] = tours[in_joined];
                            tours[in_joined] = node;
                        }
                    }


        // seconda parte: ottimizzo i depots

                    for(unsigned i = 0; i < ts.size(); ++i) {
                        auto it = next(ts.begin(), i);
                        ts[it->first] = optimizeDepot(it->first);
                    }

                } 
        }

        void repair()   {
            
            if(needs_repair) {

                //cout << "           Individual: " << this << " repair(";
                //print_tour();

                needs_repair = false;
                improved_called = false;         

                std::vector<bool> found(customers, false);

                std::unordered_set<unsigned> toReplace;
                unsigned *const tours = this->tours;

                for(unsigned i = 0; i < customers; ++i) {

                    const unsigned node = tours[i];
                    if(!found[node]) {
                        found[node] = true;
                    } else {
                        toReplace.insert(i);
                    }
                }

                bool last = false;
                if(toReplace.find(customers-1) != toReplace.end()) {
                    last = true;
                    toReplace.erase(customers-1);
                }

                auto& ts = this->tours_start;
                auto& dt = this->distance_table;

                unsigned replaced = 0;
                double best_cost;
                double cost;
                for(unsigned i = 0; i < customers; ++i) 
                if(!found[i]) {

                    unsigned toInsert = i;
                    best_cost = std::numeric_limits<double>::max();
                    replaced = cost = 0;

                    if(last) {

                        if( ts.find(customers-2) == ts.end() ) {
                            
                            cost += getCustomerCost(tours[customers-2], toInsert);
                        }

                        if(cost < best_cost) {
                            best_cost = cost;
                            replaced = customers-1;
                        }
                    }

                    for(unsigned i : toReplace) {
                        
                        cost = 0;

                        if(ts.find(i) == ts.end()) {
                                
                            cost += getCustomerCost(tours[i-1], toInsert);
                        
                            if( ts.find(i+1) == ts.end() ) {
                                
                                cost += getCustomerCost(toInsert, tours[i+1]);
                            }

                        } else {

                            const unsigned depot = optimizeDepot(toInsert);
                            cost += getDepotCost(depot, toInsert);
                        }
                        
                        if(cost < best_cost) {
                            best_cost = cost;
                            replaced = i;
                        }
                    }
                       
                    tours[replaced] = toInsert;
                    found[i] = true;
                    toReplace.erase(replaced);
                    if(replaced == customers-1)
                        last = false;
                    
                }

                //usando una map per rappresentare i veicoli, l'ordine è già mantenuto
                //bisogna limitare il numero di veicoli

                if(tours_start.size() > vehicles) {

                    //auto& ts = this->tours_start;
                    const unsigned size = ts.size();

                    vector<unsigned> pos(size, 0);
                    vector<double> len(size, 0);

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
                            len[k] = calculate_tour_cost(prev_it->first, end, false);

                            ++k;
                        }

                        for(; k < size; ++k) {
                            pos[k] = customers;
                            len[k] = std::numeric_limits<double>().max();
                        }


                        int i, j;
                        for (i = 1; i < size; i++) {
                                double tmp = len[i];
                                unsigned tmp2 = pos[i];
                                for (j = i; j >= 1 && tmp < len[j-1]; j--) {
                                    len[j] = len[j-1];
                                    pos[j] = pos[j-1];
                                }
                                len[j] = tmp;
                                pos[j] = tmp2;
                        }
                        
                        ts.erase(pos[0]);
                        len[0] = std::numeric_limits<double>::max();
                        pos[0] = customers;
                    }
                    
                    if(tours_start.find(0) == tours_start.end()) {
                        tours_start.erase(tours_start.begin());
                        tours_start[0] = optimizeDepot(0);
                    }
                }

                //print_tour();
                //cout << "           individual " << this << " repaired\n";
                //cout << ")\n";
            }
        }

        void calculate_cost()   {

            //print_tour();

            //cout << "\n       calculating all subtour costs\n";

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

            this->cost = sum;

            //cout << "           individual " << this << " cost updated from: " << oldcost <<" to: " << cost <<"\n";
        }
 
        void print_tour() {

            cout << "\ngiant tour of individual: " << this << "\n";

            for(unsigned i = 0; i < customers; ++i)
                cout << tours[i] << " ";
            cout << "\n";

            cout << "{depot : subtour => cost}\n";
            
            unsigned i = 0;
            unsigned j;
            const auto end = tours_start.end();
            for(auto it = tours_start.begin(); it != end; ++it) {

                cout << it->second << " : ";

                j = it->first;
                auto it2(it);
                const unsigned max = (++it2 == end) ? customers : it2->first;
                for(; j < max; ++j) {
                    cout << tours[j] << " ";
                }

                if(it2 == end) {
                    cout << " => " << calculate_tour_cost(it->first, customers, true);

                } else {
                    cout << " => " << calculate_tour_cost(it->first, it2->first, true);
                }

                cout << endl;
                ++i;
            }
            cout << endl;
        }

        inline double get_cost() const {

            return cost;
        }

         inline bool operator<(const Individual& o) const   {

            //cout << "           individual " << this << " operator< called versus " << &o << "\n";
            //print_tour();
            //cout << cost << "\n";
            //o.print_tour();
            //cout << o.cost << "\n";

//            if(this == &o) {
//                return true;
//            }

            return (int)cost < (int)o.cost;
        }

        inline bool operator==(const Individual& o) const   {
            
            //cout << "           individual " << this << " operator== called\n";

//            if(this == &o)
//                return true;

            return (int)cost == (int)o.cost;

        }

    private:
        //cardinalità veicoli, customers, depots
        const unsigned vehicles;
        const unsigned customers;
        const unsigned depots;
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

            //cout << "\n customer index: " << x << " " << y << " " << customers * depots + x*customers + y << "| ";
            return customers * depots + x*customers + y;
        }

        inline unsigned getDepotIndex(unsigned x, unsigned y) {
            //cout << "\n depot index: " << x << " " << y << " " << x*customers + y << "| ";
            return x*customers + y;
        }

        inline double getCustomerCost(unsigned x, unsigned y) {
            
            const unsigned index = getCustomerIndex(x, y);
            if(distance_table.find(index) == distance_table.end()) {
                distance_table[index] =
                    euclidean_distance(
                        coordinate_matrix[ x + depots ][0],
                        coordinate_matrix[ x + depots ][1],
                        coordinate_matrix[ y + depots ][0],
                        coordinate_matrix[ y + depots ][1]);
            }
            return distance_table.at(index);
        }

        inline double getDepotCost(unsigned x, unsigned y) {
            
            const unsigned index = getDepotIndex(x, y);
            if(distance_table.find(index) == distance_table.end()) {
                distance_table[index] =
                    euclidean_distance(
                        coordinate_matrix[ x ][0],
                        coordinate_matrix[ x ][1],
                        coordinate_matrix[ y + depots ][0],
                        coordinate_matrix[ y + depots ][1]);
            }
            return distance_table.at(index);
        }

        unsigned optimizeDepot(const unsigned node_position) {

            auto& dt = this->distance_table;
            auto& ts = this->tours_start;
#ifndef BASE
            auto& ac = this->activation_costs;
#endif

            const unsigned first = tours[node_position];
#ifndef BASE
            double oval = ac[0] + dt.at(  getDepotIndex(0, first ) );
#else
            double oval = getDepotCost(0, first);
#endif
            unsigned depot = 0;

            for(unsigned i = 1; i < depots; ++i) {
#ifndef BASE        

                const double nval = ac[i] + dt.at( getDepotIndex(i, first) );

                //se il costo complessivo è minore, si effettua lo swap di depot
                if(nval < oval) {
                    depot = i;
                    oval = nval;
                }
#else
                //si somma il costo di attivazione col costo dell'arco depot->primo customer del subtour
                
                const double nval = getDepotCost(i, first);

                //se il costo complessivo è minore, si effettua lo swap di depot
                if(nval < oval) {
                    depot = i;
                    oval = nval;
                }
                
#endif
            }

            return depot;

        }

        double calculate_tour_cost(const unsigned start_pos, const unsigned end_pos, bool consider_depot)   {

            //cout << "\n           starting subtour calculation...\n";
            auto& ts = this->tours_start;
            auto& dc = this->distance_table;
#ifndef BASE
            const double *const ac = this->activation_costs;
#endif

            unsigned start = start_pos;
            unsigned end = end_pos;
            if(end < start)
                std::swap(start, end);
            unsigned len = end - start;

            double sum = 0;

            if(consider_depot) {
            
                if(ts.find(start) != ts.end()) {
                    const unsigned depot = ts.at(start);
#ifndef BASE
                    sum += ac[ depot ];
#endif
                    const unsigned tsindex = tours[start];
                    sum += getDepotCost(depot, tsindex)*len;

                }
            }

            --len;

            //cout << "       depot activation + first arc costs: " << sum << endl;

            //cout << sum << endl;
            ++start;
            for(; start < end; ++start) {
                if(len <= 0)
                    cout << "ERROR";
                const unsigned c1 = tours[start-1];
                const unsigned c2 = tours[start];
                sum += getCustomerCost(c1, c2)*len;

                //cout << "\n depots: " << depots << "(" << dc.at(ci) << "*" << len <<") " << sum << endl;
                --len;
                //cout << "       sum: " << sum << endl;
            }

            
            //cout << "           individual " << this << " subtour: " << start_pos << "to: " << end_pos << "calculated\n";


            return sum;
 
        }

        void splitting_algorithm() {
            
            std::vector<double> distances(customers+1, std::numeric_limits<double>::max());
            std::vector<unsigned> predecessor(customers+1, depots);

            distances[0] = 0.0;

            auto& dt = this->distance_table;
            auto& ts = this->tours_start;
#ifndef BASE
            const double *const ac = this->activation_costs;
#endif

            for( unsigned i = 1; i < customers; ++i ) {
                
                bool improved = false;

                unsigned best_depot = optimizeDepot(i);

                double cost = calculate_tour_cost(i-1, i, false);

                cost += getDepotCost(best_depot, i-1);

#ifndef BASE
                if( distances[i-1] + cost + ac[best_depot] < distances[i] ) {
#else 
                if( distances[i-1] + cost < distances[i] ) {
#endif
                    distances[i] = distances[i-1] + cost;
                    predecessor[i] = best_depot;
                    improved = true;
                } 
                
                for( unsigned j = i+1; j < customers; ++j ) {
                    
                    double scost = calculate_tour_cost(i-1, j, false);

#ifdef BASE
                    if( distances[i-1] + scost < distances[j] ) {
#else
                    if( distances[i-1]*(j-i) + scost + ac[best_depot] < distances[j] ) {
#endif 
                        distances[j] = distances[i-1]*(j-i) + scost;
                        predecessor[j] = depots;
                        improved = true;
                    }
                            
                }

                if(!improved)
                    break;
            }

            //traduzione del risultato dell'algoritmo in termini di depots
            //auto& ts = this->tours_start;
            
        //#ifndef BASE
            ts.clear();
            ts[0] = optimizeDepot(0);
        //#endif

            for(unsigned i = 1; ts.size() < vehicles && i < customers; ++i) {
                if(predecessor[i] < depots) {
                    tours_start[ i-1 ] = predecessor[i];
                }
            }
        }
};

#endif