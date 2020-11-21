#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <iostream>
#include <algorithm>
#include <random>
#include <set>
#include <chrono>
#include <ctime>
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
        Individual():
            cost( 0 ),
            needs_repair(true),
            improved_called(false),
            needs_to_update_cost(true),
            tours( new unsigned[customers] ),
	        tours_start()
        {
            //cout << "           default constructor called\n";
            for(unsigned i = 0; i < customers; ++i)
                tours[i] = i;
            
            tours_start[0] = 0;
        }

        Individual(const Individual& o):
            cost( o.cost ),
            needs_repair( o.needs_repair ),
            improved_called( o.improved_called ),
            needs_to_update_cost( o.needs_to_update_cost),     
            tours( new unsigned[o.customers] ),
            tours_start(o.tours_start)
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

            needs_repair = m.needs_repair;
            improved_called = m.improved_called;
            needs_to_update_cost = m.needs_to_update_cost;

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
    	
        void swap3()   {

            const auto& tours_start = this->tours_start;

            std::uniform_int_distribution<unsigned> route(0, tours_start.size()-1);

            auto it = next(tours_start.begin(), route(mt));
            const unsigned route_c = it->first; 
            ++it;
            const unsigned route_e = it == tours_start.end() ? customers : it->first;

            if(route_e - route_c >= 3) {
            
                std::uniform_int_distribution<unsigned> choice(route_c, route_e-3);
                const unsigned fst = choice(mt);
                
                std::uniform_int_distribution<unsigned> choice2(fst+1, route_e-2);
                const unsigned snd = choice2(mt);

                std::uniform_int_distribution<unsigned> choice3(snd+1, route_e-1);
                const unsigned trd = choice3(mt);

                const unsigned node = tours[fst];
                tours[fst] = tours[snd];
                tours[snd] = tours[trd];
                tours[trd] = node;     
            
            } else {

                std::uniform_int_distribution<unsigned> random_cell(0, customers-3);
                const unsigned fst = random_cell(mt);

                std::uniform_int_distribution<unsigned> random_cell2(fst+1, customers-2);
                const unsigned snd = random_cell2(mt);

                std::uniform_int_distribution<unsigned> random_cell3(snd+1, customers-1);
                const unsigned trd = random_cell3(mt);

                const unsigned node = tours[fst];
                tours[fst] = tours[snd];
                tours[snd] = tours[trd];
                tours[trd] = node;
            }

            //needs_repair = true;
            //improved_called = false;
            //needs_to_update_cost = true;

            //cout << "           individual " << this << ": swap2 executed\n";
        }
		
        void swap5()   {

            std::uniform_int_distribution<unsigned> random_cell(0, customers-5);
            const unsigned fst = random_cell(mt);

            std::uniform_int_distribution<unsigned> random_cell2(fst+1, customers-4);
            const unsigned snd = random_cell2(mt);

            std::uniform_int_distribution<unsigned> random_cell3(snd+1, customers-3);
            const unsigned trd = random_cell3(mt);

            std::uniform_int_distribution<unsigned> random_cell4(trd+1, customers-2);
            const unsigned frt = random_cell4(mt);

            std::uniform_int_distribution<unsigned> random_cell5(frt+1, customers-1);
            const unsigned fft = random_cell5(mt);
            
            const unsigned node = tours[fst];
            tours[fst] = tours[snd];
            tours[snd] = tours[trd];
            tours[frt] = tours[fft];
            tours[fft] = node;

            //needs_repair = true;
            //improved_called = false;
            //needs_to_update_cost = true;
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

            //needs_repair = true;
            //improved_called = false;
            //needs_to_update_cost = true;
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

            //needs_repair = true;
            //improved_called = false;
            //needs_to_update_cost = true;

        }

        Individual one_point_cross_over(const Individual& p1, const Individual& p2)   {
            
            Individual spare_son;

            unsigned *const tours = this->tours;
		    unsigned *const stours = spare_son.tours;
            const unsigned *const tours_p1 = p1.tours;
		    const unsigned *const tours_p2 = p2.tours;
            std::uniform_int_distribution<unsigned> random_cell(0, customers-1);      

            const unsigned cutting_point = random_cell(mt);

            unsigned index_child = 0;
            for(; index_child < cutting_point; ++index_child) {
                tours[index_child] = tours_p1[index_child];
                stours[index_child] = tours_p2[index_child];
            }

            for(; index_child < customers; ++index_child) {
                tours[index_child] = tours_p2[index_child];
                stours[index_child] = tours_p1[index_child];
            }

            auto& vmap = this->tours_start;
            auto& svmap = spare_son.tours_start;
            auto& vmap_p1 = p1.tours_start;
            auto& vmap_p2 = p2.tours_start;
            const unsigned bound = min(vmap_p1.size(), vmap_p2.size());
            std::uniform_int_distribution<unsigned> random_vehicle(0, bound-1);
            const unsigned vcutting_point = random_vehicle(mt);

            vmap.clear();
            svmap.clear();
            
            auto p1b = vmap_p1.begin();
            auto vcp_1 = next(p1b, vcutting_point); 
            vmap.insert(p1b, vcp_1);
           
            auto sp1b = vmap_p2.begin();
            auto svcp_1 = next(sp1b, vcutting_point);
            svmap.insert(sp1b, svcp_1);

            auto p2b = next(vmap_p2.begin(), vcutting_point);
            auto vcp_2 = vmap_p2.end();
            vmap.insert(p2b, vcp_2);

            auto sp2b = vcp_1;
            auto svcp_2 = vmap_p1.end();
            svmap.insert(sp2b, svcp_2);
        
            //needs_repair = true;
            //improved_called = false;
            //needs_to_update_cost = true;

            return spare_son;
        }
		
        Individual two_point_cross_over(const Individual& p1, const Individual& p2) {

            Individual spare_son;

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
                
                tours[i] = toursp1[i];
                spare_son.tours[i] = toursp2[i];
            }

            for(; i < second_cutting_point; ++i) {
                
                tours[i] = toursp2[i];
                spare_son.tours[i] = toursp1[i];
            }
            
            for(; i < customers; ++i) {
                
                tours[i] = toursp1[i];
                spare_son.tours[i] = toursp2[i];
            }

            auto& ts = this->tours_start;
            auto& ts2 = spare_son.tours_start;

            const auto& tsp1 = p1.tours_start;
            auto itp1 = tsp1.begin();
            const auto endp1 = tsp1.end();

            while(itp1 != endp1) {
                
                if(itp1->first < first_cutting_point || itp1->first > second_cutting_point)
                    ts.insert(*itp1);
                else
                    ts2.insert(*itp1);
                ++itp1;
            }

            const auto& tsp2 = p2.tours_start;
            auto itp2 = tsp2.begin();
            const auto endp2 = tsp2.end();
            while(itp2 != endp2) {
                
                if((itp2->first > first_cutting_point && itp2->first < second_cutting_point))
                    ts.insert(*itp2);
                else
                    ts2.insert(*itp2);
                ++itp2;
            }

            //needs_repair = true;
            //improved_called = false;
            //needs_to_update_cost = true;
        
            return spare_son;
        }
		
        Individual uniform_cross_over(const Individual &p1, const Individual &p2)   {

            Individual spare_son;

            std::uniform_int_distribution<unsigned> random_bit(0,1);
            unsigned index[2] = {0, 0};
            const Individual*const parents[2] = {&p1, &p2};

            for(unsigned i = 0; i < customers; ++i) {
                
                const unsigned r = random_bit(mt);
                tours[i] = parents[r]->tours[i];
                spare_son.tours[i] = parents[!r]->tours[i];
            }

            const unsigned vdim = min(p1.tours_start.size(), p2.tours_start.size());

            tours_start.clear();
            spare_son.tours_start.clear();
            for(unsigned i = 0; i < vdim; ++i) {

                const unsigned rand3 = random_bit(mt);
                const auto it = parents[rand3]->tours_start.begin();
                const auto it2 = parents[!rand3]->tours_start.begin();
                tours_start.insert( *(next(it, i)) );
                spare_son.tours_start.insert( *(next(it2, i)) );
            }
        
            //needs_repair = true;
            //improved_called = false;
            //needs_to_update_cost = true;

            return spare_son;
        }

        void random_initialize()   {
            
            //euristica costruttiva: si parte da un nodo random e da lì si prende ogni volta la scelta migliore

            tours[0] = std::uniform_int_distribution<unsigned>(0,customers-1)(mt);
            tours_start[0] = best_depots[ tours[0] ];

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
                        if(best_cost - cost > std::numeric_limits<double>::epsilon()) {
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

        void random_restart() {

            //inserisco i customers in maniera sequenziale
            for(unsigned i = 0; i < customers; ++i) {
                tours[i] = i;
            }

            //randomizzo le posizioni di ogni customer
            std::shuffle(tours, tours + customers, mt);

            std::uniform_int_distribution<unsigned> dp(0, depots-1);
            std::uniform_int_distribution<unsigned> v(1, customers-vehicles-2);
            tours_start[0] = dp(mt);
            unsigned i = v(mt);
            
            //cout << "until here";

            for(; i < vehicles;) {
                tours_start[i] = dp(mt);
                const unsigned val = v(mt);
                i = val+i > customers-vehicles-2 ? i+1 : val+i;
            }

            //needs_repair = true;
            //improved_called = false;
            //needs_to_update_cost = true;
        }
        
        void local_search()   {
        
            //if(!improved_called) {

                //cout << this << endl;

                improved_called = true;

        // prima parte: riscrittura iterata tramite 2-swap (contigui)

                    auto& dt = this->distance_table;
                    auto& ts = this->tours_start;

                    bool improved_u = true;
                    unsigned g = 0;
                    //la local search viene ripetuta ogni volta che si verifica un improvement,
                    //ma senza superare un certo numero di iterazioni
                    while(improved_u && g < 50) {

                        ++g;

                        improved_u = false;
                        const int customers_n = customers - 3;
                        int not_used = 0;

                        //cout << "c1 ";
                        while(not_used < customers_n) {

                            double old_cost = 0;
                            double new_cost = 0;

                            const unsigned s1 = not_used + 1;
                            const unsigned s2 = not_used + 2;
                            evaluate(old_cost, new_cost, not_used, 0);

                            if(old_cost - new_cost > std::numeric_limits<double>::epsilon()) {

                                if(ts.find(s1) != ts.end()) {
                                    ts[s1] = best_depots[ tours[s2] ];
                                }

                                if(ts.find(s2) != ts.end()) {
                                    ts[s2] = best_depots[ tours[s1] ];
                                }

                                const unsigned node = tours[s1]; 

                                tours[s1] = tours[s2];
                                tours[s2] = node;
                                {
                                    not_used = 0;
                                    improved_u = true;
                                }
                                //cout << not_used << " " << cost << endl; 
                            } else {
                                ++not_used;
                            }
                        }
            
            // prima parte-bis: riscrittura iterata tramite 2-swap (non contigui)

                        //cout << "c2 ";
                        not_used = 0;
                        while(not_used < customers_n) {

                            bool improved = false;
                            const unsigned s1 = not_used+1;

                            for(int offset = 1; offset < customers_n-not_used; ++ offset) {

                                double old_cost = 0;
                                double new_cost = 0;

                                const unsigned s2 = not_used+2+offset;

                                evaluate(old_cost, new_cost, not_used, offset);

                                if(old_cost - new_cost > std::numeric_limits<double>::epsilon()) {

                                    if(ts.find(s1) != ts.end()) {
                                        ts[s1] = best_depots[ tours[s2] ];
                                    }

                                    if(ts.find(s2) != ts.end()) {
                                        ts[s2] = best_depots[ tours[s1] ];
                                    }

                                    const unsigned node = tours[s1]; 

                                    tours[s1] = tours[s2];
                                    tours[s2] = node;

                                    //improved = true;
                                    {
                                        improved = true;
                                        improved_u = true;
                                    }
                                    //cout << not_used << " " << offset << " " << cost << endl;
                                }
                            }

                            if(!improved)
                                ++not_used;
                            else
                                not_used = 0;
                        }

                        //valuto lo swap del primo customer
                        const unsigned first = 0;
                        //se parto da 1 non funziona, non so perchè
                        for(unsigned off = 1; off < customers-1; ++off) {
                            
                            double old_cost = 0.0;
                            double new_cost = 0.0;

                            evaluate_first(old_cost, new_cost, off);

                            if(old_cost - new_cost > std::numeric_limits<double>::epsilon()) {

                                ts[first] = best_depots[ tours[off] ];

                                if(ts.find(off) != ts.end()) {
                                    ts[off] = best_depots[ tours[first] ];
                                }

                                const unsigned node = tours[first]; 

                                tours[first] = tours[off];
                                tours[off] = node;

                                //improved = true;
                                {
                                    improved_u = true;
                                }
                                //cout << not_used << " " << offset << " " << cost << endl;
                            }
                        }
                    
                        //valuto lo swap dell'ultimo customer
                        const unsigned last = customers-1;
                        //dovrebbe fermarsi a < customers-1, ma non so perchè non funziona
                        for(unsigned off = 1; off < customers-1; ++off) {
                            
                            double old_cost = 0.0;
                            double new_cost = 0.0;

                            evaluate_last(old_cost, new_cost, off);

                            if(old_cost - new_cost > std::numeric_limits<double>::epsilon()) {

                                if(ts.find(last) != ts.end()) {
                                    ts[last] = best_depots[ tours[off] ];
                                }

                                if(ts.find(off) != ts.end()) {
                                    ts[off] = best_depots[ tours[last] ];
                                }

                                const unsigned node = tours[last]; 

                                tours[last] = tours[off];
                                tours[off] = node;

                                improved_u = true;
                            }
                        }

                        //valuta se swappare primo e ultimo
                        {
                            double old_cost = 0.0;
                            double new_cost = 0.0;

                            old_cost += getDepotCost(ts.at(0), tours[first]);
                            new_cost += getDepotCost(best_depots[ tours[last] ], tours[last]);
#ifndef BASE
                            old_cost += activation_costs[ts.at(0)];
                            new_cost += activation_costs[ best_depots[tours[last] ] ];
#endif
                            if(ts.find(1) == ts.end()) {
                                old_cost += getCustomerCost(tours[first], tours[1]);
                                new_cost += getCustomerCost(tours[last], tours[1]);
                            }

                            if(ts.find(last) == ts.end()) {
                                old_cost += getCustomerCost(tours[customers-2], tours[last]);
                                new_cost += getCustomerCost(tours[customers-2], tours[first]);
                            } else {
                                old_cost += getDepotCost(ts.at(last), tours[last]);
                                new_cost += getDepotCost(best_depots[ tours[first] ], tours[first]);
#ifndef BASE
                                old_cost += activation_costs[ts.at(last)];
                                new_cost += activation_costs[ best_depots[tours[first] ] ];
#endif                                
                            }

                            if(old_cost - new_cost > std::numeric_limits<double>::epsilon()) {

                                ts[first] = best_depots[ tours[last] ];
                                if(ts.find(last) != ts.end()) {
                                    ts[last] = best_depots[ tours[first] ];
                                }

                                const unsigned node = tours[first];
                                tours[first] = tours[last];
                                tours[last] = node;

                                improved_u = true;
                            }
                        }
                    

                    //cout << endl;

            // seconda parte: cerco di bilanciare i subtours
#ifndef BASE
                        const double *const ac = Individual::activation_costs;
#endif   
                        vector<unsigned> p_pos(ts.size(), 0);
                        vector<double> p_len(ts.size(), 0);

                        for(unsigned ind = 0; ind < ts.size(); ++ind) {

                            unsigned k = 0;
                            for(auto it = ts.begin(); it != ts.end(); ++it) {
                                
                                auto next_it = it;
                                ++next_it;

                                const unsigned start = it->first;
                                const unsigned end = next_it == ts.end() ? customers : next_it->first;

                                p_pos[k] = start;
                                p_len[k] = calculate_tour_cost(it->first, end, false);
#ifndef BASE
                                p_len[k] += ac[it->second];
#endif
                                ++k;
                            }

                            int i, j;
                            for (i = 1; i < ts.size(); i++) {
                                    double tmp = p_len[i];
                                    unsigned tmp2 = p_pos[i];
                                    for (j = i; j >= 1 && tmp > p_len[j-1]; j--) {
                                        p_len[j] = p_len[j-1];
                                        p_pos[j] = p_pos[j-1];
                                    }
                                    p_len[j] = tmp;
                                    p_pos[j] = tmp2;
                            }

                            auto f = ts.find(p_pos[0]);
                            auto s = next(f,1);

                            if(s == ts.end()) {
                                s = f = ts.find(p_pos[0]);
                                if(f != ts.begin()) {
                                    --f;
                                } else {
                                    continue;
                                }
                            }

                            if(s != ts.end()) {

                                const unsigned p_end = next(s,1) == ts.end() ? customers : next(s,1)->first;

                                double old_cost = 
                                    calculate_tour_cost(p_pos[0], s->first, false) + best_depots[ tours[ p_pos[0] ] ] + 
                                    calculate_tour_cost(s->first, p_end, false) + best_depots[ tours[s->first] ];
                                
                                unsigned best_new_start = s->first;
                                for(unsigned new_start = p_pos[0]+1; new_start < p_end; ++new_start) {
                                    
                                    double new_cost = 
                                        calculate_tour_cost(p_pos[0], new_start, false) + best_depots[ tours[p_pos[0] ] ] +
                                        calculate_tour_cost(new_start, p_end, false) + best_depots[ tours[new_start] ];

                                    if(old_cost - new_cost > std::numeric_limits<double>::epsilon() && new_start != s->first) {
                                        best_new_start = new_start;
                                        old_cost = new_cost;
                                    }
                                }

                                if(best_new_start != s->first) {

                                    improved_u = true;
                                    ts.erase(s->first);
                                    ts[best_new_start] = best_depots[ tours[best_new_start] ];
                                }
                            }
                        }
            
            // seconda parte bis: tento di aumentare il numero di subtours

                        if(ts.size() < vehicles) {
                            
                            const unsigned size = vehicles;

                            vector<unsigned> pos(size, 0);
                            vector<double> len(size, 0);

                            for(unsigned not_used = 0; 
                                ts.size() < vehicles && not_used < ts.size(); 
                                ++not_used) {
                                
                                unsigned k = 0;
                                for(auto it = ts.begin(); it != ts.end(); ++it) {
                                    
                                    auto next_it = it;
                                    ++next_it;

                                    const unsigned start = it->first;
                                    const unsigned end = next_it == ts.end() ? customers : next_it->first;

                                    pos[k] = start;
                                    len[k] = calculate_tour_cost(it->first, end, false);
#ifndef BASE
                                    len[k] += ac[it->second];
#endif
                                    ++k;
                                }

                                int i, j;
                                for (i = 1; i < ts.size(); i++) {
                                        double tmp = len[i];
                                        unsigned tmp2 = pos[i];
                                        for (j = i; j >= 1 && tmp > len[j-1]; j--) {
                                            len[j] = len[j-1];
                                            pos[j] = pos[j-1];
                                        }
                                        len[j] = tmp;
                                        pos[j] = tmp2;
                                }
                                
                                auto it_end = ++ts.find(pos[0]);
                                unsigned end = it_end == ts.end() ? customers : it_end->first;

                                unsigned new_start = best_new_start;
                                double best_new_cost = 
                                    calculate_tour_cost(new_start, end, false) + best_depots[ tours[new_start] ] +
                                    calculate_tour_cost(pos[0], new_start, false) + best_depots[ tours[ pos[0] ] ];
#ifndef BASE    
                                best_new_cost += ac[ ( ts.find(pos[0]) )->second ];
                                best_new_cost += ac[ best_depots[ tours[new_start] ] ];
#endif                                  
                                bool split_improving = false;
                                
                                while(new_start < end) {

                                    double new_cost = 
                                        calculate_tour_cost(new_start, end, false) + best_depots[ tours[new_start] ] +
                                        calculate_tour_cost(pos[0], new_start, false) + best_depots[ tours[ pos[0] ] ];
#ifndef BASE    
                                    new_cost += ac[ ( ts.find(pos[0]) )->second ];
                                    new_cost += ac[ best_depots[ tours[new_start] ] ];
#endif                          
                                    if(best_new_cost - new_cost > std::numeric_limits<double>::epsilon()) {
                                        split_improving = true;
                                        best_new_start = new_start;
                                        best_new_cost = new_cost;
                                    }
                                    ++new_start;
                                }
                                
                                if(split_improving) {
                                    improved_u = true;
                                    ts[best_new_start] = best_depots[ tours[best_new_start] ];
                                }
                            }

                        }

            // terza parte: ottimizzo i depots

                        for(unsigned i = 0; i < ts.size(); ++i) {
                            auto it = next(ts.begin(), i);
                            const unsigned dep = best_depots[ tours[it->first] ];
                            if(it->second != dep) {

                                improved_u = true;
                                ts[it->first] = dep;
                            }
                        }
                    }
                    //needs_to_update_cost = true;
                //} 
        }

        void repair()   {
            
            //if(needs_repair) {

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

                        if( ts.find(customers-1) == ts.end() ) {
                            
                            cost += getCustomerCost(tours[customers-2], toInsert);
                        } else {
                            
                            cost += getDepotCost(best_depots[ toInsert ], toInsert);
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
                        
                        } else {

                            const unsigned depot = best_depots[ toInsert ];
                            cost += getDepotCost(depot, toInsert);

                        }
                        
                        if( i < customers-1 && ts.find(i+1) == ts.end() ) {
                                
                            cost += getCustomerCost(toInsert, tours[i+1]);
                        }
                        
                        if(best_cost - cost > std::numeric_limits<double>::epsilon()) {
                            best_cost = cost;
                            replaced = i;
                        }
                    }
                       
                    tours[replaced] = toInsert;
                    
                    if(ts.find(replaced) != ts.end())
                        ts[replaced] = best_depots[ toInsert ]; 

                    found[i] = true;
                    toReplace.erase(replaced);
                    if(replaced == customers-1)
                        last = false;
                    
                }

                //usando una map per rappresentare i veicoli, l'ordine è già mantenuto
                //bisogna limitare il numero di veicoli

                if(ts.size() > vehicles) {
#ifndef BASE
                    const double *const ac = activation_costs;
#endif
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
                            len[k] += best_depots[ tours[prev_it->first] ];
#ifndef BASE
                            len[k] += ac[it->second];
#endif
                            ++k;
                        }

                        for(; k < size; ++k) {
                            pos[k] = customers;
                            len[k] = std::numeric_limits<double>::max();
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
                        tours_start[0] = best_depots[ tours[0] ];
                    }
                }

                //print_tour();
                //cout << "           individual " << this << " repaired\n";
                //cout << ")\n";
                //needs_to_update_cost = true;
           // }
        }

        void calculate_cost()   {
            
            //if(!needs_to_update_cost)
            //    return;
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

        void calculate_diversity_ratio(unsigned himself, const std::vector<Individual>& pop) {

            double sum = 0;
            //meglio se un divisore di pop.size()
            static const unsigned d = 5;

            const unsigned min_i = (pop.size()/d) * (himself / (pop.size()/d));
            unsigned max_i = min_i + (pop.size()/d);
            max_i = max_i > pop.size() ? pop.size() : max_i;

            for(unsigned i = min_i; i < himself; ++i) {
                sum += calculate_diversity(pop[i]);
            }

            for(unsigned i = himself+1; i < max_i; ++i) {
                sum += calculate_diversity(pop[i]);
            }

            sum /= (double)pop.size();

            div_ratio = sum;
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

        inline void set_normalized(double factor) {
            normalized_cost = cost / factor;
        }

        inline double get_normalized_cost() const {
            return normalized_cost;
        }

        inline bool is_feasible() const {
            return !needs_repair;
        }

        inline double get_cost() const {

            return cost;
        }

        inline double get_diversity() const {
            return div_ratio;
        }

        inline bool operator<(const Individual& o) const   {

            return (unsigned)cost < (unsigned)o.cost;
        }

        inline bool operator==(const Individual& o) const   {

            return (unsigned)cost == (unsigned)o.cost;

        }

        static void setEnv(unsigned v, unsigned c, unsigned d, double **m, distTable* t) {

            vehicles = v;
            customers = c;
            depots = d;
            
            coordinate_matrix = m;
            distance_table = t;
#ifndef BASE
            calculate_activation_costs();
#endif

            best_depots = new unsigned[c];
            for(unsigned i = 0; i < c; ++i)
                optimizeDepot(i);
        }

        static void freeEnv() {

            delete[] best_depots;
#ifndef BASE
            delete[] activation_costs;
#endif
        }

    private:
        //costo
        double cost;
        double normalized_cost;
        //diversity ratio
        double div_ratio;
        //cardinalità veicoli, customers, depots
        inline static unsigned vehicles;
        inline static unsigned customers;
        inline static unsigned depots;
        //booleani d'utilità
        bool needs_repair;
        bool improved_called;
        bool needs_to_update_cost;
#ifndef BASE    
        //costi di attivazione per ogni depot
        inline static double *activation_costs;
#endif
        //matrice delle coordinate (per calcolare le distanze)
        inline static double **coordinate_matrix;
        //tabella di hash contente le distanze node <-> customer
        inline static distTable* distance_table;
        //conserva il depot migliore per ogni cliente
        inline static unsigned* best_depots;
        //giant tour
        unsigned* tours;
        //associa ad ogni inizio subtour il suo depot
        map<unsigned, unsigned> tours_start;

        
        static inline unsigned getCustomerIndex(unsigned x, unsigned y) {

            //cout << "\n customer index: " << x << " " << y << " " << customers * depots + x*customers + y << "| ";
            return customers * depots + x*customers + y;
        }

        static inline unsigned getDepotIndex(unsigned x, unsigned y) {
            //cout << "\n depot index: " << x << " " << y << " " << x*customers + y << "| ";
            return x*customers + y;
        }

        static inline double getCustomerCost(unsigned x, unsigned y) {
            
            const unsigned index = getCustomerIndex(x, y);
            distTable& dt = *Individual::distance_table;
            if(dt.find(index) == dt.end()) {
                const unsigned specular_index = getCustomerIndex(y, x);
                double val =
                    euclidean_distance(
                        coordinate_matrix[ x + depots ][0],
                        coordinate_matrix[ x + depots ][1],
                        coordinate_matrix[ y + depots ][0],
                        coordinate_matrix[ y + depots ][1]);
                dt[index] = val;
                dt[specular_index] = val;
            }
            return dt.at(index);
        }

        static inline double getDepotCost(unsigned x, unsigned y) {
            
            const unsigned index = getDepotIndex(x, y);
            distTable& dt = *Individual::distance_table;
            if(dt.find(index) == dt.end()) {
                dt[index] =
                    euclidean_distance(
                        coordinate_matrix[ x ][0],
                        coordinate_matrix[ x ][1],
                        coordinate_matrix[ y + depots ][0],
                        coordinate_matrix[ y + depots ][1]);
            }
            return dt.at(index);
        }

        static void optimizeDepot(const unsigned node) {

            distTable& dt = *Individual::distance_table;
#ifndef BASE
            const double *const ac = Individual::activation_costs;
#endif

            const unsigned first = node;
#ifndef BASE
            double oval = ac[0] + getDepotCost(0, first );
#else
            double oval = getDepotCost(0, first);
#endif
            unsigned depot = 0;

            for(unsigned i = 1; i < depots; ++i) {
#ifndef BASE        
//si somma il costo di attivazione col costo dell'arco depot->primo customer del subtour
                const double nval = ac[i] + getDepotCost(i, first);

                //se il costo complessivo è minore, si effettua lo swap di depot
                if(nval < oval) {
                    depot = i;
                    oval = nval;
                }
#else
                
                
                const double nval = getDepotCost(i, first);

                //se il costo complessivo è minore, si effettua lo swap di depot
                if(nval < oval) {
                    depot = i;
                    oval = nval;
                }
                
#endif
            }

            best_depots[ node ] = depot;

        }

#ifndef BASE
        static void calculate_activation_costs() {

            activation_costs = new double[depots];

            //per ogni depot
            for(unsigned i = 0; i < depots; ++i) {

                double mean = 0;

                double max_cost = getDepotCost(i, 0);
                
                mean += max_cost;
                //cout << mean << " " << max_cost << "\n";
                //per ogni cliente dopo il primo
                for(unsigned j = 1; j < customers; ++j) {

                    double cost = getDepotCost(i, j);
                    
                    mean += cost;

                    //se è la distanza massima trovata finora, conservala
                    if(cost > max_cost) {
                        max_cost = cost;
                    }

                    //cout << mean << " " << max_cost << "\n";
                }

                //calcolo della media e formula del costo d'attivazione
                mean /= (double)customers;
                double result = (max_cost - mean) * (double)customers;

                //il costo d'attivazione viene salvato
                activation_costs[i] = result/(double)vehicles;

                //cout << "RESULT: " << activation_costs[i] << endl;
            }
        }
#endif

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
                    cout << "ERROR LEN < 0 ";
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

        inline double calculate_diversity(const Individual& ind) {

            double sum = 0;
            const unsigned *const tours_ind = ind.tours;
            for(unsigned i = 0; i < customers; ++i)
                sum += tours[i] != tours_ind[i] ? 1 : 0;

            const auto it = tours_start.begin();
            const auto it2 = ind.tours_start.begin();
            
            unsigned depot = it->second;
            unsigned depot2 = it2->second;
            
            auto next_it = next(it,1);
            auto next_it2 = next(it2, 1);
            
            const auto end = tours_start.end();
            const auto end2 = ind.tours_start.end();
            for(unsigned i = 0; i < customers; ++i) {
                
                sum += tours[i] != tours_ind[i] ? 1 : 0;
                sum += (depot != depot2) ? 1 : 0;
    
                if(next_it != end && next_it->first == i) {
                    depot = next_it->second;
                    ++next_it;
                }

                if(next_it2 != end && next_it2->first == i) {
                    depot2 = next_it2->second;
                    ++next_it2;
                }
            }
            
            sum /= (double)customers * 2.0;

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

                unsigned best_depot = best_depots[ tours[i] ];

                double cost = calculate_tour_cost(i-1, i, false);

                cost += getDepotCost(best_depot, i-1);

#ifndef BASE
                if( distances[i] - distances[i-1] - cost - ac[best_depot] > std::numeric_limits<double>::epsilon() ) {
#else 
                if( distances[i] - distances[i-1] - cost > std::numeric_limits<double>::epsilon() ) {
#endif
                    distances[i] = distances[i-1] + cost;
                    predecessor[i] = best_depot;
                    improved = true;
                } 
                
                for( unsigned j = i+1; j < customers; ++j ) {
                    
                    double scost = calculate_tour_cost(i-1, j, false);

#ifdef BASE
                    if( distances[j] - distances[i-1] - scost > std::numeric_limits<double>::epsilon() ) {
#else
                    if( distances[j] - distances[i-1]*(j-i) - scost - ac[best_depot] > std::numeric_limits<double>::epsilon() ) {
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
            ts[0] = best_depots[ tours[0] ];
        //#endif

            for(unsigned i = 1; ts.size() < vehicles && i < customers; ++i) {
                if(predecessor[i] < depots) {
                    tours_start[ i-1 ] = predecessor[i];
                }
            }
        }

        inline void evaluate(double& old_cost, double& new_cost, int pos, int offset) {
            
            const auto& ts = this->tours_start;
#ifndef BASE
            const double *const ac = this->activation_costs;
#endif
            //posizione del primo candidato per lo swap
            const unsigned in_broken = pos+1;
            //posizione del secondo candidato per lo swap
            const unsigned in_joined = pos+2+offset;
            //customer successivo al secondo candidato
            const unsigned other_chain = pos+3+offset;

            const unsigned chain = tours[pos];
            const unsigned candidate = tours[in_broken];
            const unsigned candidate2 = tours[in_joined];
            const unsigned chain2 = tours[other_chain];

            //ci serve la posizione dei candidati nei rispettivi subtours
            //per valutare correttamente il costo (latency)
            const auto beg = ts.begin();
            const auto end = ts.end();

            //calcoliamo la posizione del primo candidato
            auto b_end = beg;
            while(b_end != end && b_end->first <= in_broken) {
                ++b_end;
            }

            unsigned b_len = b_end == end ? customers : b_end->first;
            b_len -= in_broken;

            //calcoliamo la posizione del secondo candidato
            auto j_end = beg;
            while(j_end != end && j_end->first <= in_joined) {
                ++j_end;
            }

            unsigned j_len = j_end == end ? customers : j_end->first;
            j_len -= in_joined;

            //print_tour();
            //cout << pos << " " << offset << " " << b_len << " " << j_len << endl;

            if(ts.find(in_broken) == end) {

                old_cost += getCustomerCost(chain, candidate) * b_len;
                new_cost += getCustomerCost(chain, candidate2) * b_len;

                //cout << "if1 " << old_cost << " " << new_cost << endl;
            
            } else {

                const unsigned depot = ts.at(in_broken);
                old_cost += getDepotCost(depot, candidate) * b_len;
                new_cost += getDepotCost(best_depots[ candidate2 ], candidate2) * b_len;
#ifndef BASE
                old_cost += ac[depot];
                new_cost += ac[best_depots[candidate2] ];
#endif

                //cout << "if2 " << old_cost << " " << new_cost << endl;
            }

            if(offset > 0) {
            
                const unsigned mid1 = in_broken+1;
                if(ts.find(mid1) == end) {
                    old_cost += getCustomerCost(candidate, tours[mid1]) * (b_len-1);
                    new_cost += getCustomerCost(candidate2, tours[mid1]) * (b_len-1);

                    //cout << "if3 " << old_cost << " " << new_cost << endl;
                }

            }

            if(ts.find(in_joined) != end) {
            
                const unsigned depot = ts.at(in_joined);
                old_cost += getDepotCost(depot, candidate2) * j_len;
                new_cost += getDepotCost(best_depots[ candidate ], candidate) * j_len;
#ifndef BASE
                old_cost += ac[depot];
                new_cost += ac[best_depots[candidate] ];
#endif
                //cout << "if4 " << old_cost << " " << new_cost << endl;
            
            } else if(offset > 0) {

                const unsigned mid2 = in_joined-1;
                old_cost += getCustomerCost(tours[mid2], candidate2) * j_len;
                new_cost += getCustomerCost(tours[mid2], candidate) * j_len;      

                //cout << "if5 " << old_cost << " " << new_cost << endl;      
            }

            if(ts.find(other_chain) == end) {

                old_cost += getCustomerCost(candidate2, chain2) * (j_len-1);
                new_cost += getCustomerCost(candidate, chain2) * (j_len-1);
            
                //cout << "if6 " << old_cost << " " << new_cost << endl;
            }
        }

        inline void evaluate_first(double& old_cost, double& new_cost, unsigned new_pos) {
            
            //print_tour();

            auto& ts = this->tours_start;
#ifndef BASE
            const double *const ac = this->activation_costs;
#endif
            const unsigned first = tours[0];
            const unsigned other = tours[new_pos];

            //ci serve la posizione dei candidati nei rispettivi subtours
            //per valutare correttamente il costo (latency)
            const auto beg = ts.begin();
            const auto end = ts.end();

            auto b_end = beg; 
            ++b_end;

            unsigned f_len = b_end == end ? customers : b_end->first;

            //calcoliamo la posizione del secondo candidato
            auto j_end = beg;
            while(j_end != end && j_end->first <= new_pos) {
                ++j_end;
            }

            unsigned j_len = j_end == end ? customers : j_end->first;
            j_len -= new_pos;

            const unsigned depot = best_depots[ first ];
            old_cost += getDepotCost(depot, first) * f_len;
            new_cost += getDepotCost(best_depots[ other ], other) * f_len;
#ifndef BASE
            old_cost += ac[ depot ];
            new_cost += ac[ best_depots[other] ];
#endif

            if(ts.find(1) == ts.end()) {
                if(new_pos > 1) {
                    old_cost += getCustomerCost(first, tours[1]) * (f_len-1);
                    new_cost += getCustomerCost(other, tours[1]) * (f_len-1);
                }
            }

            if(ts.find(new_pos) == ts.end()) {
                if(new_pos > 1) {
                    const unsigned before_other = tours[new_pos-1];
                    old_cost += getCustomerCost(before_other, other) * j_len;
                    new_cost += getCustomerCost(before_other, first) * j_len;

                   // cout << "if1 " << old_cost << " " << new_cost << endl;
                }
            } else {
                const unsigned depot2 = ts.at(new_pos);
                old_cost += getDepotCost(depot2, other) * j_len;
                new_cost += getDepotCost(best_depots[ first ], first) * j_len;
#ifndef BASE
                old_cost += ac[ depot2 ];
                new_cost += ac[ best_depots[ first ] ];
#endif
                //cout << "if2 " << old_cost << " " << new_cost << endl;
            }
            
            //questo if funziona male ma non so perchè
            if(ts.find(new_pos+1) == ts.end()) {
                if(new_pos+1 < customers) {
                    const unsigned after_other = tours[new_pos+1];
                    old_cost += getCustomerCost(other, after_other) * (j_len-1);
                    new_cost += getCustomerCost(first, after_other) * (j_len-1);
                
                    //cout << "if3 " << old_cost << " " << new_cost << " " << other << " " << after_other << " " << first << endl;      
                }
            }
        }

        inline void evaluate_last(double& old_cost, double& new_cost, unsigned new_pos) {
            
            auto& ts = this->tours_start;
#ifndef BASE
            const double *const ac = this->activation_costs;
#endif
            const unsigned last = tours[customers-1];
            const unsigned other = tours[new_pos];
        
            //ci serve la posizione dei candidati nei rispettivi subtours
            //per valutare correttamente il costo (latency)
            const auto beg = ts.begin();
            const auto end = ts.end();

            //calcoliamo la posizione del secondo candidato
            auto j_end = beg;
            while(j_end != end && j_end->first <= new_pos) {
                ++j_end;
            }

            unsigned j_len = j_end == end ? customers : j_end->first;
            j_len -= new_pos;

            if(ts.find(new_pos) == ts.end()) {
                if(new_pos > 0) {
                    const unsigned before_other = tours[new_pos-1];
                    old_cost += getCustomerCost(before_other, other) * j_len;
                    new_cost += getCustomerCost(before_other, last) * j_len;
                }
            } else {
                const unsigned depot = ts.at(new_pos);
                old_cost += getDepotCost(depot, other) * j_len;
                new_cost += getDepotCost(best_depots[ last ], last) * j_len;
            }

            if(ts.find(new_pos+1) == ts.end()) {
                if(new_pos+1 < customers-1) {
                    const unsigned after_other = tours[new_pos+1];
                    old_cost += getCustomerCost(other, after_other) * (j_len-1);
                    new_cost += getCustomerCost(last, after_other) * (j_len-1);
                }
            }

            if(ts.find(customers-1) == ts.end()) {
                if(new_pos < customers-2) {
                    const unsigned before_last = tours[customers-2];
                    old_cost += getCustomerCost(before_last, last);
                    new_cost += getCustomerCost(before_last, other);
                }
            } else {
                const unsigned depot = ts.at(customers-1);
                old_cost += getCustomerCost(depot, last);
                new_cost += getCustomerCost(best_depots[ other ], other);
            }
        }
            
};

#endif