#ifndef MISC_H
#define MISC_H

#include<vector>
#include<utility>
#include"../test.h"
#include"../utils.h"
#include"individual.h"
#include <algorithm>

//classe che ordina gli indici degli individui in base al costo (unbiased fitness ranking)
class keySort {

    private:
        std::vector<Individual>* a;
        std::vector<Individual>* b;
    public:
        keySort(std::vector<Individual>& _a, 
                std::vector<Individual>& _b)
        : a(&_a), b(&_b) {}

        inline bool operator()(unsigned A, unsigned B) {
            return ( (*a)[A].is_feasible() && !(*a)[B].is_feasible() )
                || ( (*a)[A].is_feasible() && (*a)[B].is_feasible() 
                &&   (*a)[A].get_cost() < (*a)[B].get_cost() );
        }

        inline void _swap() { std::swap(a, b); }
};

//classe che ordina gli indici degli individui in base alla misura di Vidal (biased fitness ranking)
class keyDiversitySort {

    private:
        std::vector<Individual>* a;
        std::vector<Individual>* b;
        double ratio;
        
    public:
        keyDiversitySort(std::vector<Individual>& _a, 
                    std::vector<Individual>& _b,
                    double ratio = 1.0)
        : a(&_a), b(&_b), ratio(ratio) {}

        inline double biased_rank(const Individual& ind) {
            return ind.get_normalized_cost() + ratio*ind.get_diversity();
        }

        inline bool operator()(unsigned A, unsigned B) {
            return ( (*a)[A].is_feasible() && !(*a)[B].is_feasible() )
                || ( (*a)[A].is_feasible() && (*a)[B].is_feasible() 
                &&   biased_rank((*a)[A]) < biased_rank((*a)[B])
                    );
        }

        inline void _swap() { std::swap(a, b); }
};

//classe che memorizza i tempi di esecuzione di uno o più metodi
class Timer {

    private:
        double time = 0;
        double calls = 0;
        std::chrono::steady_clock::time_point begin;
        std::chrono::steady_clock::time_point end;

    public:
        Timer() {}

        void measure_time(Individual& i, 
            void(Individual::*function)() ) {
            
            begin = std::chrono::steady_clock::now();

            (i.*function)();
            
            end = std::chrono::steady_clock::now();
            
            time += (double)std::chrono::duration_cast
                <std::chrono::nanoseconds>(end - begin)
                .count() / (double)1000000000;
            
            calls += 1;
        }

        void measure_time(Individual& i, void(Individual::*function)(unsigned, std::vector<Individual>&),
                    unsigned pos, std::vector<Individual>& population) {

            begin = std::chrono::steady_clock::now();

            (i.*function)(pos, population);
            
            end = std::chrono::steady_clock::now();
            
            time += (double)std::chrono::duration_cast
                <std::chrono::nanoseconds>(end - begin)
                .count() / (double)1000000000;
            
            calls += 1;
        }

        Individual measure_time(Individual& i, Individual(Individual::*function)(const Individual&, const Individual&),
                    const Individual& p1, const Individual& p2 
             )  {
            
            begin = std::chrono::steady_clock::now();

            Individual ind( (i.*function)(p1, p2) );
            
            end = std::chrono::steady_clock::now();
            
            time += (double)std::chrono::duration_cast
                <std::chrono::nanoseconds>(end - begin)
                .count() / (double)1000000000;
            
            calls += 1;

            return ind;
        }

        Individual measure_time(Individual& i, Individual(Individual::*function)(const Individual&, const Individual&, const Individual&),
                        const Individual& p1, const Individual& p2, const Individual& best_individual
              ) {
            
            begin = std::chrono::steady_clock::now();

            Individual ind( (i.*function)(p1, p2, best_individual) );
            
            end = std::chrono::steady_clock::now();
            
            time += (double)std::chrono::duration_cast
                <std::chrono::nanoseconds>(end - begin)
                .count() / (double)1000000000;
            
            calls += 1;

            return ind;
        }

        void measure_time(std::vector<unsigned>& D, keyDiversitySort sorts2) {
            
            begin = std::chrono::steady_clock::now();

            std::sort(D.begin(), D.end(), sorts2);
            sorts2._swap();
            
            end = std::chrono::steady_clock::now();
            
            time += (double)std::chrono::duration_cast
                <std::chrono::nanoseconds>(end - begin)
                .count() / (double)1000000000;
            
            calls += 1;

        }

        double getTotalTime() {
            
            return time;
        }

        void reset() {
            time = 0;
            calls = 0;
        }
};

//aggregatore di timer e responsabile delle stampe sulla performance dell'algoritmo
class Analyzer {

    public:
        Timer initialize;
        Timer crossover;
        Timer mutation;
        Timer improvement;
        Timer repair;
        Timer costs;
        Timer metrics;
        Timer sort;

        void operator()() {

            cout <<"\n    PERFORMANCE ANALYSIS (seconds):\n";
            cout <<"    initialize/restart: " << initialize.getTotalTime() << endl;
            cout <<"    repair: " << repair.getTotalTime() << endl;
            cout <<"    local search: " << improvement.getTotalTime() << endl;
            cout <<"    crossover: "  << crossover.getTotalTime() << endl;
            cout <<"    mutation: " << mutation.getTotalTime() << endl;
            cout <<"    calculate_cost: " << costs.getTotalTime() << endl;
            cout <<"    metrics: " << metrics.getTotalTime() << endl;
            cout <<"    sorting: " << sort.getTotalTime() << endl;
        }
};

#endif