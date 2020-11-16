#ifndef MISC_H
#define MISC_H

#include<vector>
#include<utility>
#include"../test.h"
#include"../utils.h"
#include"individual.h"
#include <algorithm>

class keySort {

    private:
        std::vector<Individual>* a;
        std::vector<Individual>* b;
    public:
        keySort(std::vector<Individual>& _a, 
                std::vector<Individual>& _b)
        : a(&_a), b(&_b) {}

        inline bool operator()(unsigned A, unsigned B) {
            return (*a)[A].get_cost() < (*a)[B].get_cost();
        }

        inline void _swap() { std::swap(a, b); }
};

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

        double getTotalTime() {
            
            return time;
        }

        void reset() {
            time = 0;
            calls = 0;
        }
};

class Analyzer {

    public:
        Timer initialize;
        Timer crossover;
        Timer mutation;
        Timer improvement;
        Timer repair;
        Timer costs;

        void operator()() {

            cout <<"\n    PERFORMANCE ANALYSIS (seconds):\n";
            cout <<"    random_initalize: " << initialize.getTotalTime() << endl;
            cout <<"    repair: " << repair.getTotalTime() << endl;
            cout <<"    improvement algorithm: " << improvement.getTotalTime() << endl;
            cout <<"    crossover: "  << crossover.getTotalTime() << endl;
            cout <<"    mutation: " << mutation.getTotalTime() << endl;
            cout <<"    calculate_cost: " << costs.getTotalTime() << endl;
        }
};

class Metrics {

    private:
        const double popsize;
        const unsigned max_prob;
        double mean;
        double accumulator;
        unsigned repeated;
        std::uniform_int_distribution<unsigned> u;
        std::geometric_distribution<unsigned> p;

    public:
        Metrics(double p): 
            popsize(p), 
            max_prob(std::numeric_limits<unsigned>::max()),
            mean(0),
            accumulator(0),
            repeated(0),
            p()
        {
            u = std::uniform_int_distribution<unsigned>(1, popsize-3);
        }

        inline unsigned getFirstParent() {
            
            if(repeated < 500) {
                return u(mt);
            } else {
                if(repeated % 500 == 0) {
                    p = std::geometric_distribution<unsigned>( 1.0 - (repeated)*0.000475 );
                }
                
                return (p(mt)/max_prob) * (popsize-1);
            }
        }

        inline unsigned getSecondParent(unsigned first) {

            if(repeated < 500) {
                std::uniform_int_distribution<unsigned> left(1, popsize-1-first);
                return left(mt) + first;
            
            } else {
                unsigned ret = (p(mt)/max_prob) * (popsize-2-first) + 1;

                return ret + first;
            }
        }

        inline void addToMean(double x) {
            
            accumulator += x;
        }

        inline double updateMean() {
            
            mean = accumulator / popsize;
            accumulator = 0;
            
            return mean;
        }

        inline void increaseRepeated() {
            ++repeated;
        }

        inline void resetRepeated() {
            repeated = 0;
        }
}; 

#endif