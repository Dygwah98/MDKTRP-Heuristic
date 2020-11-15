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

        void measure_time(Individual& i, void(Individual::*function)(const Individual&, const Individual&),
                    const Individual& p1, const Individual& p2 
             )  {
            
            begin = std::chrono::steady_clock::now();

            (i.*function)(p1, p2);
            
            end = std::chrono::steady_clock::now();
            
            time += (double)std::chrono::duration_cast
                <std::chrono::nanoseconds>(end - begin)
                .count() / (double)1000000000;
            
            calls += 1;
        }

        void measure_time(Individual& i, void(Individual::*function)(const Individual&, const Individual&, const Individual&),
                        const Individual& p1, const Individual& p2, const Individual& best_individual
              ) {
            
            begin = std::chrono::steady_clock::now();

            (i.*function)(p1, p2, best_individual);
            
            end = std::chrono::steady_clock::now();
            
            time += (double)std::chrono::duration_cast
                <std::chrono::nanoseconds>(end - begin)
                .count() / (double)1000000000;
            
            calls += 1;
        }
/*
        void random_initialize(Individual& i) {
            measure_time(i, &Individual::random_initialize);
        }

        void repair(Individual& i) {
            measure_time(i, &Individual::repair);
        }

        void calculate_cost(Individual& i) {
            measure_time(i, &Individual::calculate_cost);
        }

        void swap2(Individual& i) {
            measure_time(i, &Individual::swap2);
        }

        void swap3(Individual& i) {
            measure_time(i, &Individual::swap3);
        }

        void scramble(Individual& i) {
            measure_time(i, &Individual::scramble);
        }

        void inversion(Individual& i) {
            measure_time(i, &Individual::inversion);
        }

        void insertion(Individual& i) {
            measure_time(i, &Individual::insertion);
        }

        void one_point_crossover(Individual& i, const Individual& p1, const Individual& p2) {
            measure_time(i, p1, p2, &Individual::one_point_cross_over);
        }

        void position_based_crossover(Individual& i, const Individual& p1, const Individual& p2) {
            measure_time(i, p1, p2, &Individual::position_based_cross_over);
        }

        void two_point_crossover(Individual& i, const Individual& p1, const Individual& p2) {
            measure_time(i, p1, p2, &Individual::two_point_cross_over);
        }

        void best_order_crossover(Individual& i, const Individual& p1, const Individual& p2, const Individual& best_individual) {
            measure_time(i, p1, p2, best_individual, &Individual::best_order_cross_over);
        }

        void uniform_crossover(Individual& i, const Individual& p1, const Individual& p2) {
            measure_time(i, p1, p2, &Individual::uniform_cross_over);
        }

        void improvement_algorithm(Individual& i) {
            measure_time(i, &Individual::improvement_algorithm);
        }
*/
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


#endif