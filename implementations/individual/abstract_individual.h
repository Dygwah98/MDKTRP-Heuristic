#ifndef AB_INDIVIDUAL_H
#define AB_INDIVIDUAL_H

class AbstractIndividual {

	public:
		virtual void initialize() = 0;
		virtual void swap2() = 0;
		virtual void swap3() = 0;
		virtual void inversion() = 0;
		virtual void scrumble() = 0;
		virtual void insertion() = 0;
		virtual void insertion_repeated() = 0;
		virtual double calculate_tour_latency(const unsigned ) = 0;
		virtual void one_point_cross_over(const AbstractIndividual&, const AbstractIndividual&) = 0;
		virtual void two_point_cross_over(const AbstractIndividual&, const AbstractIndividual&) = 0;
		virtual void best_order_cross_over(const AbstractIndividual&, const AbstractIndividual&) = 0;
		virtual void position_base_cross_over(const AbstractIndividual&, const AbstractIndividual&, const AbstractIndividual&) = 0;
		virtual void uniform_cross_over(const AbstractIndividual&, const AbstractIndividual&) = 0;
		virtual void calculate_cost() = 0;
		virtual double get_cost() const = 0;
				
		virtual bool operator<(const AbstractIndividual&) = 0;
		virtual bool operator==(const AbstractIndividual&) = 0;

		AbstractIndividual() {}
		AbstractIndividual& operator=(const AbstractIndividual&o) { return *this; }
		virtual ~AbstractIndividual() {}

};

#endif