#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <iostream>
#include <algorithm>
#include <random>
#include <set>
#include <chrono>
using namespace std;

std::uint_least64_t seed(11363585407706497491);     
std::mt19937 mt(seed);

class Individual
{
public:

	Individual(const unsigned vehicles_, const unsigned depots_, const unsigned customers_, const double *const *const distance_matrix_) : vehicles(vehicles_), depots(depots_), customers(customers_), cost(0), distance_matrix(distance_matrix_)
	{
		tours = new unsigned[vehicles + customers];
		depot_positions = new unsigned[vehicles];
	}

	Individual(const Individual &o) : vehicles(o.vehicles), depots(o.depots), customers(o.customers), cost(o.cost), distance_matrix(o.distance_matrix)
	{
		tours = new unsigned[vehicles + customers];
		depot_positions = new unsigned[vehicles];

		for (unsigned i = 0; i < vehicles + customers; i++)
		{
			tours[i] = o.tours[i];
		}

		for (unsigned v = 0; v < vehicles; v++)
		{
			depot_positions[v] = o.depot_positions[v];
		}
	}

	void increase_age()
	{
		// age++;
		// if(age == 10000)
		// {
		// 	//cout<<"TROPPO VECCHIO!\n";
		// 	//cost = 999999;
		// 	// random_inizialize();
		// 	// calculate_cost();
		// 	// age = 0;
		// }
	}

	void random_inizialize()
	{
		std::uniform_int_distribution<unsigned> dist(0, depots - 1);
		unsigned *const tours = this->tours;

		//piazza casualmente i depots per i veicoli
		for (unsigned v = 0; v < vehicles; v++)
		{
			tours[v] = dist(mt);
		}

		//piazza i customer in modo sequenziale
		for (unsigned c = 0; c < customers; c++)
		{
			tours[c + vehicles] = c + depots;
		}

		std::shuffle(&tours[0], &tours[vehicles + customers], mt);

		//se il primo elemento è un customer swappalo con il primo depot che trovi
		if (tours[0] >= depots)
		{
			for (int i = 0; i < vehicles + customers; i++)
			{
				if (tours[i] < depots)
				{
					const unsigned tmp = tours[0];
					tours[0] = tours[i];
					tours[i] = tmp;
					break;
				}
			}
		}

		//indica le posizioni dei vari depot
		unsigned v = 0;
		for (unsigned i = 0; i < vehicles + customers; i++)
		{
			if (tours[i] < depots)
			{
				//start_routes_positions.insert(i);
				depot_positions[v] = i;
				v++;
			}
		}
	}

	void swap2()
	{
		std::uniform_int_distribution<unsigned> random_cell(0, vehicles + customers - 1);
		std::uniform_int_distribution<unsigned> random_depot(0, depots - 1);
		//unsigned *const tours = this->tours;

		//scegli la prima cella da swappare
		const unsigned first = random_cell(mt);
		if (tours[first] < depots)
		{
			//scegli casualmente un nuovo depot
			tours[first] = random_depot(mt);
		}
		else
		{
			//scegli una seconda cella
			unsigned second = random_cell(mt);
			//la seconda cella non può essere uguale 0 e né uguale a first
			while (second == 0 || second == first)
			{
				second = random_cell(mt);
			}

			//se il secondo è un depot aggiorna la posizione dei depot
			if (tours[second] < depots)
			{
				unsigned v = 0;
				while (depot_positions[v] != second)
					v++;
				depot_positions[v] = first;

				//elimino la posizione precedente
				//start_routes_positions.erase(second);

				//inserisco la nuova posizione
				//start_routes_positions.insert(first);
			}

			//swap tra le due celle
			const unsigned tmp = tours[first];
			tours[first] = tours[second];
			tours[second] = tmp;
		}
	}

	void swap3()
	{
		std::uniform_int_distribution<unsigned> random_cell(0, vehicles + customers - 1);
		std::uniform_int_distribution<unsigned> random_depot(0, depots - 1);
		unsigned *const tours = this->tours;

		// std::shuffle(positions,positions+vehicles+customers-1,mt);
		// const unsigned first = positions[0];
		// const unsigned second = positions[1];
		// const unsigned third = positions[2];

		//std::sort(depot_positions, depot_positions+vehicles);

		//scegli la prima cella da swappare
		const unsigned first = random_cell(mt);
		if (tours[first] < depots)
		{
			//scegli casualmente un nuovo depot
			tours[first] = random_depot(mt);
		}
		else
		{
			//scegli una seconda cella
			unsigned second = random_cell(mt);
			//la seconda cella non può essere uguale 0 e né uguale a first
			while (second == 0 || second == first)
			{
				second = random_cell(mt);
			}

			//scegli una terza cella
			unsigned third = random_cell(mt);
			//la terza cella non può essere uguale 0, né uguale a first, né a second
			while (third == 0 || third == first || third == second)
			{
				third = random_cell(mt);
			}

			//se il secondo è un depot aggiorna la posizione dei depot
			if (tours[second] < depots)
			{
				unsigned i = 0;
				while (depot_positions[i] != second)
				{
					i++;
				}
				depot_positions[i] = first;
			}

			//se il terzo è un depot aggiorna la posizione dei depot
			if (tours[third] < depots)
			{
				unsigned i = 0;
				while (depot_positions[i] != third)
					i++;
				depot_positions[i] = second;
			}

			const unsigned tmp3 = tours[third];
			tours[third] = tours[first];
			tours[first] = tours[second];
			tours[second] = tmp3;
		}
	}

	void inversion()
	{
		std::uniform_int_distribution<unsigned> random_cell(0, vehicles + customers - 1);
		std::uniform_int_distribution<unsigned> random_depot(0, depots - 1);
		unsigned *const tours = this->tours;

		//scegli la prima cella
		unsigned first = random_cell(mt);
		if (tours[first] < depots)
		{
			//scegli casualmente un nuovo depot
			tours[first] = random_depot(mt);
		}
		else
		{
			//metto anche first == 0 perché se lo estrae casualmente la seconda volta può uscire 0
			//first non può essere 0 né uguale alla fine dell'array
			while (first == vehicles + customers - 1 || first == 0)
			{
				first = random_cell(mt);
			}

			//scegli una seconda cella maggiore della prima
			std::uniform_int_distribution<unsigned> random_cell2(first + 1, vehicles + customers - 1);
			unsigned second = random_cell2(mt);
			// while (second <= first)
			// {
			// 	second = random_cell(mt);
			// }
			//inverti questo segmento di array
			std::reverse(&tours[first], &tours[second + 1]);

			//aggiorno la posizione degli altri depot
			unsigned i = first;
			for (unsigned v = 0; v < vehicles; v++)
			{
				//controllo se il depot è nel range dell'inversione
				if (depot_positions[v] >= first && depot_positions[v] <= second)
				{
					//cerco il depot all'interno dell'inversione
					while (tours[i] >= depots)
					{
						i++;
					}
					depot_positions[v] = i;
					i++;
				}
			}
		}
	}

	void scrumble()
	{
		std::uniform_int_distribution<unsigned> random_cell(0, vehicles + customers - 1);
		std::uniform_int_distribution<unsigned> random_depot(0, depots - 1);
		unsigned *const tours = this->tours;

		//scegli la prima cella
		unsigned first = random_cell(mt);
		if (tours[first] < depots)
		{
			//scegli casualmente un nuovo depot
			tours[first] = random_depot(mt);
		}
		else
		{
			while (first == vehicles + customers - 1 || first == 0)
			{
				first = random_cell(mt);
			}

			//scegli una seconda cella maggiore della prima
			std::uniform_int_distribution<unsigned> random_cell2(first + 1, vehicles + customers - 1);
			unsigned second = random_cell2(mt);
			//unsigned second = random_cell(mt);
			// while (second <= first)
			// {
			// 	second = random_cell(mt);
			// }
			std::shuffle(&tours[first], &tours[second + 1], mt);

			//aggiorno la posizione degli altri depot
			unsigned i = first;
			for (unsigned v = 0; v < vehicles; v++)
			{
				//controllo se il depot è nel range dello scrumble
				if (depot_positions[v] >= first && depot_positions[v] <= second)
				{
					//cerco il depot all'interno dello scrumble
					while (tours[i] >= depots)
					{
						i++;
					}
					depot_positions[v] = i;
					i++;
				}
			}
		}
	}

	bool insertion()
	{
		//il mutatore sceglie un elemento a caso e se è un depot sceglie il depot migliore per il tour,
		//se è un customer, sceglie la posizione migliore per il customer all'interno del suo tour
		bool improved = false;

		//scegli un elemento casualmente
		std::uniform_int_distribution<unsigned> random_position(0, vehicles + customers - 1);
		const unsigned position = random_position(mt);
		const unsigned node = tours[position];
		if (node < depots)
		{
			std::uniform_int_distribution<unsigned> random_bit(0, 1);
			if (random_bit(mt) == 0)
			{
				//prova tutti i depot e vedi quale fornisce la latenza migliore
				double min_latency = calculate_tour_latency(position);
				unsigned best_depot = node;
				for (unsigned d = 0; d < depots; d++)
				{
					if (d != node)
					{
						tours[position] = d;

						//calcola la nuova latenza
						double latency = calculate_tour_latency(position);
						if (latency < min_latency)
						{
							min_latency = latency;
							best_depot = node;

							improved = true;
						}
					}
				}

				tours[position] = best_depot;
			}
			//bisogna escludere il depot in posizione 0 perché non si può spostare
			else if (position != 0)
			{
				//cerca il depot precedente a questo
				unsigned last_tour_position = position - 1;
				while (last_tour_position >= depots)
				{
					last_tour_position--;
				}

				//cerca la posizione migliore del depot all'interno del tour
				const unsigned N = vehicles + customers;
				double min_latency = calculate_tour_latency(position) + calculate_tour_latency(last_tour_position);
				unsigned p = position + 1;
				unsigned best_position = position;
				while (p < N && tours[p] >= depots)
				{
					//swap sulla nuova posizione
					tours[position] = tours[p];
					tours[p] = node;

					double latency = calculate_tour_latency(p) + calculate_tour_latency(last_tour_position);
					if (latency < min_latency)
					{
						min_latency = latency;
						best_position = p;

						improved = true;
					}

					//reversiamo lo swap
					tours[p] = tours[position];
					tours[position] = node;

					p++;
				}

				p = position - 1;
				while (p != last_tour_position)
				{
					//swap sulla nuova posizione
					tours[position] = tours[p];
					tours[p] = node;

					double latency = calculate_tour_latency(p) + calculate_tour_latency(last_tour_position);
					if (latency < min_latency)
					{
						min_latency = latency;
						best_position = p;

						improved = true;
					}

					//reversiamo lo swap
					tours[p] = tours[position];
					tours[position] = node;

					p--;
				}

				//facciamo lo swap migliore
				tours[position] = tours[best_position];
				tours[best_position] = node;

				//aggiorniamo il depot
				unsigned v = 0;
				while (depot_positions[v] != position)
				{
					v++;
				}
				depot_positions[v] = best_position;
			}
		}
		else
		{
			//trova il depot più a sinistra
			unsigned start_tour = position - 1;
			while (tours[start_tour] >= depots)
			{
				start_tour--;
			}

			double min_latency = calculate_tour_latency(start_tour);
			unsigned p = position - 1;
			unsigned best_position = position;
			while (p != start_tour)
			{
				//swap sulla nuova posizione
				tours[position] = tours[p];
				tours[p] = node;

				double latency = calculate_tour_latency(start_tour);
				if (latency < min_latency)
				{
					min_latency = latency;
					best_position = p;

					improved = true;
				}

				//reversiamo lo swap
				tours[p] = tours[position];
				tours[position] = node;

				p--;
			}

			const unsigned N = customers + vehicles;
			p = position + 1;
			while (p < N && tours[p] >= depots)
			{
				//swap sulla nuova posizione
				tours[position] = tours[p];
				tours[p] = node;

				double latency = calculate_tour_latency(start_tour);
				if (latency < min_latency)
				{
					min_latency = latency;
					best_position = p;

					improved = true;
				}

				//reversiamo lo swap
				tours[p] = tours[position];
				tours[position] = node;

				p++;
			}

			//swap nella posizione migliore trovata
			tours[position] = tours[best_position];
			tours[best_position] = node;
		}

		return improved;
		//sanity_check();
	}

	void insertion_repeated()
	{
		//fa un tot di tentativi con insertion
		unsigned tries = 0;
		while (tries < 5 && !insertion())
		{
			tries++;
		}

		//sanity_check();
	}

	double calculate_tour_latency(const unsigned position_vehicle)
	{
		//restituisce la latenza del tour, dato il veicolo di partenza

		const unsigned N = vehicles + customers;
		//trova il prossimo depot
		unsigned p1 = position_vehicle + 1;
		unsigned len = 0;
		while (p1 < N && tours[p1] >= depots)
		{
			len++;

			p1++;
		}

		double latency = 0;
		unsigned last = position_vehicle;
		unsigned next = position_vehicle + 1;
		while (next < N && tours[next] >= depots)
		{
			latency += distance_matrix[tours[last]][tours[next]] * len;

			len--;
			last = next;
			next++;
		}

		return latency;
	}

	void one_point_cross_over(const Individual &p1, const Individual &p2)
	{
		unsigned *const tours = this->tours;
		const unsigned *const tours_p1 = p1.tours;
		const unsigned *const tours_p2 = p2.tours;

		const unsigned l = customers + depots;
		bool *customers_not_inserted = new bool[l];
		for (unsigned i = depots; i < l; i++)
		{
			customers_not_inserted[i] = true;
		}

		std::uniform_int_distribution<unsigned> random_cell(0, vehicles + customers - 1);
		//scegli cutting point
		const unsigned cutting_point = random_cell(mt);

		//copia interamente la prima metà di elementi
		unsigned index_child = 0;
		unsigned vehicles_inserted = 0;
		for (; index_child < cutting_point; index_child++)
		{
			tours[index_child] = tours_p1[index_child];
			if (tours[index_child] < depots)
			{
				depot_positions[vehicles_inserted] = index_child;
				vehicles_inserted++;
			}
			else
			{
				customers_not_inserted[tours[index_child]] = false;
			}
		}

		//inserisci il resto degli elementi in ordine di apparizione del secondo genitore
		const unsigned N = vehicles + customers;
		unsigned index_parent = 0;
		while (index_child < N)
		{
			const unsigned node = tours_p2[index_parent];
			//se è un veicolo e ancora devo inserire veicoli, lo inserisco direttamente
			if (node < depots)
			{
				if (vehicles_inserted < vehicles)
				{
					tours[index_child] = node;
					depot_positions[vehicles_inserted] = index_child;

					index_child++;
					vehicles_inserted++;
				}
			}
			//è un customer
			else
			{
				//dobbiamo capire se è stato già inserito oppure no
				if (customers_not_inserted[node])
				{
					tours[index_child] = node;
					index_child++;

					customers_not_inserted[node] = false;
				}
			}

			index_parent++;
		}

		delete[] customers_not_inserted;
		//sanity_check();
	}

	void two_point_cross_over(const Individual &p1, const Individual &p2)
	{
		unsigned *const tours = this->tours;
		const unsigned *const tours_p1 = p1.tours;
		const unsigned *const tours_p2 = p2.tours;

		const unsigned l = customers + depots;
		bool *customers_not_inserted = new bool[l];
		for (unsigned i = depots; i < l; i++)
		{
			customers_not_inserted[i] = true;
		}

		std::uniform_int_distribution<unsigned> random_cell(0, vehicles + customers - 1);
		//scegli cutting point
		unsigned cutting_point1 = random_cell(mt);
		while (cutting_point1 == vehicles + customers - 1)
		{
			cutting_point1 = random_cell(mt);
		}

		unsigned cutting_point2 = random_cell(mt);
		while (cutting_point2 <= cutting_point1)
		{
			cutting_point2 = random_cell(mt);
		}

		unsigned vehicles_inserted = 0;
		unsigned index_child = cutting_point1;
		for (; index_child < cutting_point2; index_child++)
		{
			tours[index_child] = tours_p1[index_child];
			if (tours[index_child] < depots)
			{
				depot_positions[vehicles_inserted] = index_child;
				vehicles_inserted++;
			}
			else
			{
				customers_not_inserted[tours[index_child]] = false;
			}
		}

		//inserisci il resto degli elementi in ordine di apparizione del secondo genitore
		const unsigned N = vehicles + customers;
		unsigned index_parent = 0;
		index_child = 0;
		while (index_child < cutting_point1)
		{
			const unsigned node = tours_p2[index_parent];
			//se è un veicolo e ancora devo inserire veicoli, lo inserisco direttamente
			if (node < depots)
			{
				if (vehicles_inserted < vehicles)
				{
					tours[index_child] = node;
					depot_positions[vehicles_inserted] = index_child;

					index_child++;
					vehicles_inserted++;
				}
			}
			//è un customer
			else
			{
				//dobbiamo capire se è stato già inserito oppure no
				if (customers_not_inserted[node])
				{
					tours[index_child] = node;
					index_child++;

					customers_not_inserted[node] = false;
				}
			}

			index_parent++;
		}

		//parte a destra della parte centrale
		index_child = cutting_point2;
		while (index_child < N)
		{
			const unsigned node = tours_p2[index_parent];
			//se è un veicolo e ancora devo inserire veicoli, lo inserisco direttamente
			if (node < depots)
			{
				if (vehicles_inserted < vehicles)
				{
					tours[index_child] = node;
					depot_positions[vehicles_inserted] = index_child;

					index_child++;
					vehicles_inserted++;
				}
			}
			//è un customer
			else
			{
				if (customers_not_inserted[node])
				{
					tours[index_child] = node;
					index_child++;

					customers_not_inserted[node] = false;
				}
			}

			index_parent++;
		}

		delete[] customers_not_inserted;
		//sanity_check();
	}

	void best_order_cross_over(const Individual &p1, const Individual &p2, const Individual &best)
	{
		const unsigned J = vehicles + customers;
		// int *const p1_tours = new int[J];
		// for (unsigned i = 0; i < J; i++)
		// {
		// 	p1_tours[i] = p1.tours[i];
		// }

		const unsigned l = customers + depots;
		bool *customers_in_sequence = new bool[l];
		for (unsigned i = depots; i < l; i++)
		{
			customers_in_sequence[i] = false;
		}

		//scelgo il primo cutting point
		std::uniform_int_distribution<unsigned> len_random_cutting_point(1, J / 3);
		unsigned len_1 = len_random_cutting_point(mt);
		unsigned cutting_point_1 = 0 + len_1;

		//cout << "len: " << len_1 << " J/3: " << J / 3 << "\n";

		//scegliamo una modalità per la sequenza
		std::uniform_int_distribution<unsigned> random_value_sequence(0, 2);
		const unsigned value_sequence = random_value_sequence(mt);
		unsigned v = 0;
		if (value_sequence == 0)
		{
			//copia interamente la sequenza dal genitore principale
			for (unsigned i = 0; i < cutting_point_1; i++)
			{
				const unsigned node = p1.tours[i];
				tours[i] = node;
				if (node < depots)
				{
					depot_positions[v] = i;
					v++;
				}
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
			for (unsigned i = 0; i < cutting_point_1; i++)
			{
				const unsigned node = p1.tours[i];
				if (node < depots)
				{
					vehicles_insertable++;
				}
				else
				{
					customers_in_sequence[node] = true;
				}
			}

			unsigned index_child = 0;
			for (unsigned i_p2 = 0; i_p2 < J && index_child < cutting_point_1; i_p2++)
			{
				const unsigned node = other_tours[i_p2];
				if (node < depots)
				{
					//voglio inserire un depot
					if (vehicles_insertable > 0)
					{
						vehicles_insertable--;

						tours[index_child] = node;
						depot_positions[v] = index_child;

						v++;
						index_child++;
					}
				}
				else
				{
					//voglio inserire un customer
					if (customers_in_sequence[node])
					{
						tours[index_child] = node;
						index_child++;

						customers_in_sequence[node] = false;
					}
				}
			}
		}

		unsigned last_cutting_point = cutting_point_1;
		unsigned len_before_end = J - last_cutting_point;
		while (len_before_end > J / 3)
		{
			unsigned next_cutting_point = last_cutting_point + len_random_cutting_point(mt);

			//scegli una modalità
			const unsigned value_sequence = random_value_sequence(mt);
			if (value_sequence == 0)
			{
				//copia interamente la sequenza dal genitore principale
				for (unsigned i = last_cutting_point; i < next_cutting_point; i++)
				{
					const unsigned node = p1.tours[i];
					tours[i] = node;
					if (node < depots)
					{
						depot_positions[v] = i;
						v++;
					}
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
				for (unsigned i = last_cutting_point; i < next_cutting_point; i++)
				{
					const unsigned node = p1.tours[i];
					if (node < depots)
					{
						vehicles_insertable++;
					}
					else
					{
						customers_in_sequence[node] = true;
					}
				}

				unsigned index_child = last_cutting_point;
				for (unsigned i_p2 = 0; i_p2 < J && index_child < next_cutting_point; i_p2++)
				{
					const unsigned node = other_tours[i_p2];
					if (node < depots)
					{
						//voglio inserire un depot
						if (vehicles_insertable > 0)
						{
							vehicles_insertable--;

							tours[index_child] = node;
							depot_positions[v] = index_child;

							v++;
							index_child++;
						}
					}
					else
					{
						//voglio inserire un customer
						if (customers_in_sequence[node])
						{
							tours[index_child] = node;
							index_child++;

							customers_in_sequence[node] = false;
						}
					}
				}
			}

			last_cutting_point = next_cutting_point;
			len_before_end = J - last_cutting_point;
		}

		//operiamo sull'ultimo tratto del cromosoma
		const unsigned end_sequence_value = random_value_sequence(mt);
		if (end_sequence_value == 0)
		{
			//copia interamente la sequenza dal genitore principale
			for (unsigned i = last_cutting_point; i < J; i++)
			{
				const unsigned node = p1.tours[i];
				tours[i] = node;
				if (node < depots)
				{
					depot_positions[v] = i;
					v++;
				}
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
			for (unsigned i = last_cutting_point; i < J; i++)
			{
				const unsigned node = p1.tours[i];
				if (node < depots)
				{
					vehicles_insertable++;
				}
				else
				{
					customers_in_sequence[node] = true;
				}
			}

			unsigned index_child = last_cutting_point;
			for (unsigned i_p2 = 0; i_p2 < J && index_child < J; i_p2++)
			{
				const unsigned node = other_tours[i_p2];
				if (node < depots)
				{
					//voglio inserire un depot
					if (vehicles_insertable > 0)
					{
						vehicles_insertable--;

						tours[index_child] = node;
						depot_positions[v] = index_child;

						v++;
						index_child++;
					}
				}
				else
				{
					//voglio inserire un customer
					if (customers_in_sequence[node])
					{
						tours[index_child] = node;
						index_child++;

						customers_in_sequence[node] = false;
					}
				}
			}
		}

		// p1.print_tour_matrix();
		// p2.print_tour_matrix();
		// print_tour_matrix();

		//sanity_check();

		//delete[] p1_tours;
		delete[] customers_in_sequence;
	}

	void position_base_cross_over(const Individual &p1, const Individual &p2)
	{
		const unsigned l = customers + depots;
		bool *customers_not_inserted = new bool[l];
		for (unsigned i = depots; i < l; i++)
		{
			customers_not_inserted[i] = true;
		}

		//scegli casualmente un numero di posizioni
		std::uniform_int_distribution<unsigned> random_n(1, customers + vehicles - 2);
		const unsigned n_positions = random_n(mt);

		const unsigned N = customers + vehicles - 1;

		unsigned *const positions = new unsigned[N + 1];
		for (unsigned i = 0; i <= N; i++)
		{
			positions[i] = i;
		}

		std::shuffle(positions + 1, positions + customers + vehicles, mt);

		unsigned v = 0;
		for (unsigned i = 0; i < n_positions; i++)
		{
			const unsigned position = positions[i];
			const unsigned node = p2.tours[position];

			tours[position] = node;
			if (node < depots)
			{
				depot_positions[v] = position;
				v++;
			}
			else
			{
				customers_not_inserted[node] = false;
			}
		}
		//const unsigned * const p1_tours = p1.tours;
		unsigned index_parent = 0;
		for (unsigned i = n_positions; i <= N; i++)
		{
			while (true)
			{
				const unsigned node = p1.tours[index_parent];
				if (node < depots)
				{
					if (v < vehicles)
					{
						tours[positions[i]] = node;

						depot_positions[v] = positions[i];
						v++;
						index_parent++;
						break;
					}
				}
				else
				{
					if (customers_not_inserted[node])
					{
						tours[positions[i]] = node;
						index_parent++;

						customers_not_inserted[node] = false;
						break;
					}
				}

				index_parent++;
			}
		}

		delete[] positions;
		delete[] customers_not_inserted;

		// p1.print_tour_matrix();
		// p2.print_tour_matrix();
		// print_tour_matrix();

		//sanity_check();
	}

	void uniform_cross_over(const Individual &p1, const Individual &p2)
	{
		const unsigned N = vehicles + customers;
		std::uniform_int_distribution<unsigned> random_bit(0, 1);
		unsigned index_p1 = 0;
		unsigned index_p2 = 0;
		unsigned v = 0;

		const unsigned l = customers + depots;
		bool *customers_not_inserted = new bool[l];
		for (unsigned i = depots; i < l; i++)
		{
			customers_not_inserted[i] = true;
		}

		for (unsigned i = 0; i < N; i++)
		{
			if (random_bit(mt) == 0)
			{
				//in posizione i devo mettere un valore del parent
				while (index_p1 < N)
				{
					const unsigned node = p1.tours[index_p1];
					if (node < depots)
					{
						if (v < vehicles)
						{
							tours[i] = node;
							index_p1++;

							depot_positions[v] = i;
							v++;
							break;
						}
					}
					else
					{
						if (customers_not_inserted[node])
						{
							tours[i] = node;
							index_p1++;

							customers_not_inserted[node] = false;
							break;
						}
					}

					index_p1++;
				}
			}
			else
			{
				//in posizione i devo mettere un valore del parent
				while (index_p2 < N)
				{
					const unsigned node = p2.tours[index_p2];
					if (node < depots)
					{
						if (v < vehicles)
						{
							tours[i] = node;
							index_p2++;

							depot_positions[v] = i;
							v++;
							break;
						}
					}
					else
					{
						if (customers_not_inserted[node])
						{
							tours[i] = node;
							index_p2++;

							customers_not_inserted[node] = false;
							break;
						}
					}

					index_p2++;
				}
			}
		}

		delete[] customers_not_inserted;

		// p1.print_tour_matrix();
		// p2.print_tour_matrix();
		// print_tour_matrix();

		//sanity_check();
		// cout<<"\n";
	}

	void sanity_check()
	{
		int v = 0;
		for (int i = 0; i < customers + vehicles; i++)
		{
			if (tours[i] < depots)
			{
				v++;
			}
		}
		if (v > vehicles)
		{
			cout << "TROPPI VEICOLI!!!\n";
		}

		for (int i = 0; i < vehicles; i++)
		{
			if (tours[depot_positions[i]] >= depots)
			{
				cout << "ERRORE!!!!\n";
			}
		}

		for (int i = depots; i < customers + depots; i++)
		{
			bool found = false;
			for (int j = 0; j < customers + vehicles; j++)
			{
				if (tours[j] == i)
				{
					found = true;
					break;
				}
			}
			if (!found)
			{
				cout << "CUSTOMER " << i << " NON TROVATO!!!\n";
				exit(0);
			}
		}
	}

	void calculate_cost()
	{
		const unsigned *const tours = this->tours;
		//const int *const depot_positions = this->depot_positions;
		// const std::set<unsigned> &start_routes_positions_c = start_routes_positions;

		// //cout << "ARRAY DEPOTS\n";

		// double sum = 0;

		// auto it = start_routes_positions_c.cbegin();
		// auto next_it = it;
		// ++next_it;
		// auto end = start_routes_positions_c.cend();
		// --end;
		// //calcola la latenza per le rotte tranne l'ultima
		// for (; it != end; ++it, ++next_it)
		// {
		// 	//cout<<*it<<" ";
		// 	unsigned first = *it;
		// 	unsigned second = *next_it;

		// 	int len = second - first - 1;
		// 	int last = first;
		// 	first++;
		// 	for (; first < second; first++)
		// 	{
		// 		sum += distance_matrix[tours[last]][tours[first]] * len;

		// 		len--;
		// 		last = first;
		// 	}
		// }
		// //cout<<"\n";

		// //calcolo la latenza per l'ultima rotta
		// const int N = vehicles + customers;
		// unsigned first = *it;
		// unsigned second = N;
		// int len = second - first - 1;
		// int last = first;
		// first++;
		// for (; first < second; first++)
		// {
		// 	sum += distance_matrix[tours[last]][tours[first]] * len;

		// 	len--;
		// 	last = first;
		// }

		std::sort(depot_positions, depot_positions + vehicles);

		const unsigned *const depot_positions = this->depot_positions;

		double sum = 0;
		const unsigned N = vehicles + customers;
		const unsigned vehicles_n = vehicles - 1;
		for (unsigned v = 0; v < vehicles_n; v++)
		{
			unsigned first = depot_positions[v];
			const unsigned second = depot_positions[v + 1];

			unsigned len = second - first - 1;

			unsigned last = first;
			first++;
			for (; first < second; first++)
			{
				sum += distance_matrix[tours[last]][tours[first]] * len;

				len--;
				last = first;
			}
		}

		unsigned first = depot_positions[vehicles_n];
		unsigned len = N - first - 1;
		unsigned last = first;
		first++;
		for (; first < N; first++)
		{
			sum += distance_matrix[tours[last]][tours[first]] * len;

			len--;
			last = first;
		}

		cost = sum;

		//sanity_check();
	}

	double get_cost() const
	{
		return cost;
	}

	void print_tour_matrix() const
	{
		//std::cout << "TOUR ARRAY\n";
		for (unsigned i = 0; i < vehicles + customers; i++)
		{
			std::cout << tours[i] << " ";
		}
		std::cout << "\n";
	}

	Individual &operator=(const Individual &m)
	{
		if (this == &m)
			return *this;

		cost = m.cost;
		for (unsigned i = 0; i < vehicles + customers; i++)
		{
			tours[i] = m.tours[i];
		}

		for (unsigned v = 0; v < vehicles; v++)
		{
			depot_positions[v] = m.depot_positions[v];
		}

		//start_routes_positions = m.start_routes_positions;
		//age = m.age;

		return *this;
	}

	bool operator<(const Individual &o) const
	{
		return cost < o.cost;
	}

	bool operator==(const Individual &o) const
	{
		if (this == &o)
			return true;

		const unsigned N = vehicles + customers;
		const unsigned *const tours = this->tours;
		const unsigned *const tours_o = o.tours;
		for (unsigned i = 0; i < N; i++)
		{
			if (tours[i] != tours_o[i])
			{
				return false;
			}
		}

		return true;
	}

	~Individual()
	{
		delete[] tours;
		delete[] depot_positions;
	}

private:
	const unsigned vehicles;
	const unsigned customers;
	const unsigned depots;
	double cost;
	unsigned *tours;
	const double *const *const distance_matrix;
	// è un array di n_veicoli elementi che indica dove si trovano i depot associati ai veicoli all'interno del cromosoma
	unsigned *depot_positions;
};

#endif