#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

enum {SWAP2 = 0, SWAP3, SCRUMBLE, INVERSION, INSERTION};
enum {ONE_POINT = 0, TWO_POINT, BEST_ORDER, POSITION_BASED, UNIFORM, NO_CROSSOVER};


void read_first_line_cordeau(string file, int &depots, int &vehicles, int &customers)
{
	ifstream instance;
	instance.open(file);

	string line;
	getline(instance, line);
	stringstream ss(line);
	string s;
	int i = 0;
	while (getline(ss, s, ' '))
	{
		if (i == 1)
		{
			vehicles = stoi(s);
		}
		else if (i == 2)
		{
			customers = stoi(s);
		}
		else if (i == 3)
		{
			depots = stoi(s);
		}
		i++;
	}

	instance.close();
}

void load_coordinate_matrix(string file, double **matrix, const int depots, const int customers)
{
	ifstream instance;
	instance.open(file);

	string line;
	//skippa le prime depots+1 righe
	for (int i = 0; i < depots + 1; i++)
	{
		getline(instance, line);
	}

	for (int i = 0; i < depots + customers; i++)
	{
		getline(instance, line);
		stringstream ss(line);
		string s;
		int j = 0;
		while (getline(ss, s, ' '))
		{
			if (j == 1)
			{
				matrix[i][0] = stoi(s);
			}
			else if (j == 2)
			{
				matrix[i][1] = stoi(s);
			}
			j++;
		}
	}

	instance.close();
}

void swap_depots(double **matrix, const int depots, const int customers)
{
	for (int i = 0; i < depots; i++)
	{
		double *tmp = matrix[0 + i];
		matrix[0 + i] = matrix[customers + i];
		matrix[customers + i] = tmp;
	}
}

void test(int depots, int customers)
{
	// 	int N = depots + customers;
	// for (int i = 0; i < N; i++)
	// {
	// 	for (int j = 0; j < N; j++)
	// 	{
	// 		for (int k = 0; k < N; k++)
	// 		{
	// 			double d = 0;
	// 			d += 2 * euclidean_distance(coordinate_matrix[i][0], coordinate_matrix[i][1], coordinate_matrix[j][0], coordinate_matrix[j][1]);
	// 			d += euclidean_distance(coordinate_matrix[j][0], coordinate_matrix[j][1], coordinate_matrix[k][0], coordinate_matrix[k][1]);
	// 			if (d >= 23.49 && d <= 23.51)
	// 			{
	// 				cout << "TROVATA: " << d << "\n";
	// 			}
	// 		}
	// 	}
	// }
}

double **allocate_coordinate_matrix(const int rows, const int columns)
{
	double **matrix = new double *[rows];
	for (int i = 0; i < rows; i++)
	{
		matrix[i] = new double[columns];
	}
	return matrix;
}

void read_file_cordeau(string file, int &depots, int &vehicles, int &customers, double **&coordinate_matrix)
{
	read_first_line_cordeau(file, depots, vehicles, customers);

	coordinate_matrix = allocate_coordinate_matrix(depots + customers, 2);

	load_coordinate_matrix(file, coordinate_matrix, depots, customers);

	//mi interessa che i depots siano nelle prime righe della matrice e non alla fine
	swap_depots(coordinate_matrix, depots, customers);

	//print_matrix(coordinate_matrix, depots + customers, 2);
}

void read_file_ruiz(string file, unsigned &depots, unsigned &customers, double **&coordinate_matrix)
{
	ifstream instance;
	instance.open(file);

	string line;
	getline(instance, line);
	customers = stoi(line);

	getline(instance, line);
	depots = stoi(line);

	getline(instance, line);
	getline(instance, line);

	coordinate_matrix = allocate_coordinate_matrix(depots + customers, 2);

	double **customers_matrix = allocate_coordinate_matrix(customers, 2);
	for (int c = 0; c < customers; c++)
	{
		getline(instance, line);
		stringstream ss(line);
		string s;
		int j = 0;
		while (getline(ss, s, ' '))
		{
			if (j == 0)
			{
				customers_matrix[c][0] = stod(s);
			}
			else if (j == 1)
			{
				customers_matrix[c][1] = stod(s);
			}
			j++;
		}
	}

	double **depots_matrix = allocate_coordinate_matrix(depots, 2);
	for (int d = 0; d < depots; d++)
	{
		getline(instance, line);
		stringstream ss(line);
		string s;
		int j = 0;
		while (getline(ss, s, ' '))
		{
			if (j == 0)
			{
				depots_matrix[d][0] = stod(s);
			}
			else if (j == 1)
			{
				depots_matrix[d][1] = stod(s);
			}
			j++;
		}
	}

	for (int d = 0; d < depots; d++)
	{
		coordinate_matrix[d][0] = depots_matrix[d][0];
		coordinate_matrix[d][1] = depots_matrix[d][1];
	}

	for (int c = 0; c < customers; c++)
	{
		coordinate_matrix[c + depots][0] = customers_matrix[c][0];
		coordinate_matrix[c + depots][1] = customers_matrix[c][1];
	}

	for (int d = 0; d < depots; d++)
	{
		delete[] depots_matrix[d];
	}
	delete[] depots_matrix;

	for (int c = 0; c < customers; c++)
	{
		delete[] customers_matrix[c];
	}
	delete[] customers_matrix;

	instance.close();
}

inline double euclidean_distance(const int x1, const int y1, const int x2, const int y2)
{
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

void calculate_distance_matrix(double **distance_matrix, const double *const *const coordinate_matrix, const int depots, const int customers)
{
	for (int i = 0; i < depots + customers; i++)
	{
		for (int j = 0; j < depots + customers; j++)
		{
			distance_matrix[i][j] = euclidean_distance(coordinate_matrix[i][0], coordinate_matrix[i][1],
													   coordinate_matrix[j][0], coordinate_matrix[j][1]);
			//distance_matrix[i][j] = get_distance(coordinate_matrix[i][0], coordinate_matrix[i][1],
			//										   coordinate_matrix[j][0], coordinate_matrix[j][1]);

			// if (i == 1 && j == 11)
			// {
			// 	cout << "x1 =" << coordinate_matrix[i][0] << " y1 = " << coordinate_matrix[i][1] << " x2=" << coordinate_matrix[j][0] << " y2=" << coordinate_matrix[j][1] << "\n";
			// 	cout << "distanza: " << distance_matrix[i][j] << "\n";
			// }
		}
	}
}

void print_matrix(double **matrix, int R, int C)
{
	for (int i = 0; i < R; i++)
	{
		for (int j = 0; j < C; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << "\n";
	}
}

double** get_distance_matrix(const int depots, const int customers, double** coordinate_matrix)
{
    double** matrix = new double *[depots + customers];
	for (int i = 0; i < depots + customers; i++)
	{
		matrix[i] = new double[depots + customers];
	}
	calculate_distance_matrix(matrix, coordinate_matrix, depots, customers);

    return matrix;
}

void free_matrix(double **& matrix, const int rows) {
	for (int i = 0; i < rows; i++)
	{
		delete[] matrix[i];
	}
    delete[] matrix;
}

#endif