#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <valarray> 

using namespace std;

vector<int> LoadFromFile();
void PrintChrono(int &nodes, size_t qty, double &duration);
void CreateGraph();

int _nodes;


int main(int argc, char* argv[])
{
	int mpiRank;
	int mpiSize;

	//initialize MPI, global variables and start chrono
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	double startTime = MPI_Wtime();

	vector<int> _graph = LoadFromFile();

	if (mpiRank == 0)
	{
		//Stop chrono and print results
		double endTime = MPI_Wtime();
		double duration = endTime - startTime;
		cout << "test3" << endl;
		cout << "Duration: " << duration << " seconds" << endl;
		//		PrintChrono(mpiSize, graph.Count(), duration);
	}

	MPI_Finalize();
	return 0;
}

void PrintChrono(int &nodes, size_t qty, double &duration)
{
	ofstream myfile;
	myfile.open("./FloydOneToAllChrono.csv", ios::app);
	myfile << nodes << "\t" << qty << "\t" << duration << "\n";
	myfile.close();
}

vector<int> LoadFromFile()
{
	ifstream input("./input.txt", ios::in);

	int number;
	vector<int> numbers;

	if (input.is_open())
	{
		while (input >> number)
		{
			numbers.push_back(number);
			input.get();
		}
		input.close();
	}
	return numbers;
}