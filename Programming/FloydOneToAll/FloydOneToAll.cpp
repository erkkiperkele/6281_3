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
vector<int> GetInitialSequence(vector<int> &initialDistance);
void PrintChrono(int &nodes, size_t qty, double &duration);
void CreateGraph();

int _nodes;
int _nodesCount;


int main(int argc, char* argv[])
{
	int mpiRank;
	int mpiSize;

	//initialize MPI, global variables and start chrono
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	double startTime = MPI_Wtime();

	vector<int> distanceMatrix = LoadFromFile();
	vector<int> sequenceMatrix = GetInitialSequence(distanceMatrix);		//TODO

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
			number = number == 0
				? -1
				: number;
			numbers.push_back(number);
			input.get();
		}
		input.close();
	}

	_nodesCount = sqrt(numbers.size());

	int i = 0;
	while (i < numbers.size())
	{
		//Replace all distances to self with 0
		numbers[i] = 0;
		i += (_nodesCount +1);
	}

	return numbers;
}

vector<int> GetInitialSequence(vector<int> &initialDistance)
{
	vector<int> initialSequence;

	int i = 0;
	while (i < initialDistance.size())
	{
		//Initial sequence filled with arrival node (direct link). If no connection between two nodes, then sequence set to -1
		int toPush = initialDistance[i] == -1
			? -1
			: i % _nodesCount;
		initialSequence.push_back(toPush);
		++i;
	}
	return initialSequence;
}