#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <valarray> 

using namespace std;

vector<int> LoadInitialDistances();
vector<int> GetInitialSequences(vector<int> &initialDistance);
vector<int> SubDiviseAndReorganizeMatrix(int * matrix);

void PrintChrono(int &nodes, size_t qty, double &duration);
void CreateGraph();

void MpiGroupInit();


int _nodesCount;		//Total number of nodes (sqrt(_pairs))
int _pairs;				//Total number of pairs (matrix size)
int p;					//Total number of active processes
int k = 0;				//Current iteration

MPI_Comm _mpiCommActiveProcesses;
int mpiRank;
int mpiSize;

int main(int argc, char* argv[])
{

	vector<int> distanceMatrix;
	vector<int> sequenceMatrix;

	//initialize MPI, global variables and start chrono
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	double startTime = MPI_Wtime();

	if (mpiRank == 0)
	{
		distanceMatrix = LoadInitialDistances();
		sequenceMatrix = GetInitialSequences(distanceMatrix);
	}

	MpiGroupInit();

	if (mpiRank == 0)
	{
		distanceMatrix = SubDiviseAndReorganizeMatrix(&distanceMatrix[0]);
		sequenceMatrix = SubDiviseAndReorganizeMatrix(&sequenceMatrix[0]);
	}

	//Broadcast size of data
	MPI_Bcast(&_pairs, 1, MPI_INT, 0, _mpiCommActiveProcesses);

	//Divide both matrices across processes
	int *subdistance = new int[_pairs / p];				//MUST BE A BETTER WAY
	int *subdsequence = new int[_pairs / p];			//MUST BE A BETTER WAY
	MPI_Scatter(&distanceMatrix[0], _pairs / p, MPI_INT, &subdistance, _pairs / p, MPI_INT, 0, _mpiCommActiveProcesses);
	MPI_Scatter(&sequenceMatrix[0], _pairs / p, MPI_INT, &subdsequence, _pairs / p, MPI_INT, 0, _mpiCommActiveProcesses);




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

//PERF: Return pointer only to avoid copying it multiple times?
vector<int> SubDiviseAndReorganizeMatrix(int * matrix)
{
	vector<int> matrixTemp;
	matrixTemp.assign(matrix, matrix + _pairs);

	p = 4;		//TO REMOVE: TEST ONLY
	int row = 0;
	int col = 0;
	int subRowSize = sqrt(_pairs / p);
	int rowSize = _nodesCount;
	int prowSize = sqrt(p);

	int subcol = 0;
	int subrow = 0;
	

	int i = 0;
	int j = 0;

	//add each row of processes to the array
	int prow = 0;
	while (prow < prowSize)
	{
		//add each process of a row
		int pcol = 0;
		while (pcol < prowSize)
		{
			//add each row of a process
			int subrowCount = 0;
			while (subrowCount < subRowSize)
			{
				//add each element of a process row
				int subcolCount = 0;
				while (subcolCount < subRowSize)
				{
					matrixTemp[j] = matrix[i];
					++i;
					++j;
					++subcolCount;
				}
				i += rowSize - subRowSize;
				++subrowCount;
			}
			++pcol;
			i = (subRowSize * pcol) + (prow * rowSize * subRowSize);		//6 (p7 is in col 3 ) + 14 (p7 is in prow 1   
		}
		++prow;
		i = prow * rowSize * subRowSize;
	}
	return matrixTemp;
}

void PrintChrono(int &nodes, size_t qty, double &duration)
{
	ofstream myfile;
	myfile.open("./FloydOneToAllChrono.csv", ios::app);
	myfile << nodes << "\t" << qty << "\t" << duration << "\n";
	myfile.close();
}

vector<int> LoadInitialDistances()
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

	_pairs = numbers.size();
	_nodesCount = sqrt(_pairs);

	int i = 0;
	while (i < _pairs)
	{
		//Replace all distances to self with 0
		numbers[i] = 0;
		i += (_nodesCount +1);
	}

	return numbers;
}

vector<int> GetInitialSequences(vector<int> &initialDistance)
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

void MpiGroupInit()
{
	MPI_Group world_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	p = mpiSize > _pairs
		? _pairs
		: pow((int)((int)sqrt(mpiSize) - (_nodesCount % (int)sqrt(mpiSize))), 2);


	// Remove all unnecessary ranks
	if (mpiSize - p > 0)
	{
		MPI_Group newGroup;
		
		vector<int> toExclude;
		int i = p;
		while (i < mpiSize)
		{
			toExclude.push_back(i);
			++i;
		}

		MPI_Group_excl(world_group, mpiSize-p, &toExclude[0], &newGroup);

		// Create a new communicator
		MPI_Comm_create(MPI_COMM_WORLD, newGroup, &_mpiCommActiveProcesses);
	}
	else
	{
		MPI_Comm_dup(MPI_COMM_WORLD, &_mpiCommActiveProcesses);
	}
}