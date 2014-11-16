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
void DivideMatrix(int * matrix, int matrixSize, int submatricesCount);

void PrintChrono(int &nodes, size_t qty, double &duration);
void CreateGraph();

void MpiGroupInit();


int _nodesCount;		//Total number of nodes (sqrt(_pairs))
int _pairs;				//Total number of pairs (matrix size)
int _pRows;				//Number of processes row (same as cols)
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

	//Broadcast size of data
	MPI_Bcast(&_pairs, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//kill unnecessary processes. Create new group with active nodes only
	MpiGroupInit();

	if (mpiRank == 0)
	{
		//Divide both matrices into submatrices to send to processes
		DivideMatrix(&distanceMatrix[0], _pairs, p);
		DivideMatrix(&sequenceMatrix[0], _pairs, p);
	}

	//Send submatrices to processes
	int *subdistance = new int[_pairs / p];				
	int *subsequence = new int[_pairs / p];			
	MPI_Scatter(&distanceMatrix[0], _pairs / p, MPI_INT, subdistance, _pairs / p, MPI_INT, 0, _mpiCommActiveProcesses);
	MPI_Scatter(&sequenceMatrix[0], _pairs / p, MPI_INT, subsequence, _pairs / p, MPI_INT, 0, _mpiCommActiveProcesses);
	
	// COUT PROPERLY FORMATTED TO KEEP
	//	cout << "rank " << mpiRank << " - distance: " << subdistance[i] << endl;

	_nodesCount = sqrt(_pairs);
	_pRows = sqrt(p);

	//Process position in row (x) and col (y)
	int prow = mpiRank / _pRows;
	int pcol = mpiRank % _pRows;

	//Send data in columns
	if (pcol == k)
	{
		cout << "rank " << mpiRank << " - pcol " << pcol << endl;
		//Send col 0 to columnGroup
	}

	//Send data in rows
	if (prow == k)
	{
		cout << "rank " << mpiRank << " - prow " << prow << endl;
		//Send row 0 to rowGroup
	}

	MPI_Barrier(_mpiCommActiveProcesses);
	if (mpiRank == 0)
	{
		//Stop chrono and print results
		double endTime = MPI_Wtime();
		double duration = endTime - startTime;
		cout << "Duration: " << duration << " seconds" << endl;
		//		PrintChrono(mpiSize, graph.Count(), duration);
	}

	MPI_Finalize();

	return 0;
}

//																			0	1	2	3		0	1	4	5	
//Divides the given matrix into submatrices									4	5	6	7	->	2	3	6	7
// and reorganize it accordingly.											8	9	10	11		8	9	12	13
//Submatrices are ordered by row (versus column)							12	13	14	15		10	11	14	15
//
void DivideMatrix(int * matrix, int matrixSize, int submatricesCount)
{
	vector<int> matrixTemp;
	matrixTemp.assign(matrix, matrix + matrixSize);

	//all matrices are square here therefore row and col have the same size
	int subRowSize = sqrt(matrixSize / submatricesCount);		
	int rowSize = sqrt(matrixSize);
	int matricesPerRow = sqrt(submatricesCount);
	int i = 0;
	int j = 0;

	//add each row of processes to the array
	int prow = 0;
	while (prow < matricesPerRow)
	{
		//add each process of a row
		int pcol = 0;
		while (pcol < matricesPerRow)
		{
			//add each row of a process
			int subrowCount = 0;
			while (subrowCount < subRowSize)
			{
				//add each element of a process row
				int subcolCount = 0;
				while (subcolCount < subRowSize)
				{
					matrix[j] = matrixTemp[i];
					++i;
					++j;
					++subcolCount;
				}
				i += rowSize - subRowSize;
				++subrowCount;
			}
			++pcol;
			i = (subRowSize * pcol) + (prow * rowSize * subRowSize);
		}
		++prow;
		i = prow * rowSize * subRowSize;
	}
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
				? -1					//weight of 0 means there's no connection between those nodes
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
		//Replace all distances to self with 0 instead of -1
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
			: i % _nodesCount;				//arrival node id
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
		cout << "rank " << mpiRank << " - just copying MPI_WORLD" << endl;
	}
}