#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <valarray> 
#include <algorithm> 

using namespace std;

vector<int> LoadInitialDistances();
vector<int> GetInitialSequences(vector<int> &initialDistance);
void DivideMatrix(int * matrix, int matrixSize, int submatricesCount);

void BroadcastRow(int *receivedRowDistance);
void BroadcastCol(int *receivedColDistance);


void PrintChrono(int &nodes, size_t qty, double &duration);
void CreateGraph();

void MpiGroupInit();


int _nodes;				//Total number of nodes (sqrt(_pairs))
int _pairs;				//Total number of pairs (matrix size)
int _pRows;				//Number of processes row (same as cols)
int _subPairs;			//Number of pairs per submatrix
int _subNodes;			//Number of nodes per submatrix
int p;					//Total number of active processes
int k = 0;				//Current iteration (one iteration per row/col of processor)

MPI_Comm _mpiCommActiveProcesses;
vector<MPI_Comm> _mpiColProc;
vector<MPI_Comm> _mpiRowProc;

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
	_subPairs = _pairs / p;
	_subNodes = sqrt(_subPairs);
	int *subdistance = new int[_subPairs];
	int *subsequence = new int[_subPairs];
	MPI_Scatter(&distanceMatrix[0], _subPairs, MPI_INT, subdistance, _subPairs, MPI_INT, 0, _mpiCommActiveProcesses);

	// COUT PROPERLY FORMATTED TO KEEP
	//	cout << "rank " << mpiRank << " - distance: " << subdistance[i] << endl;

	_nodes = sqrt(_pairs);
	_pRows = sqrt(p);

	//Calculate the submatrices
	while (k < _pRows)
	{
		////Process position in row (x) and col (y)
		int prow = mpiRank / _pRows;
		int pcol = mpiRank % _pRows;

		//Broadcast in rows
		vector<int> receivedRowDistance;
		receivedRowDistance.assign(subdistance, subdistance + _subPairs);
		BroadcastRow(&receivedRowDistance[0]);

		//Broadcast in cols
		vector<int> receivedColDistance;
		receivedColDistance.assign(subdistance, subdistance + _subPairs);
		BroadcastCol(&receivedColDistance[0]);

		//Calculate shortest path
		//TODO: rows and cols and subrow and subcols are all messed up for sure
		if (pcol != k && prow != k)
		{
			int subNode = 0;
			int subk = 0;
			int currentSubNode = 0;
			//Updating the full submatrix at each sub iteration.
			while (subk < _subNodes)
			{
				//Calculating the distance to intermediate node
				int col = subk;
				while (col < _subPairs)
				{
					//calculating the distance from intermediate node to destination
					int row = subk * _subNodes;
					while (row < (_subNodes * (subk + 1)))
					{
						int subcol = subNode / _subNodes;
						int subrow = subNode % _subNodes;
						bool isSelf = (pcol * _subNodes + subcol) == (prow * _subNodes + subrow);

						if (!isSelf)
						{
							int newPathDistance = receivedColDistance[col] + receivedRowDistance[row];	//TODO: I must have swapped them.

							if (newPathDistance != subdistance[currentSubNode])
							{
								int intermediateNodeAddress = (col % _subNodes) + (pcol * _subNodes);	//Need absolute address of intermediate node
								subdistance[currentSubNode] = newPathDistance;
								subsequence[currentSubNode] = intermediateNodeAddress;
							}
						}
						++row;
						++currentSubNode;
					}
					col += _subNodes;
				}
				++subk;
			}
		}
		++k;
		MPI_Barrier(_mpiCommActiveProcesses);
	}

	//TODO: Gather all submatrices into one and reassemble it

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

//PERF: those calls are synchronized, when they should be parallel...
//Send data in rows
void BroadcastRow(int *receivedRowDistance)
{
	int i = 0;
	int col = 0;

	while (i < _pRows)
	{
		int broadcaster = col + (k * _pRows);
		MPI_Bcast(receivedRowDistance, _nodes / _pRows, MPI_INT, broadcaster, _mpiRowProc[i]);
		++col;
		++i;
	}
}

//PERF: those calls are synchronized, when they should be parallel...
//Send data in columns
void BroadcastCol(int *receivedColDistance)
{
	int j = 0;
	int row = 0;
	while (j < _pRows)
	{
		int broadcaster = (j * _pRows) + k;
		MPI_Bcast(receivedColDistance, _nodes / _pRows, MPI_INT, broadcaster, _mpiColProc[j]);
		++row;
		++j;
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
	_nodes = sqrt(_pairs);

	int i = 0;
	while (i < _pairs)
	{
		//Replace all distances to self with 0 instead of -1
		numbers[i] = 0;
		i += (_nodes + 1);
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
			: i % _nodes;				//arrival node id
		initialSequence.push_back(toPush);
		++i;
	}
	return initialSequence;
}

void MpiGroupInit()
{
	MPI_Group world_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	//TO VERIFY!!
	p = mpiSize > _pairs
		? _pairs
		: pow((int)((int)sqrt(mpiSize) - (_nodes % (int)sqrt(mpiSize))), 2);


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

		MPI_Group_excl(world_group, mpiSize - p, &toExclude[0], &newGroup);

		// Create a new communicator
		MPI_Comm_create(MPI_COMM_WORLD, newGroup, &_mpiCommActiveProcesses);
	}
	else
	{
		MPI_Comm_dup(MPI_COMM_WORLD, &_mpiCommActiveProcesses);
		cout << "rank " << mpiRank << " - just copying MPI_WORLD" << endl;
	}

	//LEGIT?
	int colSize = (int)sqrt(p);
	MPI_Comm *mpiColProc = new MPI_Comm[colSize];		//DELETE!
	MPI_Comm *mpiRowProc = new MPI_Comm[colSize];		//DELETE!

	int i = 0;
	while (i < sqrt(p))
	{
		vector<int> colP;
		vector<int> rowP;

		int col = i;
		while (col < p)
		{
			colP.push_back(col);
			col += colSize;
		}

		int row = i * colSize;
		while (row < i + colSize)
		{
			rowP.push_back(row);
			++row;
		}

		MPI_Group world_group_col;
		MPI_Comm_group(MPI_COMM_WORLD, &world_group_col);

		MPI_Group world_group_row;
		MPI_Comm_group(MPI_COMM_WORLD, &world_group_row);

		MPI_Group colGroup;
		MPI_Group_incl(world_group_col, colSize, &colP[0], &colGroup);

		MPI_Group rowGroup;
		MPI_Group_incl(world_group_row, colSize, &rowP[0], &rowGroup);

		MPI_Comm_create(_mpiCommActiveProcesses, colGroup, &mpiColProc[i]);
		MPI_Comm_create(_mpiCommActiveProcesses, rowGroup, &mpiRowProc[i]);

		++i;
	}
	_mpiColProc.assign(mpiColProc, mpiColProc + colSize);
	_mpiRowProc.assign(mpiRowProc, mpiRowProc + colSize);


	//MPI_Comm_create(_mpiCommActiveProcesses, newGroup, &col);

}