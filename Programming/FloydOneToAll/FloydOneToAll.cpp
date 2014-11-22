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
void DivideOrUnifyMatrix(int * matrix, int matrixSize, int submatricesCount);

void GetRowMatrix(int *receivedRowMatrix);
void GetColMatrix(int *receivedColMatrix);


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
MPI_Comm _mpiCommGrid;
MPI_Comm _mpiCommRow;
MPI_Comm _mpiCommCol;
int _cartRank;
int _pCoords[2];

int mpiRank;
int mpiSize;

int main(int argc, char* argv[])
{

	vector<int> distanceMatrix;

	//initialize MPI, global variables and start chrono
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	double startTime = MPI_Wtime();

	if (mpiRank == 0)
	{
		distanceMatrix = LoadInitialDistances();
	}

	//Broadcast size of data
	MPI_Bcast(&_pairs, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//kill unnecessary processes. Create new group with active nodes only
	MpiGroupInit();

	if (mpiRank == 0)
	{
		//Divide both matrices into submatrices to send to processes
		DivideOrUnifyMatrix(&distanceMatrix[0], _pairs, p);
	}

	//Send submatrices to processes
	
	int *subdistance = new int[_subPairs];
	MPI_Scatter(&distanceMatrix[0], _subPairs, MPI_INT, subdistance, _subPairs, MPI_INT, 0, _mpiCommActiveProcesses);

	// COUT PROPERLY FORMATTED TO KEEP
	//	cout << "rank " << mpiRank << " - distance: " << subdistance[i] << endl;

	//Calculate the submatrices
	while (k < _pRows)
	{
		////Process position in row (x) and col (y)
		int prow = mpiRank / _pRows;
		int pcol = mpiRank % _pRows;

		//Broadcast in rows
		int rowRank;
		MPI_Comm_rank(_mpiCommCol, &rowRank);
		vector<int> receivedRowDistance(_subPairs);

		if (rowRank == k)
		{
			cout << "rank " << mpiRank << " - init submatrix to send - rowRank: " << rowRank << endl;
			//Initialize row submatrix to broadcast
			receivedRowDistance.assign(subdistance, subdistance + _subPairs);
		}

		GetRowMatrix(&receivedRowDistance[0]);

		//Broadcast in cols
		int colRank;
		MPI_Comm_rank(_mpiCommRow, &colRank);
		vector<int> receivedColDistance(_subPairs);

		if (colRank == k)	//TEMP
		{
			cout << "rank " << mpiRank << " - init submatrix to send - colRank: " << colRank << endl;
			receivedColDistance.assign(subdistance, subdistance + _subPairs);
		}
		GetColMatrix(&receivedColDistance[0]);

		//TOREMOVE TEST ONLY: Cutting short
		MPI_Barrier(_mpiCommActiveProcesses);
		MPI_Finalize();
		return 0;
		//END OF TOREMOVE TEST ONLY: Cutting short

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

							if (newPathDistance < subdistance[currentSubNode])
							{
								//int intermediateNodeAddress = (col % _subNodes) + (pcol * _subNodes);	//Need absolute address of intermediate node
								subdistance[currentSubNode] = newPathDistance;
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

	//Gather all submatrices into one
	MPI_Gather(subdistance, _subPairs, MPI_INT, &distanceMatrix[0], _subPairs, MPI_INT, 0, _mpiCommActiveProcesses);

	//Reunify submatrices into an ordered one
	if (mpiRank == 0)
	{
		DivideOrUnifyMatrix(&distanceMatrix[0], _pairs, p);
	}

	//Finalization (stop chrono, print time etc.
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

//																			0	1	2	3		0	1	4	5			0	1	2	3
//Divides the given matrix into submatrices									4	5	6	7	->	2	3	6	7		->	4	5	6	7
// and reorganize it accordingly.											8	9	10	11		8	9	12	13			8	9	10	11
//Submatrices are ordered by row (versus column)							12	13	14	15		10	11	14	15			12	13	14	15
//
void DivideOrUnifyMatrix(int * matrix, int matrixSize, int submatricesCount)
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

void GetRowMatrix(int *receivedRowMatrix)
{
	MPI_Bcast(receivedRowMatrix, _subPairs, MPI_INT, k, _mpiCommCol);
}

void GetColMatrix(int *receivedColMatrix)
{
	MPI_Bcast(receivedColMatrix, _subPairs, MPI_INT, k, _mpiCommRow);
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
			//weight of 0 means there's no connection between those nodes
			number = number == 0
				? -1					
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


void MpiGroupInit()
{
	MPI_Group world_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	//TO VERIFY!!
	p = mpiSize > _pairs
		? _pairs
		: pow((int)((int)sqrt(mpiSize) - (_nodes % (int)sqrt(mpiSize))), 2);

	_pRows = sqrt(p);
	_subPairs = _pairs / p;
	_subNodes = sqrt(_subPairs);

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
	}

	//Create Cartesian Grid comm
	int pDim[2] = { _pRows, _pRows };
	int periods[2] = { 0, 0 };
	MPI_Cart_create(MPI_COMM_WORLD, 2, pDim, periods, 0, &_mpiCommGrid);
	MPI_Comm_rank(_mpiCommGrid, &_cartRank);
	MPI_Cart_coords(_mpiCommGrid, _cartRank, 2, _pCoords);

	//Create Rows comm
	int keepRows[2] = { 0, 1 };
	int keepCols[2] = { 1, 0 };
	MPI_Cart_sub(_mpiCommGrid, keepRows, &_mpiCommRow);

	//Create Cols comm
	MPI_Cart_sub(_mpiCommGrid, keepCols, &_mpiCommCol);
}