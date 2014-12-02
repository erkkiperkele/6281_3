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

//TODO:
//MISMATCH!!! Output file is similar to the 2 others, but not exactly the same. Concurrency?
//Don't send the whole submatrix, but just the rows and cols...


vector<int> LoadInitialDistances();
void DivideOrUnifyMatrix(int * matrix, int matrixSize, int submatricesCount, bool isDividing);

void FloydPipeline(int * subdistance);
void PropagateRow(int * subdistance, vector<int> &row);
void PropagateCol(int * subdistance, vector<int> &col);

void PrintChrono(double &duration);
void PrintResults(int* distanceMatrix);
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
int _rowRank;
int _colRank;

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
    if (mpiRank >= p)
    {
        return 0;
    }
    
    if (_cartRank == 0)
    {
        //Divide the matrix into submatrices to send to processes
        DivideOrUnifyMatrix(&distanceMatrix[0], _pairs, p, true);
    }
    
    //Send submatrices to processes
    int *subdistance = new int[_subPairs];
    MPI_Scatter(&distanceMatrix[0], _subPairs, MPI_INT, subdistance, _subPairs, MPI_INT, 0, _mpiCommActiveProcesses);
    
    FloydPipeline(subdistance);
    
    //Gather all submatrices into one
    MPI_Gather(subdistance, _subPairs, MPI_INT, &distanceMatrix[0], _subPairs, MPI_INT, 0, _mpiCommActiveProcesses);
    
    if (_cartRank == 0)
    {
        //Reunify submatrices into an ordered one
        DivideOrUnifyMatrix(&distanceMatrix[0], _pairs, p, false);
        
        //Finalization (stop chrono, print time and results)
        double endTime = MPI_Wtime();
        double duration = endTime - startTime;
        cout << "Duration: " << duration << " seconds" << endl;
        PrintChrono(duration);
        PrintResults(&distanceMatrix[0]);
    }
    
    MPI_Finalize();
    return 0;
}

void FloydPipeline(int * subdistance)
{
    //Calculate the submatrices
    while (k < _pRows)
    {
        //Visual feedback to indicate progress every 50th node
        if (_cartRank == 0)
        {
            cout << "iteration " << k + 1 << "/" << _pRows << endl;
        }
        ////Process position in row (x) and col (y)
        int prow = _cartRank / _pRows;
        int pcol = _cartRank % _pRows;
        
        //Calculate shortest path
        //Updating the full submatrix at each sub iteration.
        int subk = 0;
        while (subk < _subNodes)
        {
            //Propagate to rows
            vector<int> receivedRowDistance(_subPairs);
            PropagateRow(&subdistance[0], receivedRowDistance);
            
            //Propagate to cols
            vector<int> receivedColDistance(_subPairs);
            PropagateCol(&subdistance[0], receivedColDistance);
            
            //Calculating the distance to intermediate node
            int col = subk;
            int currentSubNode = 0;
            while (col < _subPairs)
            {
                //calculating the distance from intermediate node to destination
                int row = subk * _subNodes;
                while (row < (_subNodes * (subk + 1)))
                {
                    int subcol = currentSubNode / _subNodes;
                    int subrow = currentSubNode % _subNodes;
                    int coordX = pcol * _subNodes + subcol;
                    int coordY = prow * _subNodes + subrow;
                    bool isSelf = coordX == coordY;
                    
                    if (!isSelf && receivedColDistance[col] > 0 && receivedRowDistance[row] > 0)
                    {
                        
                        int newPathDistance = receivedColDistance[col] + receivedRowDistance[row];
                        
                        if (newPathDistance < subdistance[currentSubNode])
                        {
                            subdistance[currentSubNode] = newPathDistance;
                        }
                        
                        //-1 means there's no connection. Therefore any other value is better
                        if (subdistance[currentSubNode] < 0)
                        {
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
        ++k;
    }
}

void PropagateRow(int * subdistance, vector<int> &receivedRowDistance)
{
    int upper = _rowRank - 1;
    int lower = _rowRank + 1;
    
    MPI_Request sendRequestNext;
    MPI_Request sendRequestPrevious;
    MPI_Status status;
    
    if (_rowRank == k)
    {
        receivedRowDistance.assign(subdistance, subdistance + _subPairs);
    }
    
    else
    {
        int sender = _rowRank > k
        ? upper
        : lower;
        MPI_Recv(&receivedRowDistance[0], _subPairs, MPI_INT, sender, 0, _mpiCommCol, &status);
    }
    
    if (_rowRank >= k && lower < _pRows)
    {
        MPI_Isend(&receivedRowDistance[0], _subPairs, MPI_INT, lower, 0, _mpiCommCol, &sendRequestNext);
    }
    if (_rowRank <= k && upper >= 0)
    {
        MPI_Isend(&receivedRowDistance[0], _subPairs, MPI_INT, upper, 0, _mpiCommCol, &sendRequestPrevious);
    }
}

//PERF: don't propagate full submatrix
void PropagateCol(int * subdistance, vector<int> &receivedColDistance)
{
    int left = _colRank - 1;
    int right = _colRank + 1;
    
    MPI_Request sendRequestNext;
    MPI_Request sendRequestPrevious;
    MPI_Status status;
    
    if (_colRank == k)
    {
        receivedColDistance.assign(subdistance, subdistance + _subPairs);
    }
    
    else
    {
        int sender = _colRank > k
        ? left
        : right;
        MPI_Recv(&receivedColDistance[0], _subPairs, MPI_INT, sender, 0, _mpiCommRow, &status);
    }
    
    //Current rank sends on both directions
    if (_colRank >= k && right < _pRows)
    {
        MPI_Isend(&receivedColDistance[0], _subPairs, MPI_INT, right, 0, _mpiCommRow, &sendRequestNext);
    }
    
    //Other ranks send only to the direction opposite to the current rank (propagating the data)
    if (_colRank <= k && left >= 0)
    {
        MPI_Isend(&receivedColDistance[0], _subPairs, MPI_INT, left, 0, _mpiCommRow, &sendRequestPrevious);
    }
}

//													0	1	2	3		0	1	4	5			0	1	2	3
//Divides the given matrix into submatrices			4	5	6	7	->	2	3	6	7		->	4	5	6	7
// and reorganize it accordingly.					8	9	10	11		8	9	12	13			8	9	10	11
//Submatrices are ordered by row (versus column)	12	13	14	15		10	11	14	15			12	13	14	15
//
void DivideOrUnifyMatrix(int * matrix, int matrixSize, int submatricesCount, bool isDividing)
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
                    int posMatrix = isDividing
                    ? j
                    : i;
                    int posTemp = isDividing
                    ? i
                    : j;
                    matrix[posMatrix] = matrixTemp[posTemp];
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

void PrintChrono(double &duration)
{
    ofstream myfile;
    myfile.open("./FloydOneToAllChrono.csv", ios::app);
    myfile << p << "\t" << _nodes << "\t" << duration << "\n";
    myfile.close();
}

void PrintResults(int* distanceMatrix)
{
    ofstream myfile;
    myfile.open("./output.txt", ios::out);
    
    //Print matrix:
    myfile << "Final" << endl;
    for (int i = 0; i < _pairs; ++i)
    {
        myfile << distanceMatrix[i] << "\t";
        if (i % _nodes == (_nodes - 1))
            myfile << endl;
    }
    
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
    _nodes = sqrt(_pairs);
    _pRows = sqrt(mpiSize);
    
    //PERF: There's a better way
    while (_nodes %  _pRows > 0)
    {
        --_pRows;
    }
    p = pow(_pRows, 2);
    
    _subPairs = _pairs / p;
    _subNodes = sqrt(_subPairs);
    
    // Remove all unnecessary ranks
    if (mpiSize - p > 0)
    {
        vector<int> toExclude;
        int i = p;
        while (i < mpiSize)
        {
            toExclude.push_back(i);
            ++i;
        }
        
        MPI_Group newGroup;
        MPI_Group_excl(world_group, toExclude.size(), &toExclude[0], &newGroup);
        
        // Create a new communicator
        MPI_Comm_create(MPI_COMM_WORLD, newGroup, &_mpiCommActiveProcesses);
        
        //Abort excluded processes
        if (mpiRank >= p)
        {
            cout << "aborting process: " << mpiRank << endl;
            MPI_Finalize();
            return;
        }
    }
    else
    {
        MPI_Comm_dup(MPI_COMM_WORLD, &_mpiCommActiveProcesses);
    }
    
    if (mpiRank == 0)
    {
        cout << "remaining active processes: " << p << " / " << mpiSize << endl;
    }
    
    //Create Cartesian Grid comm
    int pDim[2] = { _pRows, _pRows };
    int periods[2] = { 0, 0 };
    MPI_Cart_create(_mpiCommActiveProcesses, 2, pDim, periods, 1, &_mpiCommGrid);
    MPI_Comm_rank(_mpiCommGrid, &_cartRank);
    MPI_Cart_coords(_mpiCommGrid, _cartRank, 2, _pCoords);
    
    //Create Rows comm
    int keepRows[2] = { 0, 1 };
    MPI_Cart_sub(_mpiCommGrid, keepRows, &_mpiCommRow);
    MPI_Comm_rank(_mpiCommRow, &_colRank);
    
    //Create Cols comm
    int keepCols[2] = { 1, 0 };
    MPI_Cart_sub(_mpiCommGrid, keepCols, &_mpiCommCol);
    MPI_Comm_rank(_mpiCommCol, &_rowRank);
}
