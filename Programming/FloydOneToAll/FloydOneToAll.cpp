#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <valarray> 

using namespace std;

void PrintChrono(int &nodes, size_t qty, double &duration);
void CreateGraph();

int _nodes;
valarray<int> _graph;


int main(int argc, char* argv[])
{
    int mpiRank;
    int mpiSize;
    
	//initialize MPI, global variables and start chrono
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	double startTime = MPI_Wtime();
	
    CreateGraph();

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

void CreateGraph()
{
    int _nodes = 36;
    
    int graphFlat[] =
    {
        0,1,3,0,0,3,
        1,0,5,1,0,0,
        3,5,0,2,1,0,
        0,1,2,0,4,0,
        0,0,1,4,0,5,
        3,0,0,0,5,0
    };
    
    valarray<int> graph(graphFlat, _nodes);
    _graph = graph;
    
//    int i = 0;
//    while (i < _nodes)
//    {
//        cout << "array: " << graphFlat[i] << endl;
//        cout << "valarray: " << graph[i] << endl;
//        ++i;
//    }

}

void PrintChrono(int &nodes, size_t qty, double &duration)
{
    ofstream myfile;
	myfile.open("./FloydOneToAllChrono.csv", ios::app);
	myfile << nodes << "\t" << qty << "\t" << duration << "\n";
	myfile.close();
}