#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <chrono>

using namespace std;
using namespace std::chrono;

vector<int> LoadInitialDistances();

void PrintChrono(double &duration);
void PrintResults(int* distanceMatrix);

int _nodes;				//Total number of nodes (sqrt(_pairs))
int _pairs;				//Total number of pairs (matrix size)

int main(int argc, char* argv[])
{
	high_resolution_clock::time_point startTime = high_resolution_clock::now();

	vector<int> distanceMatrix;
	distanceMatrix = LoadInitialDistances();

	int i, j, k = 0;
	int ij, ik, kj;
	for (k = 0; k < _nodes; k++)
	{
		//Visual feedback to indicate progress every 50th node
		if (k % 100 == 0)
		{
			cout << "iteration " << k << "/" << _nodes << endl;
		}

		for (i = 0; i < _nodes; i++)
		{
			for (j = 0; j < _nodes; j++)
			{
				ij = i * _nodes + j;
				ik = i * _nodes + k;
				kj = k * _nodes + j;

				// Only consider positive values since 0 means it's travelling to itself, and -1 means nodes are not connected
				if (i != j && distanceMatrix[ik] > 0 && distanceMatrix[kj] > 0)
				{
					int newPathDistance = distanceMatrix[ik] + distanceMatrix[kj];
					int oldPathDistance = distanceMatrix[ij];
					if (newPathDistance < oldPathDistance)
					{
						distanceMatrix[ij] = newPathDistance;
					}
					if (oldPathDistance == -1)
					{
						distanceMatrix[ij] = newPathDistance;
					}
				}
			}
		}
	}

	//Finalization (stop chrono, print time and results)
	high_resolution_clock::time_point endTime = high_resolution_clock::now();
	double time_span = duration_cast< duration<double> >(endTime - startTime).count();
	cout << "Duration: " << time_span << " seconds" << endl;
	PrintChrono(time_span);
	PrintResults(&distanceMatrix[0]);

	return 0;
}


void PrintChrono(double &duration)
{
	ofstream myfile;
	myfile.open("./FloydSerial.csv", ios::app);
	myfile << 1 << "\t" << _nodes << "\t" << duration << "\n";
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