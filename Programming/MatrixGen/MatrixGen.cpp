#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <climits>
#include <vector>

using namespace std;

void PrintResults(vector<int> pairs, int nodes);

int main(int argc, char* argv[])
{
	int nodes = argc < 2
		? 100
		: atoi(argv[1]);;

	vector<int> pairs(pow(nodes, 2));
	cout << "number of elements generated: " << pairs.size() << endl;

	int row = 0;
	while (row < nodes)
	{
		int col = row + 1;									// distance to itself is 0
		pairs[row * nodes + row] = 0;
		while (col < nodes)
		{
			int value = (rand() % INT_MAX + 1) % 2 == 0		// a 7th of pairs are unconnected
				? -1
				: ((rand() % INT_MAX + 1) % 100) + 1;		// distance between 1 and 100
			pairs[(col + (row*nodes))] = value;
			pairs[(row + (col*nodes))] = value;
			++col;
		}
		++row;
	}

	PrintResults(pairs, nodes);
	return 0;
}

void PrintResults(vector<int> pairs, int nodes)
{
	ofstream myfile;
	myfile.open("./input.txt", ios::out);

	for (int i = 0; i < pairs.size(); ++i)
	{
		myfile << pairs[i] << "\t";
		if (i % nodes == (nodes - 1))
			myfile << endl;
	}

	myfile.close();
}
