

#include <iostream>
#include <string>
#include <ctime>
#include <chrono>

#include "LoadData.h"
using namespace std;
using namespace std::chrono;


int main()
{

	/*int series_num = 5;

	fstream exp1ExecTime;
	exp1ExecTime.open("Results//Exp1//ExecTime.txt", fstream::out);

	double execTime = 0.0;

	for(int Nf = 15; Nf <= 25; Nf += 10)
	{
		for(int Ni = 100; Ni <= 160; Ni += 10)
		{
			for(K = 20; K <= 100; K += 20)
			{
				for(int series = 1; series <= series_num; series++)
				{
				string path = "Exp1//DataSe" + to_string(series) + "Ps5Pn10Ni" + to_string(Ni) + "Nf" + to_string(Nf) +".txt";
				cout << "Dataset " << path  << " " << K << " ";

				LoadDataset(path);
				TransformData();
				SortDataset();

				high_resolution_clock::time_point t1 = high_resolution_clock::now();
				Miner();
				high_resolution_clock::time_point t2 = high_resolution_clock::now();
				auto duration = duration_cast<seconds>( t2 - t1 ).count();

				execTime += duration;
				cout << duration << endl;

				ClearStructures();
				theta = 1.0;
				}

				exp1ExecTime << execTime/(double)series_num << "\t";
				execTime = 0.0;

			}

			exp1ExecTime << endl;
		}

		exp1ExecTime << endl << "#########" << endl;
	}

	exp1ExecTime.close();*/

	int series_num = 5;

	fstream exp2ExecTime;
	exp2ExecTime.open("Results//Exp2//ExecTime.txt", fstream::out);

	double execTime = 0.0;

	for(int Nf = 15; Nf <= 15; Nf += 10)
	{
		for(int Ni = 10; Ni <= 60; Ni += 10)
		{
			for(K = 20; K <= 100; K += 20)
			{
				for(int series = 1; series <= series_num; series++)
				{
				string path = "Exp2//DataSe" + to_string(series) + "Ps5Pn5Ni" + to_string(Ni) + "Nf" + to_string(Nf) +".txt";
				cout << "Dataset " << path  << " " << K << " ";

				LoadDataset(path);
				TransformData();
				SortDataset();

				high_resolution_clock::time_point t1 = high_resolution_clock::now();
				Miner();
				high_resolution_clock::time_point t2 = high_resolution_clock::now();
				auto duration = duration_cast<seconds>( t2 - t1 ).count();

				execTime += duration;
				cout << duration << endl;

				ClearStructures();
				theta = 1.0;
				}

				exp2ExecTime << execTime/(double)series_num << "\t";
				execTime = 0.0;
			}

			exp2ExecTime << endl;
			}

			exp2ExecTime << endl << "#########" << endl;
		}

		exp2ExecTime.close();
}
