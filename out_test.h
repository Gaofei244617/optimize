#pragma once
#include <iostream>
#include "Individual.h"
#include <vector>
#include <fstream>

using namespace std;

void printIndivs(opt::Individual* ptr, const int N,vector<vector<double>> vec)
{
	fstream file;
	file.open("mydata.txt", ios::app);
	file << "Indivs:";
	for (int i = 0; i < N; i++)
	{
		file << (ptr + i)->fitness << "(" << (ptr + i)->vars[0] << "," << (ptr + i)->vars[1] <<"," << vec[i][0] <<","<< vec[i][1]  <<"," << vec[i][2] << ")" << "  ";
	}
	
	file << endl;
	file.close();
}
