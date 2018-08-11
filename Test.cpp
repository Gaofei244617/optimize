#include <iostream>
#include <cstdlib>
#include "GA.h"
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include "Range_Random.h"
#include "GA_Time.h"
#include <chrono>
#include <windows.h>

using namespace std;
using namespace opt;

double fitTest(double x, double y, double z)
{
	return x * x - y * z + z;
}

long double getSec(const Second& s)
{
	return s.value;
}

// 测试用例，最优解(x,y) = (13,27)
double test_Func(double x, double y)
{
	double a = sqrt((x - 13)*(x - 13) + (y - 27)*(y - 27)) + 2;
	double b = cos(a - 2) + 2;
	return b * (50 - a);
}

// 测试用例，最优解(x,y) = (13,17)
double test_Func2(double x, double y)
{
	return -0.2*((x - 13)*(x - 13) + (y - 17)*(y - 17)) + 21;
}


int main()
{	
	//cout << getSec(Second(41.2)) << endl;
	//cout << getSec(Minute(2)) << endl;
	//cout << getSec(Hour(1)) << endl;
	//cout << "Size of GA_State is: " << sizeof(GA_State) << endl;
	//cout << "Size of Individual is: " << sizeof(Individual) << endl;

	auto a = opt::createGAGroup(test_Func2, 10000);
	a.setBoundary({ {0, 25}, {0, 35} });
	a.setMaxGeneration(2000);
	a.setThreadNum(4);

	// profile
	DWORD start = GetTickCount();
	//a.single_start();
	a.start();
	DWORD end = GetTickCount();
	cout << "The run time is:" << (end - start) / 1000.0 << " s" << endl;

	a.test();

	//decltype(a) b(a);
	
	system("pause");
	return 0;
}
