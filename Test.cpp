#include <iostream>
#include <cstdlib>
#include "GA.h"
//#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include "Range_Random.h"
#include "GA_Time.h"
#include <chrono>
#include <windows.h>

////////////
//#include <crtdbg.h>
#include <stdio.h>
#include <string.h>

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
	{
		auto a = opt::createGAGroup(test_Func2, 100);
		a.setBoundary({ {0, 25}, {0, 35} });
		
		a.setMaxGeneration(2000);
		//a.setMaxRuntime(Second(0.6));
		
		a.setThreadNum(2);
		//a.setThreadNum(1);

		// profile
		DWORD start = GetTickCount();
		a.start();
		a.wait_result();
		DWORD end = GetTickCount();
		cout << "The run time is:" << (end - start) / 1000.0 << " s" << endl;
		a.test();
	}

	//decltype(a) b(a);
	
	system("pause");
	return 0;
}
