#include <iostream>
#include <cstdlib>
#include "optimize.h"
//#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include "range_random.h"
#include "opt_time.h"
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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{	
	auto a = opt::createGAGroup(test_Func2,5);
	a.setBoundary({ {0, 25}, {0, 35} });
	a.setCrossProb(0.95);
		
	a.setMaxGeneration(20);
	//a.setMaxRuntime(Second(0.6));
		
	a.setThreadNum(4);
	//a.setThreadNum(1);

	// profile
	DWORD start = GetTickCount();
	a.start();
	a.wait_result();
	DWORD end = GetTickCount();
	cout << "The run time is:" << (end - start) / 1000.0 << " s" << endl;
	cout << endl;

	//////////////////////////////// 输出计算结果 //////////////////////////////////////////////////////////	
	// 输出停止条件
	// Stop Code : -1-未停止; 0-最优解收敛于稳定值; 1-达到最大迭代次数; 2-达到最大迭代时间; 3-人为停止迭代
	switch (a.getStopCode())
	{
	case 0:
		std::cout << "Stop condition: reach the convergency." << std::endl;
		break;
	case 1:
		std::cout << "Stop condition: reach the max generation." << std::endl;
		break;
	case 2:
		std::cout << "Stop condition: reach the max time." << std::endl;
		break;
	}
	cout << endl;

	auto res = a.getBestIndivs();
	cout << "子代最优解进化过程：" << endl;
	cout << res[0].fitness << endl;
	cout << endl;
	for (size_t i = 1; i < res.size(); i++)
	{
		if (res[i].fitness > res[i - 1].fitness)
		{
			cout << res[i].fitness << endl;
			cout << endl;
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////

	//decltype(a) b(a);
	
	system("pause");
	return 0;
}
