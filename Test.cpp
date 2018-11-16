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

////////////////////////////////////////////////////////////////////////
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

void monitor(const GA_Info& indiv, void* dat)
{
	cout << "***********************************" << endl;
	cout << "NGen: " << indiv.NGen << endl;
	cout << "Time: " << indiv.time << endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
void out_res(T& a)
{
	//////////////////////////////// 输出计算结果 //////////////////////////////////////////////////////////
	// 输出停止条件
	// Stop Code : -1-未开始迭代; 0-正在迭代; 1-达到最大迭代次数; 2-达到最大迭代时间; 3-最优解收敛于稳定值; 4-人为停止迭代
	switch (a.getStopCode())
	{
	case 1:
		std::cout << "Stop condition: reach the max generation." << std::endl;
		break;
	case 2:
		std::cout << "Stop condition: reach the max time." << std::endl;
		break;
	case 3:
		std::cout << "Stop condition: reach the convergency." << std::endl;
		break;
	case 4:
		std::cout << "Stop condition: killed." << std::endl;
		break;
	}
	cout << endl;

	auto res = a.getBestIndivs();
	cout << "子代最优解进化过程：" << endl;
	cout << res[0].fitness << endl;
	cout << endl;
	for (std::size_t i = 1; i < res.size(); i++)
	{
		if (res[i].fitness > res[i - 1].fitness)
		{
			cout << res[i].fitness << endl;
			cout << endl;
		}
	}
}

int main()
{
	auto a = opt::createGAGroup(test_Func2, 500000);
	a.setBoundary({ {0, 25}, {0, 35} });
	a.setCrossProb(0.95);

	a.setMaxGeneration(6);
	//a.setMaxRuntime(Second(0.6));

	a.setThreadNum(4);
	//a.setThreadNum(1);

	a.setMonitor(monitor, nullptr);

	// profile
	DWORD start = GetTickCount();

	a.start();

	//a.pause();
	//cout << "pause..." << endl;
	//a.proceed();
	//cout << "go on..." << endl;

	//auto b = a.clone();
	//GAGroup<double(double, double)> c(b);

	//a.kill();

	a.wait_result();

	DWORD end = GetTickCount();

	cout << "The run time is:" << (end - start) / 1000.0 << " s" << endl;
	cout << endl;

	out_res(a);

	//////////////////////////////////////////////////////////////////////////////////////////////

	std::cout << std::endl;

	//auto flag = b.start();
	//c.start();

	//cout << "B = " << flag << endl;
	//b.wait_result();
	//out_res(b);

	//c.wait_result();

	//out_res(c);

	system("pause");
	return 0;
}