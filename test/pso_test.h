#ifndef _PSO_TEST_H_
#define _PSO_TEST_H_

#include <iostream>
#include <vector>
#include <windows.h>
#include "..\opt_time.h"
#include "..\optimize.h"
#include "fit_funcs.h"

using namespace std;
using namespace opt;

template<class T>
void pso_out(const PSO<T>& pso)
{
    auto best = pso.getBestIndivs();
    for (auto& indiv : best)
    {
        cout << indiv.fitness << endl;
    }
}

void pso_test()
{
    auto a = opt::createPSO(test_Func2, 5);

    a.setBoundary({ {0, 25}, {0, 35} });
    a.setMaxGeneration(10);
    //a.setMaxRuntime(Second(0.6));
    a.setThreadNum(4);
    //a.setThreadNum(1);

    a.start();
    a.wait_result();

    pso_out(a);

    // auto b = a.clone();
    // PSO<double(double, double)> c(b);
}

#endif
