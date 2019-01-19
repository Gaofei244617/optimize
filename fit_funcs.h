#ifndef  _FIT_FUNC_H_
#define _FIT_FUNC_H_

#include <cmath>

double fitTest(double x, double y, double z)
{
    return x * x - y * z + z;
}

// 测试用例1，最优解(x,y) = (13,27)
double test_Func(double x, double y)
{
    double a = sqrt((x - 13)*(x - 13) + (y - 27)*(y - 27)) + 2;
    double b = cos(a - 2) + 2;
    return b * (50 - a);
}

// 测试用例2，最优解(x,y) = (13,17)
double test_Func2(double x, double y)
{
    return -0.2*((x - 13)*(x - 13) + (y - 17)*(y - 17)) + 21;
}

#endif
