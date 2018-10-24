```
#include <iostream>
#include <vector>
#include "opt_time.h"

using namespace std;
using namespace opt;  // 优化算法所在的命名空间

// 测试用例——适应度函数，最优解：(x,y) = (13,17)
double test_Func2(double x, double y)
{
    return -0.2*((x - 13)*(x - 13) + (y - 17)*(y - 17)) + 21;
}

int main()
{
    // 创建一个一个遗传种群,包含500个个体,default = 1000
    auto a = opt::createGAGroup(test_Func2, 200);
    // 另一种方法
    // GAGroup<double(double, double)> a(test_Func2, 200);

    a.setBoundary({ {0, 25}, {0, 35} });   // 设置fitness函数变量区间
    a.setMaxGeneration(20);                // 设置最大迭代次数

    a.setCrossProb(0.95);                  // 设置交叉概率,default=0.6
    a.setThreadNum(4);                     // 设置并行计算的线程,default=1

    // 种群开始迭代进化
    a.start();

    // 阻塞当前线程，等待优化结果
    a.wait_result();

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

    // 获取历次迭代的最优解
    vector<opt::Individual> res = a.getBestIndivs();

    // 输出最优解收敛过程
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

    system("pause");
    return 0;
}

```