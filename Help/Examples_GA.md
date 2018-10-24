```
#include <iostream>
#include <vector>
#include "opt_time.h"

using namespace std;
using namespace opt;  // �Ż��㷨���ڵ������ռ�

// ��������������Ӧ�Ⱥ��������Ž⣺(x,y) = (13,17)
double test_Func2(double x, double y)
{
    return -0.2*((x - 13)*(x - 13) + (y - 17)*(y - 17)) + 21;
}

int main()
{
    // ����һ��һ���Ŵ���Ⱥ,����500������,default = 1000
    auto a = opt::createGAGroup(test_Func2, 200);
    // ��һ�ַ���
    // GAGroup<double(double, double)> a(test_Func2, 200);

    a.setBoundary({ {0, 25}, {0, 35} });   // ����fitness������������
    a.setMaxGeneration(20);                // ��������������

    a.setCrossProb(0.95);                  // ���ý������,default=0.6
    a.setThreadNum(4);                     // ���ò��м�����߳�,default=1

    // ��Ⱥ��ʼ��������
    a.start();

    // ������ǰ�̣߳��ȴ��Ż����
    a.wait_result();

    // ���ֹͣ����
    // Stop Code : -1-δֹͣ; 0-���Ž��������ȶ�ֵ; 1-�ﵽ����������; 2-�ﵽ������ʱ��; 3-��Ϊֹͣ����
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

    // ��ȡ���ε��������Ž�
    vector<opt::Individual> res = a.getBestIndivs();

    // ������Ž���������
    cout << "�Ӵ����Ž�������̣�" << endl;
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