#ifndef _CROSS_FACTOR_H_
#define _CROSS_FACTOR_H_

#include <utility>
#include <cmath>
#include "Range_Random.h"

namespace opt
{
	// SBX�����㷨(ģ������Ƶ��㽻��)
	// �Ӵ����� (1) p1 + p2 = c1 + c2, (2) beta = |(c2 - c1) / (p2 - p1)|
	// p1��p2Ϊ����, c1��c2Ϊ�Ӵ�, һ�Ը���������������Ӵ�
	// �βηֱ�Ϊ������ĸ����ͬλ�õ�һ�����򣬷��ؽ�����������������
	std::pair<double, double> cross_SBX(const double m, const double f)
	{
		double beta = 0;                              // ���ȷֲ�����
		double nc = 1.0;                              // ����ֲ�ָ��(����0)���Ƽ�Ϊ1; ncԽ���Ӵ������븸��ԽԶ

		double rand_u = random_real(0, 1);            // ��������� 
		if (rand_u <= 0.5)
		{
			beta = std::pow(2.0 * rand_u, 1.0 / (nc + 1.0));
		}
		else
		{
			beta = std::pow(1.0 / (2.0 * (1.0 - rand_u)), 1.0 / (nc + 1));
		}

		// ���򽻲�����Ӵ�����
		double mm = 0.5 * ((1 + beta)*m + (1 - beta)*f);
		double ff = 0.5 * ((1 - beta)*m + (1 + beta)*f);

		return std::make_pair(mm, ff);
	}
}

#endif // !_CROSS_FACTOR_H_

