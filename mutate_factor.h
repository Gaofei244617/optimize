#ifndef _MUTATE_FACTOR_H_
#define _MUTATE_FACTOR_H_

#include <cmath>
#include "Range_Random.h"

namespace opt
{
	// ����ʽ����
	double mutate_PM(const double gene, const double LOW, const double UP)
	{
		double delta = 0;
		double u = random_real(0, 1);
		double nc = 1;   // �ֲ�ָ��

		if (u <= 0.5)
		{
			double delta1 = (gene - LOW) / (UP - LOW);
			delta = std::pow(2 * u + (1 + 2 * u)*std::pow(1 - delta1, nc + 1), 1.0 / (nc + 1)) - 1.0;
		}
		else
		{
			double delta2 = (UP - gene) / (UP - LOW);
			delta = 1.0 - std::pow(2 * (1 - u) + 2 * (u - 0.5)*std::pow(1 - delta2, nc + 1), 1.0 / (nc + 1));
		}

		// �Ӵ�����������
		double gene_child = gene + delta * (UP - LOW);

		// ���Ӵ����򳬳��߽�
		if (gene_child < LOW || gene_child > UP)
		{
			gene_child = random_real(LOW, UP);
		}

		return gene_child;
	}
}

#endif // !_MUTATE_FACTOR_H_

