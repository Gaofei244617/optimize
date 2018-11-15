#pragma once
#include "Individual.h"

namespace opt
{
	struct GA_Info
	{
		double time;                     // ��ǰ����ʱ��(��)
		int NGen;                        // ��ǰ��������
		Individual best_indiv;           // ��ǰ���Ÿ���

		GA_Info(double t, int N, const Individual& indiv)
			:time(t),
			NGen(N),
			best_indiv(indiv)
		{}
	};
}
