#pragma once
#include "Individual.h"

namespace opt
{
	struct GA_Info
	{
		double time;                     // 当前迭代时间(秒)
		int NGen;                        // 当前迭代次数
		Individual best_indiv;           // 当前最优个体

		GA_Info(double t, int N, const Individual& indiv)
			:time(t),
			NGen(N),
			best_indiv(indiv)
		{}
	};
}
