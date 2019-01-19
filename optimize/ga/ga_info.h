#ifndef _GA_INFO_H_
#define _GA_INFO_H_

#include "ga_individual.h"

namespace opt
{
	struct GA_Info
	{
		const double time;                        // 当前迭代时间(秒)
		const int NGen;                           // 当前迭代次数
		const GA_Individual best_indiv;           // 当前最优个体

		GA_Info(double t, int N, const GA_Individual& indiv)
			:time(t),
			NGen(N),
			best_indiv(indiv)
		{}
	};
}

#endif
