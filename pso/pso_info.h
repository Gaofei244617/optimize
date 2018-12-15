#ifndef _PSO_INFO_H_
#define _PSO_INFO_H_

#include "pso_individual.h"

namespace opt
{
	struct PSO_Info
	{
		const double time;                           // 当前迭代时间(秒)
		const int NGen;                              // 当前迭代次数
		const PSO_Individual best_indiv;             // 当前最优个体

		PSO_Info(double t, int N, const PSO_Individual& indiv)
			:time(t),
			NGen(N),
			best_indiv(indiv)
		{}
	};
}

#endif
