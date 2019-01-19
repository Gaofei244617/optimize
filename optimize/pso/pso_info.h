#ifndef _PSO_INFO_H_
#define _PSO_INFO_H_

#include "pso_individual.h"

namespace opt
{
	struct PSO_Info
	{
		const double time;                           // ��ǰ����ʱ��(��)
		const int NGen;                              // ��ǰ��������
		const PSO_Individual best_indiv;             // ��ǰ���Ÿ���

		PSO_Info(double t, int N, const PSO_Individual& indiv)
			:time(t),
			NGen(N),
			best_indiv(indiv)
		{}
	};
}

#endif
