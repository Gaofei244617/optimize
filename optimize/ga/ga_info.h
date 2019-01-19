#ifndef _GA_INFO_H_
#define _GA_INFO_H_

#include "ga_individual.h"

namespace opt
{
	struct GA_Info
	{
		const double time;                        // ��ǰ����ʱ��(��)
		const int NGen;                           // ��ǰ��������
		const GA_Individual best_indiv;           // ��ǰ���Ÿ���

		GA_Info(double t, int N, const GA_Individual& indiv)
			:time(t),
			NGen(N),
			best_indiv(indiv)
		{}
	};
}

#endif
